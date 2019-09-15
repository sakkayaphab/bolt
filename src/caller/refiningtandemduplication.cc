#include "refiningtandemduplication.h"
#include "smithwaterman.h"
#include "readdepthanalysis.h"

RefiningTandemDuplication::RefiningTandemDuplication()
{
    variantresult.setVariantType("DUP");
}

void RefiningTandemDuplication::execute()
{
    variantresult.setChr(evidence.getChr());
    variantresult.setEndChr(evidence.getEndChr());

    prepareBamReader();
    ReadDepthAnalysis rda(filepath);

    first();
    if (variantresult.isQuailtyPass())
    {
        return;
    }
    // second();
    // if (variantresult.isQuailtyPass())
    // {
    //     return;
    // }
}

void RefiningTandemDuplication::first()
{
    std::string findRange = convertRangeToString(evidence.getChr(), evidence.getPos() + evidence.getCiPosLeft(), evidence.getPos() + evidence.getCiPosRight());

    const char *range = findRange.c_str();
    // const char *mChr = evidence.getChr().c_str();
    // std::cout << range << " / " << samplestat->getReadLength() << std::endl;
    refineStartToEnd(range);
}

void RefiningTandemDuplication::refineStartToEnd(const char *range)
{
    hts_itr_t *iter = NULL;

    iter = sam_itr_querys(bam_index, bam_header, range);
    if (iter == NULL)
        return;
    read = bam_init1();
    ReadParser readparser;
    readparser.setBamHeader(bam_header);
    readparser.setBamRead(read);

    int32_t positionStartReference = evidence.getEnd() + evidence.getCiEndLeft();
    int32_t positionEndReference = evidence.getEnd() + evidence.getCiEndRight();
    std::string seqrefString = fastareader.getSeqbyPosition(std::string(evidence.getEndChr()),
                                                            positionStartReference,
                                                            positionEndReference);
    // std::cout << "chr : " << evidence.getEndChr() << ", pos : " << positionStartReference << " end : " << positionEndReference << "ci end left :" << evidence.getCiEndLeft() << " / " << evidence.getCiEndRight() << std::endl;

    replaceSeqToUppercase(&seqrefString);
    const char *seqC = seqrefString.c_str();
    // std::cout << seqrefString << std::endl;

    // Alignment alignment(&seqrefString);
    // alignment.genarateMatrix("DUPSTART");
    // alignment.setPosReference(positionStartReference);
    if (seqrefString.size() < 50)
    {
        return;
    }

    StringSearchAlignment ssa;
    ssa.setReference(seqrefString);
    ssa.setSVType("DUPSTART");
    ssa.setPosReference(positionStartReference);
    ssa.buildReference();

    std::map<std::pair<int32_t, int32_t>, MatchRead> listPosition;
    std::map<int32_t, int> SCReadLists;
    while (sam_itr_next(inFile, iter, read) >= 0)
    {
        if (readparser.isUnmapped())
        {
            continue;
        }

        if (readparser.isNotPassingFilters())
        {
            continue;
        }

        if (readparser.isPCR())
        {
            continue;
        }

        if (readparser.isSupplementaryAlignment())
        {
            continue;
        }

        auto cigar = readparser.getCigar();

        if (cigar.at(0).getOperatorName() != 'S')
        {
            continue;
        }

        if (cigar.at(0).getLength() <= 10)
        {
            continue;
        }

        // if (!(cigar.at(0).getOperatorName() == 'S' && cigar.at(1).getOperatorName() == 'M'))
        // {
        //     continue;
        // }

        // if (readparser.getAlignMD().size()<=1) {
        //     continue;
        // }

        std::string fullRead = readparser.getSequence();

        // std::cout << cigar.at(cigar.size() - 1).getOperatorName() << "--->>>> : " << fullRead.size() << " = " << cigar.at(cigar.size() - 1).getLength() << std::endl;
        // std::cout << "---------- Read ----------" << std::endl;
        // readparser.getAlignMD();
        // std::cout << "pos : " << readparser.getPos() << std::endl;
        // std::cout << "end : " << readparser.getEnd() << std::endl;
        // std::cout << "sc : " << cigar.at(cigar.size() - 1).getLength() << std::endl;

        // std::cout << "seq : " << fullRead << std::endl;
        // std::cout << "---------- GGGG ----------" << cigar.at(cigar.size() - 1).getLength() <<  std::endl;

        // std::vector<SmithWaterman::ScoreAlignment> result = alignment.alignDuplicationTargetAtStart(&fullRead);
        StringSearchConfig ssc;
        ssc.setAllowMissMatch(2);
        ssc.setMaxContinueMissMatch(1);
        std::vector<StringSearch::Score> result = ssa.alignDuplicationTargetAtStart(&fullRead, &ssc);
        for (auto n : result)
        {
            if (n.matchCount <= 4)
            {
                continue;
            }

            int32_t mPos = n.endseq + readparser.getPosOfSeq();
            int32_t mEnd = n.end;

            if (mEnd <= mPos + 2)
            {
                continue;
            }

            // std::cout << "Read Pos : " << readparser.getPos() << std::endl; 
            // std::cout << "mPos : " << mPos << std::endl;
            // std::cout << "mtEnd : " << mEnd << std::endl;
            // std::cout << "seq : " << readparser.getSequence().substr(0,n.endseq) << std::endl;
            
            // std::cout << "pattern : " << n.matchSeqPattern << std::endl;
            // std::cout << "mEnd : " << mEnd << std::endl;
            // std::cout << "---------- MD TAG ----------" << std::endl;

            auto rangeMapping = n.endseq;
            listPosition[std::make_pair(mPos, mEnd)].NumberOfMatchRead++;
            listPosition[std::make_pair(mPos, mEnd)].MatchLists.push_back(rangeMapping);
            listPosition[std::make_pair(mPos, mEnd)].MapQLists.push_back(readparser.getMapQuality());

            // for (int i = 0; i < rangeMapping; i++)
            // {
            //     int32_t tempPos = mPos - i;
            //     int32_t tempEnd = mEnd - i;

            //     if (i > 4)
            //     {
            //         break;
            //     }

            //     listPosition[std::make_pair(tempPos, tempEnd)].NumberOfMatchRead++;
            //     listPosition[std::make_pair(tempPos, tempEnd)].MatchLists.push_back(rangeMapping - i);
            //     listPosition[std::make_pair(tempPos, tempEnd)].MapQLists.push_back(readparser.getMapQuality());
            //     // if (listPosition[std::make_pair(tempPos, tempEnd)].Sequence.size() < n.scorepattern - 1)
            //     // {
            //     //     listPosition[std::make_pair(tempPos, tempEnd)].maxMatchSequence = n.pattern.size() - 1;
            //     // }
            // }
        }
    }

    RefiningTandemDuplication::calculateFinalBreakpoint(&listPosition);

    hts_itr_destroy(iter);
    return;
}

void RefiningTandemDuplication::second()
{
    std::string findRange = convertRangeToString(evidence.getChr(), evidence.getEnd() + evidence.getCiEndLeft(), evidence.getEnd() + evidence.getCiEndRight() + 1000);

    const char *range = findRange.c_str();
    // const char *mChr = evidence.getChr().c_str();
    std::cout << range << " / " << samplestat->getReadLength() << std::endl;
    refineEndToStart(range);
}

void RefiningTandemDuplication::refineEndToStart(const char *range)
{
    hts_itr_t *iter = NULL;

    iter = sam_itr_querys(bam_index, bam_header, range);
    if (iter == NULL)
        return;
    read = bam_init1();
    readparser.setBamHeader(bam_header);
    readparser.setBamRead(read);

    int32_t positionStartReference = evidence.getPos() + evidence.getCiPosLeft();
    int32_t positionEndReference = evidence.getPos() + evidence.getCiPosRight();
    std::string seqrefString = fastareader.getSeqbyPosition(std::string(evidence.getEndChr()),
                                                            positionStartReference,
                                                            positionEndReference);
    replaceSeqToUppercase(&seqrefString);
    // std::cout << "ref pos : " << positionStartReference
    //           << " ref end : " << positionEndReference << std::endl;
    // std::cout << seqrefString << std::endl;

    std::map<std::pair<int32_t, int32_t>, MatchRead> listPosition;
    std::map<int32_t, int> SCReadLists;

    //  Alignment alignment(&seqrefString);
    // alignment.genarateMatrix("DUPEND");
    // alignment.setPosReference(positionStartReference);
    if (seqrefString.size() < 50)
    {
        return;
    }
    StringSearchAlignment ssa;
    ssa.setReference(seqrefString);
    ssa.setSVType("DUPEND");
    ssa.setPosReference(positionStartReference);
    ssa.buildReference();
    // std::cout << positionStartReference << " / " << positionEndReference << std::endl;
    while (sam_itr_next(inFile, iter, read) >= 0)
    {

        if (readparser.isUnmapped())
        {
            continue;
        }

        if (readparser.isNotPassingFilters())
        {
            continue;
        }

        if (readparser.isPCR())
        {
            continue;
        }

        if (readparser.isSupplementaryAlignment())
        {
            continue;
        }

        auto cigar = readparser.getCigar();
        // if (cigar.size() <= 1)
        // {
        //     continue;
        // }
        if (!(cigar.at(cigar.size() - 1).getOperatorName() == 'S'))
        {
            continue;
        }

         if (cigar.at(cigar.size() - 1).getLength() <= 2)
        {
            continue;
        }
        // if (!(cigar.at(cigar.size() - 1).getOperatorName() == 'S'))
        // {
        //     continue;
        // }

        // if (readparser.getAlignMD().size()<=1) {
        //     continue;
        // }

        std::string fullRead = readparser.getSequence();

        // std::vector<SmithWaterman::ScoreAlignment> result = alignment.alignDuplicationTargetAtEnd(&fullRead);
        StringSearchConfig ssc;
        ssc.setAllowMissMatch(2);
        ssc.setMaxContinueMissMatch(1);
        std::vector<StringSearch::Score> result = ssa.alignDuplicationTargetAtEnd(&fullRead, &ssc);
        for (auto n : result)
        {
            if (n.matchCount <= 4)
            {
                continue;
            }

            int32_t mPos = n.pos;
            int32_t mEnd = n.posseq + readparser.getPosOfSeq();

            if (mPos >= mEnd + 2)
            {
                continue;
            }

            // std::cout << "mPos : " << mPos << std::endl;
            // std::cout << "mtEnd : " << mEnd << std::endl;
            // std::cout << "pattern : " << n.pattern << std::endl;
            // std::cout << "mEnd : " << mEnd << std::endl;
            // std::cout << "---------- MD TAG ----------" << std::endl;

            auto rangeMapping = n.endseq - n.posseq;
            listPosition[std::make_pair(mPos, mEnd)].NumberOfMatchRead++;
            listPosition[std::make_pair(mPos, mEnd)].MatchLists.push_back(rangeMapping);
            listPosition[std::make_pair(mPos, mEnd)].MapQLists.push_back(readparser.getMapQuality());

            // for (int i = 0; i < rangeMapping; i++)
            // {
            //     int32_t tempPos = mPos + i;
            //     int32_t tempEnd = mEnd + i;

            //     if (i < 10)
            //     {
            //         break;
            //     }

            //     listPosition[std::make_pair(tempPos, tempEnd)].NumberOfMatchRead++;
            //     listPosition[std::make_pair(tempPos, tempEnd)].MatchLists.push_back(rangeMapping - i);
            //     listPosition[std::make_pair(tempPos, tempEnd)].MapQLists.push_back(readparser.getMapQuality());
            //     // if (listPosition[std::make_pair(tempPos, tempEnd)].Sequence.size() < n.scorepattern - 1)
            //     // {
            //     //     listPosition[std::make_pair(tempPos, tempEnd)].maxMatchSequence = n.pattern.size() - 1;
            //     // }
            // }
        }
    }

    RefiningTandemDuplication::calculateFinalBreakpoint(&listPosition);
    // std::cout << variantresult.getPos() << std::endl;
    hts_itr_destroy(iter);
    return;
}

void RefiningTandemDuplication::calculateFinalBreakpoint(std::map<std::pair<int32_t, int32_t>, RefiningSV::MatchRead> *listPosition)
{
    int32_t bPos = 0;
    int32_t bEnd = 0;
    int32_t bHit = 0;
    uint8_t bMaxQuality = 0;
    int bMaxMatchSize = 0;
    int bFrequency = 0;
    std::vector<uint8_t> bMapQList;
    int32_t svlength = evidence.getEndDiscordantRead() - evidence.getPosDiscordantRead() - samplestat->getMedianSampleStat();
    // std::cout << "svlength :" << svlength << std::endl;
    int lastscore = 0;

    for (auto const &x : *listPosition)
    {
        int maxMatchSize = getMaxIntFromVector(x.second.MatchLists);

        uint8_t maxQuality = getMaxUInt8FromVector(x.second.MapQLists);

        if (maxMatchSize == 0)
        {
            continue;
        }

        // if (maxMatchSize > 80)
        // {
        //     continue;
        // }

        // if (maxQuality == 0)
        // {
        //     continue;
        // }

        int number = x.second.NumberOfMatchRead;

        int score = (number) * (2 * maxMatchSize);

        if (score > lastscore)
        {
            lastscore = score;
            bPos = x.first.first;
            bEnd = x.first.second;
            bHit = number;
            bMaxMatchSize = maxMatchSize;
            bMapQList = x.second.MapQLists;
        }
    }

    variantresult.setPos(bPos);
    variantresult.setEnd(bEnd);
    variantresult.setFrequency(bHit);
    variantresult.setMapQList(bMapQList);
    variantresult.setChr(evidence.getChr());
    variantresult.setRPMapQ(*evidence.getMapQVector());
    variantresult.setEndChr(evidence.getEndChr());
    variantresult.LNGMATCH = bMaxMatchSize;

    if (bPos == 0)
    {
        return;
    }
    if (bEnd == 0)
    {
        return;
    }

    if (bEnd - bPos > 1000000)
    {
        return;
    }

    if (bEnd - bPos <= 0)
    {
        return;
    }

    variantresult.setQuailtyPass(true);
}