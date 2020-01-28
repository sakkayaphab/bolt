#include "refinetranslocation.h"
#include "smithwaterman.h"

RefineTranslocation::RefineTranslocation()
{
    variantresult.setVariantType("BND");
}

void RefineTranslocation::execute()
{
    variantresult.setChr(evidence.getChr());
    variantresult.setEndChr(evidence.getEndChr());

    prepareBamReader();

    // // confirmAreaBySCAtStart("chr1", 6951929 - 10, 6951929);
    // // confirmAreaBySCAtEnd("chr1",6941979,6941979+10);
    first();
    if (variantresult.isQuailtyPass())
    {
        return;
    }

    // std::cout << "---+ run complete +--- : " << std::endl;
}

void RefineTranslocation::first()
{
    std::string findRange = convertRangeToString(evidence.getChr(), evidence.getPos() + evidence.getCiPosLeft(), evidence.getPos() + evidence.getCiPosRight());

    const char *range = findRange.c_str();
    // const char *mChr = evidence.getChr().c_str();
    //         std::cout << range  << std::endl;
    refineStartToEnd(range);
}

void RefineTranslocation::refineStartToEnd(const char *range)
{
    hts_itr_t *iter = NULL;

    iter = sam_itr_querys(bam_index, bam_header, range);
    if (iter == NULL)
        return;
    read = bam_init1();
    readparser.setBamHeader(bam_header);
    readparser.setBamRead(read);

    int32_t positionStartReference = evidence.getEnd() + evidence.getCiEndLeft();
    int32_t positionEndReference = evidence.getEnd() + evidence.getCiEndRight();
    std::string seqrefString = fastareader.getSeqbyPosition(std::string(evidence.getEndChr()),
                                                            positionStartReference,
                                                            positionEndReference);
    replaceSeqToUppercase(&seqrefString);
    // const char *seqC = seqrefString.c_str();
    // std::cout << "ref pos : " << positionStartReference
    //           << " ref end : " << evidence.getEnd() + evidence.getCiEndRight() << std::endl;
    // std::cout << seqrefString << std::endl;

    std::map<std::pair<int32_t, int32_t>, MatchRead> listPosition;
    std::map<int32_t, int> SCReadLists;

    StringSearchAlignment ssa;
    ssa.setReference(seqrefString);
    ssa.setSVType("TRASTART");
    ssa.setPosReference(positionStartReference);
    ssa.buildReference();

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
        // if (!cigar.size() == 2)
        // {
        //     // continue;
        // }

        // if ((cigar.at(0).getOperatorName() == 'S' && cigar.at(1).getOperatorName() == 'M'))
        if (cigar.at(0).getOperatorName() == 'S')
        {
            // continue;
            std::string fullRead = readparser.getSequence();

            // std::vector<SmithWaterman::ScoreAlignment> result = alignment.alignDuplicationTargetAtStart(&fullRead);
            StringSearchConfig ssc;
            std::vector<StringSearch::Score> result = ssa.alignDuplicationTargetAtStart(&fullRead, &ssc);
            for (auto n : result)
            {

                int32_t mPos = n.endseq + readparser.getPosOfSeq();
                int32_t mEnd = n.end;

                // std::cout << "mPos : " << mPos << std::endl;
                // std::cout << "mtEnd : " << mEnd << std::endl;
                // std::cout << "pattern : " << n.pattern << std::endl;
                auto rangeMapping = n.endseq;

                listPosition[std::make_pair(mPos, mEnd)].NumberOfMatchRead++;
                    listPosition[std::make_pair(mPos, mEnd)].MatchLists.push_back(rangeMapping);
                    listPosition[std::make_pair(mPos, mEnd)].MapQLists.push_back(readparser.getMapQuality());

                // for (int i = 0; i < rangeMapping; i++)
                // {
                //     int32_t tempPos = mPos - i;
                //     int32_t tempEnd = mEnd - i;

                //     if (i > rangeMapping - 6)
                //     {
                //         break;
                //     }

                //     listPosition[std::make_pair(tempPos, tempEnd)].NumberOfMatchRead++;
                //     listPosition[std::make_pair(tempPos, tempEnd)].MatchLists.push_back(rangeMapping - i);
                //     listPosition[std::make_pair(tempPos, tempEnd)].MapQLists.push_back(readparser.getMapQuality());
                // }
            }
        }
        else if ((cigar.at(cigar.size() - 1).getOperatorName() == 'S'))
        //  if (true)
        {
            // continue;
            std::string fullRead = readparser.getSequence();
            // std::cout << "-------------" << std::endl;
            // std::cout << " sc : " << cigar.at(cigar.size() - 1).getLength() << std::endl;
            // std::vector<SmithWaterman::ScoreAlignment> result = alignment.alignTranslocationTargetAtStartSCE(&fullRead);
            StringSearchConfig ssc;
            std::vector<StringSearch::Score> result = ssa.alignTranslocationTargetAtStartSCE(&fullRead, &ssc);
            // result = alignment.alignTranslocationTargetAtStartSCE(&fullRead);
            for (auto n : result)
            {

                int32_t mPos = n.posseq + readparser.getPosOfSeq();
                int32_t mEnd = n.pos;

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

                //     if (i > rangeMapping - 6)
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
    }

    RefineTranslocation::calculateFinalBreakpoint(&listPosition);

    hts_itr_destroy(iter);
    return;
}


void RefineTranslocation::calculateFinalBreakpoint(std::map<std::pair<int32_t, int32_t>, RefineSV::MatchRead> *listPosition)
{

    int32_t bPos = 0;
    int32_t bEnd = 0;
    int32_t bHit = 0;
    uint8_t bMaxQuality = 0;
    int bMaxMatchSize = 0;
    int bFrequency = 0;
    std::vector<uint8_t> bMapQList;
    int32_t svlength = evidence.getEndDiscordantRead() - evidence.getPosDiscordantRead() - samplestat->getAverageSampleStat();
    // std::cout << "svlength :" << svlength << std::endl;
    int lastscore = 0;

    for (auto const &x : *listPosition)
    {
        int maxMatchSize = getMaxIntFromVector(x.second.MatchLists);

        uint8_t maxQuality = getMaxUInt8FromVector(x.second.MapQLists);

        if (maxMatchSize < 25)
        {
            continue;
        }
        // std::cout << "maxMatchSize : " << maxMatchSize << std::endl;
        if (maxMatchSize > 80)
        {
            continue;
        }

        if (maxQuality == 0)
        {
            continue;
        }

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
    variantresult.setEndChr(evidence.getEndChr());
    variantresult.setRPMapQ(*evidence.getMapQVector());
    variantresult.LNGMATCH = bMaxMatchSize;

    if (bPos == 0)
    {
        return;
    }

    if (bEnd == 0)
    {
        return;
    }

     if (evidence.getMark()=="SR") {
        variantresult.setMark("SR");
    }

    variantresult.setQuailtyPass(true);
}