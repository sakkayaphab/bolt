#include "refiningdeletion.h"
#include <fasta/fastareader.h>
#include <string.h>
#include "smithwaterman.h"
#include "alignment.h"
#include "readdepthanalysis.h"

RefiningDeletion::RefiningDeletion()
{
    variantresult.setVariantType("DEL");
}

RefiningDeletion::~RefiningDeletion()
{

}

void RefiningDeletion::execute()
{
    ReadDepthAnalysis rda(filepath);

    variantresult.setChr(evidence.getChr());
    variantresult.setEndChr(evidence.getEndChr());

    // std::cout << "evidence : " << evidence.getPos() << " - " << evidence.getEnd() << std::endl;

    prepareBamReader();
    first();
    resultFromStart = true;
    if (variantresult.isQuailtyPass())
    {
        return;
    }
    resultFromStart = false;
    second();
    if (variantresult.isQuailtyPass())
    {
        return;
    }

    // if (!variantresult.isQuailtyPass()) {
    //     approximate();
    // }
}

void RefiningDeletion::approximate()
{
    variantresult.setPos(evidence.getLastPosDiscordantRead());
    variantresult.setEnd(evidence.getEndDiscordantRead());
    variantresult.setFrequency(evidence.getMapQVector()->size());
    variantresult.setChr(evidence.getChr());
    variantresult.setEndChr(evidence.getEndChr());
    if (variantresult.getEnd() - variantresult.getPos() > 10000)
    {
        variantresult.setQuailtyPass(true);
    }
}

void RefiningDeletion::first()
{
    std::string findRange = convertRangeToString(evidence.getChr(), evidence.getPos() + evidence.getCiPosLeft(), evidence.getPos() + evidence.getCiPosRight());

    const char *range = findRange.c_str();
    // const char *mChr = evidence.getChr().c_str();
    //         std::cout << range  << std::endl;
    refineStartToEnd(range);
}

void RefiningDeletion::second()
{
    std::string findRange = convertRangeToString(evidence.getEndChr(), evidence.getEnd() + evidence.getCiEndLeft(),
                                                 evidence.getEnd() + evidence.getCiEndRight());
    const char *range = findRange.c_str();

    // std::cout << range << "/" << samplestat->getReadLength() << std::endl;
    refineEndToStart(range);
}

void RefiningDeletion::refineStartToEnd(const char *range)
{
    hts_itr_t *iter = NULL;
    // std::cout << "refineStartToEnd range :" << range << std::endl;

    iter = sam_itr_querys(bam_index, bam_header, range);
    if (iter == NULL)
        return;
    read = bam_init1();
    readparser.setBamHeader(bam_header);
    readparser.setBamRead(read);

        // std::cout << "refineStartToEnd NULL" << std::endl;


    int32_t positionStartReference = evidence.getEnd() + evidence.getCiEndLeft();
    std::string seqrefString = fastareader.getSeqbyPosition(std::string(evidence.getEndChr()),
                                                            positionStartReference,
                                                            evidence.getEnd() + evidence.getCiEndRight());
    replaceSeqToUppercase(&seqrefString);
    // const char *seqC = seqrefString.c_str();
    // std::cout << "ref pos : " << positionStartReference
    //           << " ref end : " << evidence.getEnd() + evidence.getCiEndRight() << std::endl;
    // std::cout << seqrefString << std::endl;

    std::map<std::pair<int32_t, int32_t>, MatchRead> listPosition;
    std::map<int32_t, int> SCReadLists;

    // Alignment alignment(&seqrefString);
    // alignment.genarateMatrix("DELSTART");

    // alignment.setPosReference(positionStartReference);
    if (seqrefString.size() < 50)
    {
        std::cout << "back" << std::endl;
        return;
    }

    StringSearchAlignment ssa;
    ssa.setReference(seqrefString);
    ssa.setSVType("DELSTART");
    ssa.setPosReference(positionStartReference);
    ssa.buildReference();
    std::vector<ReadParser::Cigar> cigar;

    // std::cout << "positionStartReference : " << positionStartReference << std::endl;

    bool SCRead;
    int32_t SCsize = 0;
    int32_t AlterSCsize = 0;

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
        // if (readparser.getMapQuality()<15) {
        //             continue;
        //         }
        std::string fullRead = readparser.getSequence();
        AlterSCsize = 0;
        // optimize read
        std::vector<StringSearch::Score> result;
        cigar = readparser.getCigar();
        // std::cout << cigar.size() << std::endl;

        if (cigar.size() <= 1)
        {
            continue;

            // if (haveIndel(cigar)) {
                
            // }else {
            //     continue;
            // }


            SCsize = 0;
            // if (evidence.getSvLength()<250) {
            //     continue;
            // }

            if (readparser.getAlignMD().size() <= 1)
            {
                continue;
            }

            StringSearchConfig ssc;
            ssc.setAllowMissMatch(0);
            ssc.setMaxContinueMissMatch(0);
            result = ssa.alignDeletionTargetAtStart(&fullRead, &ssc);

            // continue;
            SCRead = false;
        }
        else
        {
            if (!(cigar.at(cigar.size() - 1).getOperatorName() == 'S'))
            {
                continue;
            }

            SCsize = cigar.at(cigar.size() - 1).getLength();

            if (SCsize <= 2)
            {
                continue;
            }

            if (cigar.at(0).getOperatorName() == 'S')
            {
                AlterSCsize = cigar.at(0).getLength();
            }

            StringSearchConfig ssc;
            ssc.setAllowMissMatch(2);
            ssc.setMaxContinueMissMatch(1);
            ssc.setMaxAllowAlign(SCsize);
            result = ssa.alignDeletionTargetAtStart(&fullRead, &ssc);

            SCRead = true;
        }

        // if (readparser.getPos()>230471084) {
        //     continue;
        // }

        //  std::cout << "pos : " << readparser.getPos() << std::endl;
        //         continue;

        // if (readparser.getEnd()!=56129321) {
        //     continue;
        // }
        //      std::cout << "---------- New Read ----------" << std::endl;
        //     std::cout << readparser.getPos() << " = " << readparser.getEnd() << " sc:" << cigar.at(cigar.size() - 1).getLength() << std::endl;
        // std::cout << fullRead << std::endl;

        for (auto n : result)
        {

            // continue;
            if (n.matchCount <= 4)
            {
                continue;
            }

            if (SCRead)
            {
                // if (n.matchCount < cigar.at(0).getLength())
                // {
                //     continue;
                // }
            }
            else
            {
                // continue;
                // if (n.missmatchCount >= 1)
                // {
                //     continue;
                // }
                if (readparser.getLastToStartMissMatchPosMD() == 0)
                {
                    continue;
                }

                if (!(readparser.getLastToStartMissMatchPosMD() <= n.matchCount))
                {
                    continue;
                }

                if (n.matchCount<20) {
                    continue;
                }

                if (n.missmatchCount >= 1)
                {
                    continue;
                }

            }

            // if (n.matchCount <= 20) {
            //     if (n.missmatchCount >= 1 ) {
            //         continue;
            //     }
            // }

            int32_t mPos = n.posseq + readparser.getPosOfSeq();
            int32_t mEnd = n.pos;

            if (mEnd - mPos <= 10)
            {
                continue;
            }

            // if (mEnd <= mPos + 2)
            // {
            //     continue;
            // }

            // std::cout << "mPos : " << mPos << std::endl;
            // std::cout << "mtEnd : " << mEnd<< std::endl;
            //  std::cout << "pattern : " << n. << std::endl;
            // std::cout << "mEnd : " << mEnd << std::endl;
            // std::cout << "---------- MD TAG ----------" << std::endl;
            auto rangeMapping = n.endseq - n.posseq;

            listPosition[std::make_pair(mPos, mEnd)].NumberOfMatchRead++;
            listPosition[std::make_pair(mPos, mEnd)].MatchLists.push_back(rangeMapping);
            listPosition[std::make_pair(mPos, mEnd)].MapQLists.push_back(readparser.getMapQuality());
            if (listPosition[std::make_pair(mPos, mEnd)].maxSC < SCsize)
            {
                listPosition[std::make_pair(mPos, mEnd)].maxSC = SCsize;
            }
            if (listPosition[std::make_pair(mPos, mEnd)].maxAlterSC < AlterSCsize)
            {
                listPosition[std::make_pair(mPos, mEnd)].maxAlterSC = AlterSCsize;
            }

            // for (int i = 0; i < rangeMapping; i++)
            // {

            //     int32_t tempPos = mPos + i;
            //     int32_t tempEnd = mEnd + i;

            //     if (i > 6)
            //     {
            //         break;
            //     }
            //     if (SCRead) {
            //         listPosition[std::make_pair(tempPos, tempEnd)].alignWithSoftClipped = true;
            //     }
            //     listPosition[std::make_pair(tempPos, tempEnd)].NumberOfMatchRead++;
            //     listPosition[std::make_pair(tempPos, tempEnd)].MatchLists.push_back(rangeMapping + i);
            //     listPosition[std::make_pair(tempPos, tempEnd)].MapQLists.push_back(readparser.getMapQuality());
            //     // if (listPosition[std::make_pair(tempPos, tempEnd)].Sequence.size() < n.scorepattern - 1)
            //     // {
            //     //     listPosition[std::make_pair(tempPos, tempEnd)].maxMatchSequence = n.pattern.size() - 1;
            //     // }
            //     // break;
            // }
        }
    }

    RefiningDeletion::calculateFinalBreakpoint(&listPosition);

    hts_itr_destroy(iter);
    return;
}

void RefiningDeletion::refineEndToStart(const char *range)
{
    hts_itr_t *iter = NULL;

    iter = sam_itr_querys(bam_index, bam_header, range);
    if (iter == NULL)
        return;
    read = bam_init1();
    readparser.setBamHeader(bam_header);
    readparser.setBamRead(read);

    int32_t positionStartReference = evidence.getPos() + evidence.getCiPosLeft();
    std::string seqrefString = fastareader.getSeqbyPosition(std::string(evidence.getEndChr()),
                                                            positionStartReference,
                                                            evidence.getPos() + evidence.getCiPosRight());
    replaceSeqToUppercase(&seqrefString);
    // std::cout << "ref pos : " << positionStartReference
    //           << " ref end : " << evidence.getEnd() + evidence.getCiEndRight() << std::endl;
    // std::cout << seqrefString << std::endl;

    std::map<std::pair<int32_t, int32_t>, MatchRead> listPosition;
    std::map<int32_t, int> SCReadLists;

    // Alignment alignment(&seqrefString);
    // alignment.genarateMatrix("DELEND");

    // alignment.setPosReference(positionStartReference);
    if (seqrefString.size() < 50)
    {
        return;
    }
    StringSearchAlignment ssa;
    ssa.setReference(seqrefString);
    ssa.setSVType("DELEND");
    ssa.setPosReference(positionStartReference);
    ssa.buildReference();
    std::vector<ReadParser::Cigar> cigar;
    int32_t svlength = evidence.getEndDiscordantRead() - evidence.getPosDiscordantRead() - samplestat->getMedianSampleStat();
    bool SCRead = false;
    int32_t SCsize = 0;
    int32_t AlterSCsize = 0;

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

        // if (readparser.getMapQuality()<15) {
        //     continue;
        // }
        AlterSCsize = 0;
        std::string fullRead = readparser.getSequence();
        std::vector<StringSearch::Score> result;

        cigar = readparser.getCigar();
        if (cigar.size() <= 1)
        {
            continue;
            // if (haveIndel(cigar)) {
                
            // }else {
            //     continue;
                
            // }

            SCsize = 0;

            //  if (evidence.getSvLength()<250) {
            //     continue;
            // }

            if (readparser.getAlignMD().size() <= 1)
            {
                continue;
            }
            // continue;
            SCRead = false;

            StringSearchConfig ssc;
            ssc.setAllowMissMatch(0);
            ssc.setMaxContinueMissMatch(0);

            result = ssa.alignDeletionTargetAtEnd(&fullRead, &ssc);
        }
        else
        {
            if (!(cigar.at(0).getOperatorName() == 'S'))
            {
                continue;
            }
            SCsize = cigar.at(0).getLength();

            if (SCsize <= 2)
            {
                continue;
            }

            if (cigar.at(cigar.size() - 1).getOperatorName() == 'S')
            {
                AlterSCsize = cigar.at(cigar.size() - 1).getLength();
            }

            SCRead = true;

            StringSearchConfig ssc;
            ssc.setAllowMissMatch(2);
            ssc.setMaxContinueMissMatch(1);
            ssc.setMaxAllowAlign(SCsize);

            result = ssa.alignDeletionTargetAtEnd(&fullRead, &ssc);
        }

        for (auto n : result)
        {
            // continue;
            if (n.matchCount <= 4)
            {
                continue;
            }

            if (SCRead)
            {
                // if (n.matchCount < cigar.at(0).getLength())
                // {
                //     continue;
                // }
            }
            else
            {
                // continue;

                if (readparser.getStartToEndMissMatchPosMD()==0) {
                    continue;
                }

                if (!(readparser.getStartToEndMissMatchPosMD() <= n.matchCount))
                {
                    continue;
                }

                if (n.matchCount<20) {
                    continue;
                }

                if (n.missmatchCount >= 1)
                {
                    continue;
                }
                // continue;
            }

            // if (n.matchCount <= 20) {
            //     if (n.missmatchCount >= 1 ) {
            //         continue;
            //     }
            // }

            // if (evidence.getSvLength() < 300)
            // {
            //     if (n.missmatchCount >= 1)
            //     {
            //         continue;
            //     }
            // }

            int32_t mPos = n.end;
            int32_t mEnd = readparser.getPosOfSeq() + n.endseq;

            if (mEnd - mPos <= 10)
            {
                continue;
            }

            // if (mPos >= mEnd + 2)
            // {
            //     continue;
            // }
                //  std::cout << "> n.pos : " << mPos << " n.end : " << mEnd << " n.matchCount" << n.matchCount << std::endl;
            //  std::cout << "mPos : " << mPos << std::endl;
            //  std::cout << "mEnd : " << mEnd << std::endl;
            //      std::cout << "pattern : " << n.pattern << std::endl;
            // std::cout << "---------- MD TAG ----------" << std::endl;
            listPosition[std::make_pair(mPos, mEnd)].NumberOfMatchRead++;
            listPosition[std::make_pair(mPos, mEnd)].MatchLists.push_back(n.endseq);
            listPosition[std::make_pair(mPos, mEnd)].MapQLists.push_back(readparser.getMapQuality());
            if (listPosition[std::make_pair(mPos, mEnd)].maxSC < SCsize)
            {
                listPosition[std::make_pair(mPos, mEnd)].maxSC = SCsize;
            }

            if (listPosition[std::make_pair(mPos, mEnd)].maxAlterSC < AlterSCsize)
            {
                listPosition[std::make_pair(mPos, mEnd)].maxAlterSC = AlterSCsize;
            }

            // for (int i = 0; i < n.endseq; i++)
            // {

            //     int32_t tempPos = mPos - i;
            //     int32_t tempEnd = mEnd - i;

            //     if (i > 6)
            //     {
            //         break;
            //     }

            //     if (SCRead) {
            //         listPosition[std::make_pair(tempPos, tempEnd)].alignWithSoftClipped = true;
            //     }

            //     listPosition[std::make_pair(tempPos, tempEnd)].NumberOfMatchRead++;
            //     listPosition[std::make_pair(tempPos, tempEnd)].MatchLists.push_back(n.endseq - i);
            //     listPosition[std::make_pair(tempPos, tempEnd)].MapQLists.push_back(readparser.getMapQuality());
            //     // if (listPosition[std::make_pair(tempPos, tempEnd)].Sequence.size() < n.scorepattern - 1)
            //     // {
            //     //     listPosition[std::make_pair(tempPos, tempEnd)].maxMatchSequence = n.pattern.size() - 1;
            //     // }
            //     // break;
            // }
        }
    }

    RefiningDeletion::calculateFinalBreakpoint(&listPosition);
    hts_itr_destroy(iter);
    return;
}

void RefiningDeletion::calculateFinalBreakpoint(std::map<std::pair<int32_t, int32_t>, RefiningSV::MatchRead> *listPosition)
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

        if (maxMatchSize < 24)
        {
            continue;
        }

        // std::cout << x.first.first << " - " << x.first.second << std::endl;

        // if (maxMatchSize > 80)
        // {
        //     continue;
        // }

        if (maxQuality == 0)
        {
            continue;
        }

        // if (x.second.maxAlterSC >= x.second.maxSC && x.second.alignWithSoftClipped)
        // {
        //     continue;
        // }

        bFrequency = x.second.NumberOfMatchRead;

        if (evidence.getSvLength() < 500)
        {
            
            if (bFrequency <= 1)
            {
                continue;
            }
        }
        else
        {

            if (bFrequency <= 1)
            {
                continue;
            }
        }

        // if (isMatchRef(evidence.getChr(), x.first.first - 1, x.first.first + 20 - 1, evidence.getEndChr(), x.first.second, x.first.second + 20))
        // {
        //     continue;
        // }

        // // confirm
        // if (isMatchRef(evidence.getChr(), x.first.first - 1, x.first.first + 20 - 1, evidence.getEndChr(), x.first.second + 1, x.first.second + 20 + 1))
        // {
        //     continue;
        // }

        // if (isMatchRef(evidence.getChr(), x.first.first - 1 + 2, x.first.first + 20 - 1 + 2, evidence.getEndChr(), x.first.second + 2, x.first.second + 20 + 2))
        // {
        //     continue;
        // }

        // if (isMatchRef(evidence.getChr(), x.first.first + 1 - 20, x.first.first + 1, evidence.getEndChr(), x.first.second - 20, x.first.second))
        // {
        //     continue;
        // }

        // if (isMatchRef(evidence.getChr(), x.first.first + 2 - 20, x.first.first + 1 + 2, evidence.getEndChr(), x.first.second + 20 + 2, x.first.second + 2))
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

    // int CountZeroMapQ = getNumberMapQ(bMapQList,0,0);
    // int CountNotZeroMapQ = getNumberMapQ(bMapQList,1,255);

    // if (CountNotZeroMapQ>=CountZeroMapQ) {
    //     return;
    // }

    


    variantresult.setPos(bPos);
    variantresult.setEnd(bEnd);
    variantresult.setFrequency(bHit);
    variantresult.setRPMapQ(*evidence.getMapQVector());

    variantresult.setMapQList(bMapQList);
    variantresult.setChr(evidence.getChr());
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

    if (bEnd - bPos < 20)
    {
        return;
    }

    variantresult.setQuailtyPass(true);
}

int RefiningDeletion::getNumberMapQ(std::vector<uint8_t> mapqlist,uint8_t start,uint8_t end) {
    int count =0;
    for (auto n:mapqlist) {
        if (n>=start && n<= end) {
            count++;
        }
    }
    return count;
}