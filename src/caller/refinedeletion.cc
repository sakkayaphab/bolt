#include "refinedeletion.h"
#include <fasta/fastareader.h>
#include <string.h>
#include "smithwaterman.h"
#include "alignment.h"
#include "readdepthanalysis.h"

RefineDeletion::RefineDeletion()
{
    variantresult.setVariantType("DEL");
}

RefineDeletion::~RefineDeletion()
{
}

void RefineDeletion::execute()
{
    // std::cout << evidence.getResultVcfFormatString() << std::endl;

    ReadDepthAnalysis rda(filepath);

    variantresult.setChr(evidence.getChr());
    variantresult.setEndChr(evidence.getEndChr());

    if (evidence.getMark() == "SDEL")
    {
        variantresult = evidence;
        variantresult.setQuailtyPass(true);
        return;
    }

    prepareBamReader();

    first();
    second();

    variantresult = getBestResult(resultFirst, resultSecond);
    // std::cout << variantresult.getResultVcfFormatString() << std::endl;
}

Evidence RefineDeletion::getBestResult(Evidence r1, Evidence r2)
{

    if (r1.isQuailtyPass() == false && r2.isQuailtyPass() == false)
    {
        Evidence result;
        return result;
    }

    if (r1.getFrequency() > r2.getFrequency())
    {
        return r1;
    }

    return r2;
}

void RefineDeletion::approximate()
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

void RefineDeletion::first()
{
    std::string findRange = convertRangeToString(evidence.getChr(), evidence.getPos() + evidence.getCiPosLeft(), evidence.getPos() + evidence.getCiPosRight());

    const char *range = findRange.c_str();
    // const char *mChr = evidence.getChr().c_str();
    // std::cout << range  << std::endl;
    refineStartToEnd(range);
}

void RefineDeletion::second()
{
    std::string findRange = convertRangeToString(evidence.getEndChr(), evidence.getEnd() + evidence.getCiEndLeft(),
                                                 evidence.getEnd() + evidence.getCiEndRight());
    const char *range = findRange.c_str();

    // std::cout << range << "/" << samplestat->getReadLength() << std::endl;
    refineEndToStart(range);
}

void RefineDeletion::refineStartToEnd(const char *range)
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

    // Alignment alignment(&seqrefString);
    // alignment.genarateMatrix("DELSTART");

    // alignment.setPosReference(positionStartReference);
    if (seqrefString.size() < 50)
    {
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

            if (haveIndel(cigar))
            {
            }
            else
            {
                continue;
            }

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
            if (n.matchCount <= 3)
            {
                continue;
            }

            if (SCRead)
            {
                if (n.matchCount + n.missmatchCount + 4 < cigar.at(cigar.size() - 1).getLength())
                {
                    continue;
                }
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

                if (n.matchCount < 20)
                {
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

            std::vector<ReadParser::SATag> satag = readparser.getSATag();
            for (ReadParser::SATag sa : satag)
            {
                AlternativeSA tempAltSA;
                tempAltSA.chr = readparser.getChromosomeNameString();
                tempAltSA.pos = sa.pos;
                listPosition[std::make_pair(mPos, mEnd)].AltSA.push_back(tempAltSA);
            }

            // std::cout << "mPos : " << mPos << std::endl;
            // std::cout << "mtEnd : " << mEnd<< std::endl;
            // std::cout << mPos << " " << mEnd << " = " << SCsize << " > " << readparser.getEnd() << std::endl;
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
        }
    }

    resultFirst = RefineDeletion::calculateFinalBreakpoint(&listPosition);
    resultFirst.setEvidenceFrom(resultFirst.getPos());

    hts_itr_destroy(iter);
    return;
}

void RefineDeletion::refineEndToStart(const char *range)
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
    int32_t svlength = evidence.getEndDiscordantRead() - evidence.getPosDiscordantRead() - samplestat->getAverageSampleStat();
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
            if (n.matchCount <= 3)
            {
                continue;
            }

            if (SCRead)
            {
                if (n.matchCount + n.missmatchCount + 4 < cigar.at(0).getLength())
                {
                    continue;
                }
            }
            else
            {
                // continue;

                if (readparser.getStartToEndMissMatchPosMD() == 0)
                {
                    continue;
                }

                if (!(readparser.getStartToEndMissMatchPosMD() <= n.matchCount))
                {
                    continue;
                }

                if (n.matchCount < 15)
                {
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

            // std::cout << mPos << " " << mEnd << " = " << SCsize << std::endl;

            std::vector<ReadParser::SATag> satag = readparser.getSATag();
            for (ReadParser::SATag sa : satag)
            {
                AlternativeSA tempAltSA;
                tempAltSA.chr = readparser.getChromosomeNameString();
                tempAltSA.pos = sa.pos;
                listPosition[std::make_pair(mPos, mEnd)].AltSA.push_back(tempAltSA);
            }

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
        }
    }

    resultSecond = RefineDeletion::calculateFinalBreakpoint(&listPosition);
    resultSecond.setEvidenceFrom(resultSecond.getEnd());
    hts_itr_destroy(iter);
    return;
}

Evidence RefineDeletion::calculateFinalBreakpoint(std::map<std::pair<int32_t, int32_t>, RefineSV::MatchRead> *listPosition)
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
    std::vector<AlternativeSA> BAltSA;

    int lastscore = 0;

    for (auto const &x : *listPosition)
    {
        int maxMatchSize = getMaxIntFromVector(x.second.MatchLists);

        uint8_t maxQuality = getMaxUInt8FromVector(x.second.MapQLists);

        if (maxMatchSize < 20)
        {
            continue;
        }

        bFrequency = x.second.NumberOfMatchRead;

        if (bFrequency <= 1)
        {
            continue;
        }

        if (x.first.second - x.first.first < 20)
        {
            continue;
        }

        if (evidence.getMark() == "SR")
        {
            if (evidence.getMaxMapQ() < 60)
            {
                continue;
            }
        }

        if (getMaxUInt8FromVector(x.second.MapQLists) < 30)
        {
            continue;
        }

        if (evidence.getMaxMapQ() < 30)
        {
            continue;
        }

        if (maxMatchSize < getDivider(samplestat->getReadLength(), 1, 5, 1))
        {
            continue;
        }

        // std::cout << x.first.first << " = " << x.first.second << std::endl;

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
            BAltSA = x.second.AltSA;
        }
    }

    Evidence result;

    if (evidence.getMark() == "SR")
    {
        result.setMark("SR");
        result.setMapQList(bMapQList);

        if (evidence.getFrequency() >= 2 && (bPos == 0 || bEnd == 0))
        {
            if (bPos == 0 || bEnd == 0)
            {
                bPos = evidence.getPos();
                bEnd = evidence.getEnd();
                bHit = evidence.getFrequency();
                result.setMapQList(*evidence.getMapQVector());
            }
        }
    }
    else
    {
        result.setMapQList(bMapQList);
    }

    result.setPos(bPos);
    result.setEnd(bEnd);
    result.setFrequency(bHit);
    result.setRPMapQ(*evidence.getMapQVector());
    result.setChr(evidence.getChr());
    result.setEndChr(evidence.getEndChr());
    result.LNGMATCH = bMaxMatchSize;
    result.setVariantType("DEL");
    // result.AltSA = BAltSA;
    for (auto n : BAltSA)
    {
        result.addAlterSA(n.chr, n.pos);
    }

    if (bHit <= 1)
    {
        return result;
    }

    if (bPos == 0)
    {
        return result;
    }
    if (bEnd == 0)
    {
        return result;
    }

    if (bEnd - bPos > 1000000)
    {
        return result;
    }

    result.setQuailtyPass(true);

    return result;
}

int RefineDeletion::getNumberMapQ(std::vector<uint8_t> mapqlist, uint8_t start, uint8_t end)
{
    int count = 0;
    for (auto n : mapqlist)
    {
        if (n >= start && n <= end)
        {
            count++;
        }
    }
    return count;
}