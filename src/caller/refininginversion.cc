#include "refininginversion.h"
#include "smithwaterman.h"

RefiningInversion::RefiningInversion()
{
    variantresult.setVariantType("INV");
}

void RefiningInversion::execute()
{
    variantresult.setChr(evidence.getChr());
    variantresult.setEndChr(evidence.getEndChr());
    prepareBamReader();

    first();
    if (variantresult.isQuailtyPass())
    {
        return;
    }
    second();
    if (variantresult.isQuailtyPass())
    {
        return;
    }
}

void RefiningInversion::first()
{
    std::string findRange = convertRangeToString(evidence.getChr(), evidence.getPos() + evidence.getCiPosLeft(), evidence.getPos() + evidence.getCiPosRight());

    const char *range = findRange.c_str();
    // const char *mChr = evidence.getChr().c_str();
    std::cout << range << "/" << samplestat->getReadLength() << std::endl;
    refineStartToEnd(range);
}

void RefiningInversion::second()
{
    std::string findRange = convertRangeToString(evidence.getChr(), evidence.getEnd() + evidence.getCiEndLeft(), evidence.getEnd() + evidence.getCiEndRight());

    const char *range = findRange.c_str();
    // const char *mChr = evidence.getChr().c_str();
    std::cout << range << "/" << samplestat->getReadLength() << std::endl;
    refineEndToStart(range);
}

void RefiningInversion::refineStartToEnd(const char *range)
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
    // std::cout << "chr : " << evidence.getEndChr() << ", posDiscordantRead : " << evidence.getEndDiscordantRead() - 500 << ", end : " << evidence.getLastEndDiscordantRead() + 500 << std::endl;

    replaceSeqToUppercase(&seqrefString);
    std::map<std::pair<int32_t, int32_t>, MatchRead> listPosition;
    std::map<int32_t, int> SCReadLists;

    // Alignment alignment(&seqrefString);
    // alignment.genarateMatrix("INVSTART");
    // alignment.setPosReference(positionStartReference);

    if (seqrefString.size() < 50)
    {
        return;
    }

    StringSearchAlignment ssa;
    ssa.setReference(seqrefString);
    ssa.setSVType("INVSTART");
    ssa.setPosReference(positionStartReference);
    ssa.buildReference();

    while (sam_itr_next(inFile, iter, read) >= 0)
    {
        if (readparser.isUnmapped())
        {
            continue;
        }

        auto cigar = readparser.getCigar();
        // if (!cigar.size() == 2)
        // {
        //     continue;
        // }

        // if ((cigar.at(0).getOperatorName() == 'S' && cigar.at(1).getOperatorName() == 'M'))
        if (true)
        {
            std::string fullRead = readparser.getSequence();
            // std::vector<SmithWaterman::ScoreAlignment> result = alignment.alignInversionTargetAtStartSCS(&fullRead);
            StringSearchConfig ssc;
            std::vector<StringSearch::Score> result = ssa.alignInversionTargetAtStartSCS(&fullRead,&ssc);
            for (auto n : result)
            {
                if (n.matchCount <= 8)
                {
                    continue;
                }

                int32_t mPos = readparser.getPosOfSeq() + n.endseq;
                int32_t mEnd = n.pos;

                if (mEnd <= mPos + 2)
                {
                    continue;
                }
                // std::cout << "mPos : " << mPos << std::endl;
                // std::cout << "mtEnd : " << mEnd << std::endl;
                // std::cout << "pattern : " << n.pattern << std::endl;
                auto rangeMapping = n.endseq;

                for (int i = 0; i < rangeMapping; i++)
                {
                    int32_t tempPos = mPos - i;
                    int32_t tempEnd = mEnd + i;

                    if (i > 6)
                    {
                        break;
                    }

                    listPosition[std::make_pair(tempPos, tempEnd)].NumberOfMatchRead++;
                    listPosition[std::make_pair(tempPos, tempEnd)].MatchLists.push_back(rangeMapping - i);
                    listPosition[std::make_pair(tempPos, tempEnd)].MapQLists.push_back(readparser.getMapQuality());
                }
            }
        }
        // else if ((cigar.at(cigar.size() - 1).getOperatorName() == 'S'))
        if (true)
        {
            if (cigar.at(cigar.size() - 1).getLength() < 10)
            {
                continue;
            }

            std::string fullRead = readparser.getSequence();
            // std::cout << "-------------" << std::endl;
            // std::cout << fullRead << std::endl;
            // std::cout << "pos end : " << readparser.getEnd() << std::endl;
            // std::cout << " sc : " << cigar.at(cigar.size() - 1).getLength() << std::endl;
            // std::vector<SmithWaterman::ScoreAlignment> result = alignment.alignInversionTargetAtStartSCE(&fullRead);
            StringSearchConfig ssc;
            std::vector<StringSearch::Score> result = ssa.alignInversionTargetAtStartSCE(&fullRead,&ssc);
            for (auto n : result)
            {
                if (n.matchCount <= 8)
                {
                    continue;
                }

                int32_t mPos = readparser.getPosOfSeq() + n.posseq;
                int32_t mEnd = n.end;

                if (mEnd <= mPos + 2)
                {
                    continue;
                }

                // std::cout << "mPos : " << mPos << std::endl;
                // std::cout << "mEnd : " << mEnd << std::endl;
                // std::cout << "pattern : " << n.pattern << std::endl;
                // std::cout << "---------- MD TAG ----------" << std::endl;

                auto rangeMapping = n.endseq - n.posseq;

                for (int i = 0; i < rangeMapping; i++)
                {
                    int32_t tempPos = mPos + i;
                    int32_t tempEnd = mEnd - i;

                    if (i > 6)
                    {
                        break;
                    }

                    listPosition[std::make_pair(tempPos, tempEnd)].NumberOfMatchRead++;
                    listPosition[std::make_pair(tempPos, tempEnd)].MatchLists.push_back(rangeMapping - i);
                    listPosition[std::make_pair(tempPos, tempEnd)].MapQLists.push_back(readparser.getMapQuality());
                    // if (listPosition[std::make_pair(tempPos, tempEnd)].Sequence.size() < n.scorepattern - 1)
                    // {
                    //     listPosition[std::make_pair(tempPos, tempEnd)].maxMatchSequence = n.pattern.size() - 1;
                    // }
                }
            }
        }
    }

    calculateFinalBreakpoint(&listPosition);

    hts_itr_destroy(iter);
    return;
}

void RefiningInversion::refineEndToStart(const char *range)
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
    // alignment.genarateMatrix("INVEND");

    // alignment.setPosReference(positionStartReference);

    if (seqrefString.size() < 50)
    {
        return;
    }

    StringSearchAlignment ssa;
    ssa.setReference(seqrefString);
    ssa.setSVType("INVEND");
    ssa.setPosReference(positionStartReference);
    ssa.buildReference();

    while (sam_itr_next(inFile, iter, read) >= 0)
    {
        if (readparser.isUnmapped())
        {
            continue;
        }

        auto cigar = readparser.getCigar();
        // if (!cigar.size() == 2)
        // {
        //     continue;
        // }

        if ((cigar.at(0).getOperatorName() == 'S' && cigar.at(1).getOperatorName() == 'M'))
        // if (true)
        {
            // continue;
            std::string fullRead = readparser.getSequence();
            // std::vector<SmithWaterman::ScoreAlignment> result = alignment.alignInversionTargetAtStartSCS(&fullRead);
            StringSearchConfig ssc;
            std::vector<StringSearch::Score> result = ssa.alignInversionTargetAtStartSCS(&fullRead,&ssc);
            for (auto n : result)
            {
                if (n.matchCount <= 8)
                {
                    continue;
                }

                int32_t mPos = n.pos;
                int32_t mEnd = readparser.getPosOfSeq() + n.endseq;

                if (mPos >= mEnd + 2)
                {
                    continue;
                }

                // std::cout << "mPos : " << mPos << std::endl;
                // std::cout << "mtEnd : " << mEnd << std::endl;
                // std::cout << "pattern : " << n.pattern << std::endl;
                auto rangeMapping = n.endseq;

                for (int i = 0; i < rangeMapping; i++)
                {
                    int32_t tempPos = mPos - i;
                    int32_t tempEnd = mEnd + i;

                    if (i > 6)
                    {
                        break;
                    }

                    listPosition[std::make_pair(tempPos, tempEnd)].NumberOfMatchRead++;
                    listPosition[std::make_pair(tempPos, tempEnd)].MatchLists.push_back(rangeMapping - i);
                    listPosition[std::make_pair(tempPos, tempEnd)].MapQLists.push_back(readparser.getMapQuality());
                }
            }
        }

        // if ((cigar.at(cigar.size() - 1).getOperatorName() == 'S'))
        if (true)
        {
            // continue;
            // if (cigar.at(cigar.size() - 1).getLength() < 10)
            // {
            //     continue;
            // }

            std::string fullRead = readparser.getSequence();
            // std::cout << "-------------" << std::endl;
            // std::cout << fullRead << std::endl;
            // std::cout << "pos end : " << readparser.getEnd() << std::endl;
            // std::cout << " sc : " << cigar.at(cigar.size() - 1).getLength() << std::endl;
            StringSearchConfig ssc;
            std::vector<StringSearch::Score> result = ssa.alignInversionTargetAtStartSCE(&fullRead,&ssc);
            // std::vector<SmithWaterman::ScoreAlignment> result = alignment.alignInversionTargetAtStartSCE(&fullRead);
            for (auto n : result)
            {
                if (n.matchCount <= 8)
                {
                    continue;
                }

                int32_t mPos = n.end;
                int32_t mEnd = readparser.getPosOfSeq() + n.posseq;

                if (mPos >= mEnd + 2)
                {
                    continue;
                }

                // std::cout << "mPos : " << mPos << std::endl;
                // std::cout << "mEnd : " << mEnd << std::endl;
                // std::cout << "pattern : " << n.pattern << std::endl;
                // std::cout << "---------- MD TAG ----------" << std::endl;

                auto rangeMapping = n.endseq - n.posseq;

                for (int i = 0; i < rangeMapping; i++)
                {
                    int32_t tempPos = mPos - i;
                    int32_t tempEnd = mEnd + i;

                    if (i > 6)
                    {
                        break;
                    }

                    listPosition[std::make_pair(tempPos, tempEnd)].NumberOfMatchRead++;
                    listPosition[std::make_pair(tempPos, tempEnd)].MatchLists.push_back(rangeMapping - i);
                    listPosition[std::make_pair(tempPos, tempEnd)].MapQLists.push_back(readparser.getMapQuality());
                    // if (listPosition[std::make_pair(tempPos, tempEnd)].Sequence.size() < n.scorepattern - 1)
                    // {
                    //     listPosition[std::make_pair(tempPos, tempEnd)].maxMatchSequence = n.pattern.size() - 1;
                    // }
                }
            }
        }
    }

    calculateFinalBreakpoint(&listPosition);
    std::cout << variantresult.getPos() << std::endl;
    hts_itr_destroy(iter);
    return;
}