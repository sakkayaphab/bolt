#include "refiningtranslocation.h"
#include "smithwaterman.h"

RefiningTranslocation::RefiningTranslocation()
{
    variantresult.setVariantType("BND");
}

void RefiningTranslocation::execute()
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

void RefiningTranslocation::first()
{
    std::string findRange = convertRangeToString(evidence.getChr(), evidence.getPos() + evidence.getCiPosLeft(), evidence.getPos() + evidence.getCiPosRight());

    const char *range = findRange.c_str();
    // const char *mChr = evidence.getChr().c_str();
    //         std::cout << range  << std::endl;
    refineStartToEnd(range);
}

void RefiningTranslocation::refineStartToEnd(const char *range)
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
    std::cout << "ref pos : " << positionStartReference
              << " ref end : " << evidence.getEnd() + evidence.getCiEndRight() << std::endl;
    std::cout << seqrefString << std::endl;

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

        auto cigar = readparser.getCigar();
        // if (!cigar.size() == 2)
        // {
        //     // continue;
        // }

        if ((cigar.at(0).getOperatorName() == 'S' && cigar.at(1).getOperatorName() == 'M'))
        // if (true)
        {
            // continue;
            std::string fullRead = readparser.getSequence();

            // std::vector<SmithWaterman::ScoreAlignment> result = alignment.alignDuplicationTargetAtStart(&fullRead);
            std::vector<StringSearch::Score> result = ssa.alignDuplicationTargetAtStart(&fullRead);
            for (auto n : result)
            {
               

                int32_t mPos = n.endseq + readparser.getPosOfSeq();
                int32_t mEnd = n.end;

                // std::cout << "mPos : " << mPos << std::endl;
                // std::cout << "mtEnd : " << mEnd << std::endl;
                // std::cout << "pattern : " << n.pattern << std::endl;
                auto rangeMapping = n.endseq;

                for (int i = 0; i < rangeMapping; i++)
                {
                    int32_t tempPos = mPos - i;
                    int32_t tempEnd = mEnd - i;

                    if (i > rangeMapping - 6)
                    {
                        break;
                    }

                    listPosition[std::make_pair(tempPos, tempEnd)].NumberOfMatchRead++;
                    listPosition[std::make_pair(tempPos, tempEnd)].MatchLists.push_back(rangeMapping - i);
                    listPosition[std::make_pair(tempPos, tempEnd)].MapQLists.push_back(readparser.getMapQuality());
                }
            }
        }
        else if ((cigar.at(cigar.size() - 1).getOperatorName() == 'S'))
        //  if (true)
        {
            continue;
            std::string fullRead = readparser.getSequence();
            // std::cout << "-------------" << std::endl;
            // std::cout << " sc : " << cigar.at(cigar.size() - 1).getLength() << std::endl;
            // std::vector<SmithWaterman::ScoreAlignment> result = alignment.alignTranslocationTargetAtStartSCE(&fullRead);
             std::vector<StringSearch::Score> result = ssa.alignTranslocationTargetAtStartSCE(&fullRead);
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

                for (int i = 0; i < rangeMapping; i++)
                {
                    int32_t tempPos = mPos + i;
                    int32_t tempEnd = mEnd + i;

                    if (i > rangeMapping - 6)
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