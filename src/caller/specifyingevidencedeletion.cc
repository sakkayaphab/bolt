#include "specifyingevidencedeletion.h"

SpecifyingEvidenceDeletion::SpecifyingEvidenceDeletion()
{
    svtype = "DEL";
}

void SpecifyingEvidenceDeletion::updateRead()
{
    currentPos = read->core.pos + 1;
    currentMPos = read->core.mpos + 1;

    if (!readparser.isPairOnSameChromosome())
    {
        return;
    }

    if (readparser.getPos() > readparser.getMatePos())
    {
        return;
    }

    if (readparser.isMateUnmapped())
    {
        return;
    }

    if (readparser.isSecondRead())
    {
        return;
    }

    if (readparser.isReverse()) {
        return;
    }

    if (!readparser.isMateReverse())
    {
        return;
    }

    int32_t insertSizeFirstRead = (readparser.getMatePos() + readparser.getLengthSequence()) - readparser.getPos();

    if (readparser.isMateUnmapped())
    {
        return;
    }

    if (!readparser.isPairOnSameChromosome())
    {
        return;
    }

    if (insertSizeFirstRead < 0)
    {
        return;
    }

    if (readparser.isReverse() && readparser.isMateReverse())
    {
        return;
    }

    if (!readparser.isReverse() && !readparser.isMateReverse())
    {
        return;
    }

    if (insertSizeFirstRead > samplestat->getMedianSampleStat() + samplestat->getSDSampleStat())
    {
        checkRange();
        return;
    }

    // if (insertSizeFirstRead > int32_t(samplestat->getMedianSampleStat()) + int32_t(samplestat->getSDSampleStat()*0.5))
    // {
    //     if (readparser.getMapQuality() >= 20)
    //     {
    //         checkRange();
    //         return;
    //     }
    // }
}

int32_t SpecifyingEvidenceDeletion::getSVLength()
{
    return currentMPos - currentPos - (samplestat->getMedianSampleStat() + samplestat->getSDSampleStat());
}

void SpecifyingEvidenceDeletion::checkRange()
{
    bool added = false;
    // int32_t diff = currentMPos-currentPos+samplestat->getReadLength()-samplestat->getMedianSampleStat()+(2*currentPos+samplestat->getReadLength())+(2* samplestat->getSDSampleStat());

    int32_t merge = int32_t(samplestat->getMedianSampleStat()) + int32_t(samplestat->getSDSampleStat()) + (samplestat->getReadLength());
    // std::cout << "merge : " << merge  << ", " << samplestat->getMedianSampleStat() << " , " << int32_t(samplestat->getSDSampleStat()) << " , "
    // << samplestat->getReadLength()
    // << std::endl;
    // if (getSVLength() < 500)
    // {
    //     added = incrementSVFreq(merge, merge, currentPos, currentMPos);
    // }
    // else if (getSVLength() >= 500 && getSVLength() < 1000)
    // {
    //     added = incrementSVFreq(merge, merge, currentPos, currentMPos);
    // }
    // else if (getSVLength() >= 1000 && getSVLength() < 2000)
    // {
    //     added = incrementSVFreq(merge, merge, currentPos, currentMPos);
    // }
    // else if (getSVLength() >= 2000)
    // {
    added = incrementSVFreq(merge, merge, currentPos, currentMPos);
    // }

    // int positionOverlapped = findOverlapped(2000, currentPos, currentMPos);
    // if (positionOverlapped >= 0)
    // {
    //     preCollectSV.at(positionOverlapped).addAssociateRead(currentPos, currentMPos);

    //     if (preCollectSV.at(positionOverlapped).getLastPosDiscordantRead() < currentPos)
    //     {
    //         preCollectSV.at(positionOverlapped).setLastPosDiscordantRead(currentPos);
    //     }

    //     if (preCollectSV.at(positionOverlapped).getEndDiscordantRead() > currentMPos)
    //     {
    //         preCollectSV.at(positionOverlapped).setEndDiscordantRead(currentMPos);
    //     }

    //     if (preCollectSV.at(positionOverlapped).getLastEndDiscordantRead() < currentMPos)
    //     {
    //         preCollectSV.at(positionOverlapped).setLastEndDiscordantRead(currentMPos);
    //     }

    //     preCollectSV.at(positionOverlapped).addMapQ(readparser.getMapQuality());
    //     added = true;
    //     preCollectSV.at(positionOverlapped).incrementFrequency();
    // }
    // else
    // {
    //     int positionOverlappedOnlyPos = findOverlappedOnlyPos(1000, currentPos);
    //     if (positionOverlappedOnlyPos >= 0)
    //     {
    //         preCollectSV.at(positionOverlappedOnlyPos).addAssociateRead(currentPos, currentMPos);
    //     }
    // }

    checkProveEvidence();

    if (!added)
    {
        Evidence evidence;
        evidence.setVariantType(svtype);
        evidence.setChr(readparser.getChromosomeNameString());
        evidence.setEndChr(readparser.getChromosomeNameString());
        evidence.setPosDiscordantRead(currentPos);
        evidence.setEndDiscordantRead(currentMPos);
        evidence.setLastPosDiscordantRead(currentPos);
        evidence.setLastEndDiscordantRead(currentMPos);

        evidence.incrementFrequency();
        evidence.setForwardDirection(true);
        evidence.addAssociateRead(currentPos, currentMPos);

        evidence.addMapQ(readparser.getMapQuality());
        preCollectSV.push_back(evidence);
    }
}

bool SpecifyingEvidenceDeletion::incrementSVFreq(int32_t overlappedpos, int32_t overlappedsvlength, int32_t pos, int32_t mpos)
{
    bool added;
    for (int i = 0; i < preCollectSV.size(); i++)
    {
        if (checkBetween(pos, preCollectSV.at(i).getPosDiscordantRead(), overlappedpos) && checkBetween(mpos,
                                                                                                        preCollectSV.at(i).getEndDiscordantRead(),
                                                                                                        overlappedsvlength))
        {
            preCollectSV.at(i).addAssociateRead(currentPos, currentMPos);

            if (preCollectSV.at(i).getLastPosDiscordantRead() < currentPos)
            {
                preCollectSV.at(i).setLastPosDiscordantRead(currentPos);
            }

            if (preCollectSV.at(i).getEndDiscordantRead() > currentMPos)
            {
                preCollectSV.at(i).setEndDiscordantRead(currentMPos);
            }

            if (preCollectSV.at(i).getLastEndDiscordantRead() < currentMPos)
            {
                preCollectSV.at(i).setLastEndDiscordantRead(currentMPos);
            }

            preCollectSV.at(i).addMapQ(readparser.getMapQuality());
            added = true;
            preCollectSV.at(i).incrementFrequency();
        }
    }

    return added;
}

void SpecifyingEvidenceDeletion::proveEvidence(int index)
{
    int32_t plus = samplestat->getMedianSampleStat() + (samplestat->getSDSampleStat() * 2) + samplestat->getReadLength();
    if (currentPos - plus > preCollectSV.at(index).getPosDiscordantRead())
    {
        if (filterEvidence(&preCollectSV.at(index)))
        {
            calculateVCF(&preCollectSV.at(index));

            if (preCollectSV.at(index).getSvLength() > 10 && preCollectSV.at(index).getSvLength() < 1000000)
            {
                finalEvidence.push_back(preCollectSV.at(index));
            }

            preCollectSV.erase(preCollectSV.begin() + index);
            // removeDuplicateFinalEvidence();
            writeBufferEvidenceFile();
        }
        else
        {
            preCollectSV.erase(preCollectSV.begin() + index);
        }
    }
}

void SpecifyingEvidenceDeletion::removeDuplicateFinalEvidence()
{
    int number = 0;
    std::vector<Evidence> tempEvidence;
    for (auto n : finalEvidence)
    {
        bool found;
        for (auto m : finalEvidence)
        {
            if (n.getPosDiscordantRead() == m.getPosDiscordantRead())
            {
                continue;
            }

            if (n.getPosDiscordantRead() <= m.getPosDiscordantRead() && n.getLastPosDiscordantRead() >= m.getPosDiscordantRead())
            {
                found = true;
                break;
            }

            if (n.getPosDiscordantRead() <= m.getLastPosDiscordantRead() && n.getLastEndDiscordantRead() >= m.getLastPosDiscordantRead())
            {
                found = true;
                break;
            }
        }

        if (!found)
        {
            tempEvidence.push_back(n);
        }
    }

    finalEvidence = tempEvidence;
}

void SpecifyingEvidenceDeletion::calculateVCF(Evidence *evidence)
{

    int32_t firstPos = 0;
    int32_t lastPos = 0;
    int32_t avgPos = 0;
    int32_t firstEndDis = 0;
    int32_t firstEnd = 0;
    int32_t lastEnd = 0;
    int32_t avgEnd = 0;
    // int32_t svlength = evidence->getEndDiscordantRead() - evidence->getPosDiscordantRead() - samplestat->getMedianSampleStat();
    // if (evidence->getLastPosDiscordantRead() - evidence->getPosDiscordantRead() < 0)
    // {
    //     std::cout << "getLastPosDiscordantRead" << std::endl;
    //     std::cout << evidence->getLastPosDiscordantRead() << " == " << evidence->getPosDiscordantRead() << std::endl;
    // }
    // if (evidence->getLastEndDiscordantRead() - evidence->getEndDiscordantRead() < 0)
    // {
    //     std::cout << "getLastEndDiscordantRead" << std::endl;
    //     std::cout << evidence->getLastEndDiscordantRead() << " == " << evidence->getEndDiscordantRead() << std::endl;
    // }

    int32_t difflengthPos = (samplestat->getMedianSampleStat()) + (samplestat->getSDSampleStat() * 2) + (samplestat->getReadLength());
    int32_t difflengthEnd = (samplestat->getMedianSampleStat()) + (samplestat->getSDSampleStat() * 2) + (samplestat->getReadLength());
    int32_t notUsed = (samplestat->getSDSampleStat() * 2) + (samplestat->getReadLength());

    // if (difflengthEnd > 100000)
    // {
    //     // std::cout << evidence->getLastEndDiscordantRead() << " = " << evidence->getEndDiscordantRead() << std::endl;
    //     return;
    // }
    // int32_t difflengthPos =
    // int32_t difflengthEnd =

    int32_t merge = 0;
    // merge = svlength + (samplestat->getSDSampleStat()*2) + (samplestat->getReadLength()*2);
    firstPos = evidence->getLastPosDiscordantRead();
    lastPos = evidence->getLastPosDiscordantRead();
    firstEnd = evidence->getEndDiscordantRead();
    lastEnd = evidence->getEndDiscordantRead();

    // if (svlength < 500)
    // {
    //     firstPos = evidence->getPosDiscordantRead() - (samplestat->getReadLength() * 1) - samplestat->getSDSampleStat();
    //     lastPos = evidence->getLastPosDiscordantRead() + (samplestat->getReadLength() * 2) + samplestat->getSDSampleStat();
    //     firstEndDis = evidence->getEndDiscordantRead();
    //     firstEnd = firstEndDis - (samplestat->getReadLength() * 2) - samplestat->getSDSampleStat();
    //     lastEnd = evidence->getLastEndDiscordantRead() + (samplestat->getReadLength() * 1) + samplestat->getSDSampleStat();

    // // return
    // }
    // else if (svlength < 1000)
    // {
    //     firstPos = evidence->getPosDiscordantRead() - (samplestat->getReadLength() * 2) - samplestat->getSDSampleStat();
    //     lastPos = evidence->getLastPosDiscordantRead() + (samplestat->getReadLength() * 3) + samplestat->getSDSampleStat();
    //     firstEndDis = evidence->getEndDiscordantRead();
    //     firstEnd = firstEndDis - (samplestat->getReadLength() * 3) - samplestat->getSDSampleStat();
    //     lastEnd = evidence->getLastEndDiscordantRead() + (samplestat->getReadLength() * 2) + samplestat->getSDSampleStat();
    //     return;
    // }
    // else
    // {
    //     firstPos = evidence->getPosDiscordantRead() - (samplestat->getReadLength() * 4) - samplestat->getSDSampleStat();
    //     lastPos = evidence->getLastPosDiscordantRead() + (samplestat->getReadLength() * 8) + samplestat->getSDSampleStat();
    //     firstEndDis = evidence->getEndDiscordantRead();
    //     firstEnd = firstEndDis - (samplestat->getReadLength() * 8);
    //     lastEnd = evidence->getLastEndDiscordantRead() + (samplestat->getReadLength() * 4);
    //      return;
    // }

    // if (lastPos > firstEnd)
    // {
    //     lastPos = firstEndDis;
    // }

    // if (firstEnd < lastPos)
    // {
    //     firstEnd = evidence->getLastPosDiscordantRead();
    // }

    // avgPos = (firstPos + lastPos) / 2;
    // avgEnd = (firstEnd + lastEnd) / 2;
    // evidence->setPos(avgPos);
    // evidence->setCiPosLeft(firstPos - avgPos);
    // evidence->setCiPosRight(lastPos - avgPos);
    // evidence->setEnd(avgEnd);
    // evidence->setCiEndLeft(firstEnd - avgEnd);
    // evidence->setCiEndRight(lastEnd - avgEnd);

    // avgPos = (firstPos + lastPos) / 2;
    // avgEnd = (firstEnd + lastEnd) / 2;
    evidence->setPos(lastPos);
    evidence->setCiPosLeft(-notUsed);
    evidence->setCiPosRight(difflengthPos);
    evidence->setEnd(firstEnd);
    evidence->setCiEndLeft(-difflengthEnd);
    evidence->setCiEndRight(notUsed);
}

bool SpecifyingEvidenceDeletion::filterEvidence(Evidence *evidence)
{
    int32_t svLength = evidence->getEndDiscordantRead() - evidence->getPosDiscordantRead() - samplestat->getMedianSampleStat();

    if (svLength > 1000000)
    {
        return false;
    }

    if (svLength < 500)
    {
        if (evidence->getFrequency() <= 1)
        {
            return false;
        }
    }
    else if (svLength < 2000)
    {
        // if (evidence->getFrequency() <= 1)
        // {
        //     return false;
        // }
    }

    return true;
}

void SpecifyingEvidenceDeletion::checkProveEvidence()
{
    for (int i = 0; i < preCollectSV.size(); i++)
    {
        proveEvidence(i);
    }
}

void SpecifyingEvidenceDeletion::done()
{
    checkProveEvidence();
    writeFinalEvidenceAndClear();
}