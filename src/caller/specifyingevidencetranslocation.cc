#include "specifyingevidencetranslocation.h"

SpecifyingEvidenceTranslocation::SpecifyingEvidenceTranslocation()
{
    svtype = "BND";
}

void SpecifyingEvidenceTranslocation::updateRead()
{
    currentPos = read->core.pos + 1;
    currentMPos = read->core.mpos + 1;

    if (readparser.isUnmapped())
    {
        return;
    }

    if (readparser.isMateUnmapped())
    {
        return;
    }

    if (readparser.isFirstRead())
    {
        if (readparser.isReverse())
        {
            return;
        }

        if (!readparser.isMateReverse())
        {
            return;
        }
    }

    if (readparser.isSecondRead())
    {
        if (!readparser.isReverse())
        {
            return;
        }

        if (readparser.isMateReverse())
        {
            return;
        }

        return ;
    }

    if (readparser.isPairOnSameChromosome())
    {
        return;
    }

    //Pass filter
    checkRange();
}

bool SpecifyingEvidenceTranslocation::incrementSVFreq(int32_t overlappedpos, int32_t overlappedsvlength, int32_t pos, int32_t mpos)
{
    bool added = false;
    for (int positionOverlapped = 0; positionOverlapped < preCollectSV.size(); positionOverlapped++)
    {
        if (checkBetween(pos, preCollectSV.at(positionOverlapped).getPosDiscordantRead(), overlappedpos) && checkBetween(mpos,
                                                                                                                         preCollectSV.at(positionOverlapped).getEndDiscordantRead(),
                                                                                                                         overlappedsvlength))
        {
            preCollectSV.at(positionOverlapped).addAssociateRead(currentPos, currentMPos);

            if (preCollectSV.at(positionOverlapped).getLastPosDiscordantRead() < currentPos)
            {
                preCollectSV.at(positionOverlapped).setLastPosDiscordantRead(currentPos);
            }

            if (preCollectSV.at(positionOverlapped).getEndDiscordantRead() > currentMPos)
            {
                preCollectSV.at(positionOverlapped).setEndDiscordantRead(currentMPos);
            }

            if (preCollectSV.at(positionOverlapped).getLastEndDiscordantRead() < currentMPos)
            {
                preCollectSV.at(positionOverlapped).setLastEndDiscordantRead(currentMPos);
            }

            preCollectSV.at(positionOverlapped).addMapQ(readparser.getMapQuality());
            added = true;
            preCollectSV.at(positionOverlapped).incrementFrequency();
        }
    }

    return added;
}

void SpecifyingEvidenceTranslocation::checkRange()
{
    bool added;
    int32_t merge = int32_t(samplestat->getAverageSampleStat()) + int32_t(samplestat->getSDSampleStat()) + (samplestat->getReadLength());
    added = incrementSVFreq(merge, merge, currentPos, currentMPos);

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

    checkProveEvidence();

    if (!added)
    {
        Evidence evidence;
        evidence.setVariantType(svtype);
        evidence.setChr(readparser.getChromosomeNameString());
        evidence.setEndChr(readparser.getMateChromosomeNameString());
        evidence.setPosDiscordantRead(currentPos);
        evidence.setEndDiscordantRead(currentMPos);
        evidence.incrementFrequency();
        evidence.setForwardDirection(true);
        evidence.addAssociateRead(currentPos, currentMPos);
        evidence.setLastPosDiscordantRead(currentPos);
        evidence.setLastEndDiscordantRead(currentMPos);
        evidence.addMapQ(readparser.getMapQuality());
        preCollectSV.push_back(evidence);
    }
}

void SpecifyingEvidenceTranslocation::proveEvidence(int index)
{
    int32_t plus = samplestat->getAverageSampleStat() + (samplestat->getSDSampleStat()) + samplestat->getReadLength();
    if (currentPos - plus > preCollectSV.at(index).getPosDiscordantRead())
    {
        if (filterEvidence(&preCollectSV.at(index)))
        {
            calculateVCF(&preCollectSV.at(index));
            finalEvidence.push_back(preCollectSV.at(index));
            preCollectSV.erase(preCollectSV.begin() + index);

            //            filterEvidenceFinal();
            writeBufferEvidenceFile();
        }
        else
        {
            preCollectSV.erase(preCollectSV.begin() + index);
        }
    }
}

void SpecifyingEvidenceTranslocation::calculateVCF(Evidence *evidence)
{
    int32_t firstPos = 0;
    int32_t lastPos = 0;
    int32_t avgPos = 0;
    int32_t firstEndDis = 0;
    int32_t firstEnd = 0;
    int32_t lastEnd = 0;
    int32_t avgEnd = 0;
    // int32_t svlength = evidence->getEndDiscordantRead() - evidence->getPosDiscordantRead() - samplestat->getAverageSampleStat();
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

    int32_t difflengthPos = (samplestat->getAverageSampleStat()) + (samplestat->getSDSampleStat() * 2) + (samplestat->getReadLength());
    int32_t difflengthEnd = (samplestat->getAverageSampleStat()) + (samplestat->getSDSampleStat() * 2) + (samplestat->getReadLength());

    // if (difflengthEnd > 500000)
    // {
    //     std::cout << evidence->getLastEndDiscordantRead() << " = " << evidence->getEndDiscordantRead() << std::endl;
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
    evidence->setCiPosLeft(-difflengthPos);
    evidence->setCiPosRight(difflengthPos);
    evidence->setEnd(firstEnd);

    if (firstEnd-difflengthEnd<0) {
        evidence->setCiEndLeft(1);
    }else {
        evidence->setCiEndLeft(-difflengthEnd);
    }
    
    evidence->setCiEndRight(difflengthEnd);
}

bool SpecifyingEvidenceTranslocation::filterEvidence(Evidence *evidence)
{

    // if (evidence->getMaxMapQ() == 0)
    // {

    //     return false;
    // }

    if (evidence->getFrequency() >= 3)
    {
        return true;
    }

    return false;
}

void SpecifyingEvidenceTranslocation::checkProveEvidence()
{
    for (int i = 0; i < preCollectSV.size(); i++)
    {
        proveEvidence(i);
    }
}

void SpecifyingEvidenceTranslocation::done()
{
    checkProveEvidence();
    writeFinalEvidenceAndClear();
}