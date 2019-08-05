#include "specifyingevidencetandemduplication.h"

SpecifyingEvidenceTandemDuplication::SpecifyingEvidenceTandemDuplication()
{
    svtype = "DUP";
}

void SpecifyingEvidenceTandemDuplication::updateRead()
{
    currentPos = read->core.pos + 1;
    currentMPos = read->core.mpos + 1;

    if (!readparser.isPairOnSameChromosome())
    {
        return;
    }

    if (!readparser.isReverse())
    {
        return;
    }

    if (readparser.isMateReverse())
    {
        return;
    }

    if (readparser.isFirstRead()) {
        return;
    }

    //limit SVLEN
    if (readparser.getMatePos() - readparser.getPos() > 50000)
    {
        return;
    }

    if (readparser.getPos() > readparser.getMatePos())
    {
        return;
    }

    checkRange();
}

void SpecifyingEvidenceTandemDuplication::checkRange()
{
    bool added;

    int positionOverlapped = findOverlapped(2000, currentPos, currentMPos);
    if (positionOverlapped >= 0)
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
    checkProveEvidence();

    if (!added)
    {
        if (currentPos < currentMPos)
        {
            Evidence evidence;
            evidence.setVariantType(svtype);
            evidence.setChr(readparser.getChromosomeNameString());
            evidence.setEndChr(readparser.getChromosomeNameString());
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
        else
        {
            return;
        }
    }
}

void SpecifyingEvidenceTandemDuplication::proveEvidence(int index)
{
    if (currentPos - 1000 > preCollectSV.at(index).getPosDiscordantRead())
    {
        if (filterEvidence(&preCollectSV.at(index)))
        {
            calculateVCF(&preCollectSV.at(index));
            finalEvidence.push_back(preCollectSV.at(index));
            preCollectSV.erase(preCollectSV.begin() + index);
            writeBufferEvidenceFile();
        }
        else
        {
            preCollectSV.erase(preCollectSV.begin() + index);
        }
    }
}

void SpecifyingEvidenceTandemDuplication::calculateVCF(Evidence *evidence)
{
    int32_t lastPosDis = evidence->getLastPosDiscordantRead();
    int32_t posDis = evidence->getPosDiscordantRead();
    int32_t firstPos = posDis - (samplestat->getReadLength() * 6) - samplestat->getSDSampleStat();
    int32_t lastPos = posDis + (samplestat->getReadLength() * 1) + samplestat->getSDSampleStat();
    int32_t avgPos = (firstPos + lastPos) / 2;
    evidence->setPos(avgPos);
    evidence->setCiPosLeft(firstPos - avgPos);
    evidence->setCiPosRight(lastPos - avgPos);

    int32_t firstEndDis = evidence->getEndDiscordantRead();
    int32_t lastEndDis = evidence->getLastEndDiscordantRead();
    // int32_t firstEnd = firstEndDis - (samplestat->getReadLength() * 1) + samplestat->getSDSampleStat();
        int32_t firstEnd = lastEndDis - (samplestat->getReadLength() * 1) + samplestat->getSDSampleStat();

    int32_t lastEnd = lastEndDis + (samplestat->getReadLength() * 6) + samplestat->getSDSampleStat();
    int32_t avgEnd = (firstEnd + lastEnd) / 2;
    evidence->setEnd(avgEnd);
    evidence->setCiEndLeft(firstEnd - avgEnd);
    evidence->setCiEndRight(lastEnd - avgEnd);
}

bool SpecifyingEvidenceTandemDuplication::filterEvidence(Evidence *evidence)
{
    if (evidence->getEndDiscordantRead() - evidence->getLastPosDiscordantRead() < 0)
    {
        return false;
    }


    if (evidence->getFrequency() < 2)
    {
        return false;
    }

    return true;
}

void SpecifyingEvidenceTandemDuplication::checkProveEvidence()
{
    for (int i = 0; i < preCollectSV.size(); i++)
    {
        proveEvidence(i);
    }
}

void SpecifyingEvidenceTandemDuplication::done()
{
    checkProveEvidence();
    writeFinalEvidenceAndClear();
}