#include "specifyingevidencetranslocation.h"

SpecifyingEvidenceTranslocation::SpecifyingEvidenceTranslocation()
{
    svtype = "BND";
}

void SpecifyingEvidenceTranslocation::updateRead()
{
    currentPos = read->core.pos + 1;
    currentMPos = read->core.mpos + 1;

    if (readparser.isPairOnSameChromosome())
    {
        return;
    }

    if (read->core.flag & BAM_FREAD2)
    {
        return;
    }

    if (!(!(read->core.flag & BAM_FREVERSE) && (read->core.flag & BAM_FMREVERSE)))
    {
        return;
    }

    //Pass filter
    checkRange();
}

void SpecifyingEvidenceTranslocation::checkRange()
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
        else
        {
            return;
        }
    }
}

void SpecifyingEvidenceTranslocation::proveEvidence(int index)
{
    if (currentPos - 1000 > preCollectSV.at(index).getPosDiscordantRead())
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
    int32_t lastPosDis = evidence->getLastPosDiscordantRead();
    int32_t posDis = evidence->getPosDiscordantRead();
    int32_t firstPos = lastPosDis - (samplestat->getReadLength() * 4) - ((lastPosDis - posDis) / 2);
    int32_t lastPos = lastPosDis + (samplestat->getReadLength() * 4) + samplestat->getSDSampleStat();
    int32_t avgPos = (firstPos + lastPos) / 2;
    evidence->setPos(avgPos);
    evidence->setCiPosLeft(firstPos - avgPos);
    evidence->setCiPosRight(lastPos - avgPos);

    int32_t firstEndDis = evidence->getEndDiscordantRead();
    int32_t firstEnd = firstEndDis - (samplestat->getReadLength() * 4) - 300;
    int32_t lastEnd = firstEndDis + (samplestat->getReadLength() * 4);
    int32_t avgEnd = (firstEnd + lastEnd) / 2;
    evidence->setEnd(avgEnd);
    evidence->setCiEndLeft(firstEnd - avgEnd);
    evidence->setCiEndRight(lastEnd - avgEnd);
}

bool SpecifyingEvidenceTranslocation::filterEvidence(Evidence *evidence)
{

    if (evidence->getMaxMapQ() == 0)
    {

        return false;
    }

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