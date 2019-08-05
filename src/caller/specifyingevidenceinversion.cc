#include "specifyingevidenceinversion.h"

SpecifyingEvidenceInversion::SpecifyingEvidenceInversion()
{
    svtype = "INV";
}

void SpecifyingEvidenceInversion::updateRead()
{
    currentPos = read->core.pos + 1;
    currentMPos = read->core.mpos + 1;

    //    checkErrorEvidence();

    if (!readparser.isPairOnSameChromosome())
    {
        return;
    }

    if (read->core.flag & BAM_FREAD1)
    {
        return;
    }

    if (readparser.getPos() > readparser.getMatePos())
    {
        return;
    }

    int32_t diff = (readparser.getMatePos() + readparser.getLengthSequence()) - readparser.getPos();

    if (diff < 0)
    {
        return;
    }

    //Limit SVLEN
    if (diff > 100000)
    {
        return;
    }

    if (!(readparser.isReverse() && readparser.isMateReverse()))
    {
        return;
    }

    checkRange();
}

void SpecifyingEvidenceInversion::checkRange()
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
            // std::cout << "add new evidence" << std::endl;
        }
        else
        {
            return;
        }
    }
}

void SpecifyingEvidenceInversion::proveEvidence(int index)
{
    if (currentPos-1000 > preCollectSV.at(index).getPosDiscordantRead())
    {
        if (filterEvidence(&preCollectSV.at(index)))
        {
            calculateVCF(&preCollectSV.at(index));

            if (preCollectSV.at(index).getSvLength() > 1 && preCollectSV.at(index).getSvLength() < 50000)
            {
                finalEvidence.push_back(preCollectSV.at(index));
            }
            preCollectSV.erase(preCollectSV.begin() + index);
            writeBufferEvidenceFile();
        }
        else
        {
            preCollectSV.erase(preCollectSV.begin() + index);
        }
    }
}

void SpecifyingEvidenceInversion::calculateVCF(Evidence *evidence)
{
    int32_t lastPosDis = evidence->getLastPosDiscordantRead();
    int32_t posDis = evidence->getPosDiscordantRead();
    int32_t firstPos = lastPosDis - (samplestat->getReadLength() * 5) - samplestat->getSDSampleStat();
    int32_t lastPos = lastPosDis + (samplestat->getReadLength() * 2) + samplestat->getSDSampleStat();
    int32_t avgPos = (firstPos + lastPos) / 2;
    evidence->setPos(avgPos);
    evidence->setCiPosLeft(firstPos - avgPos);
    evidence->setCiPosRight(lastPos - avgPos);

    int32_t firstEndDis = evidence->getEndDiscordantRead();
    int32_t firstEnd = firstEndDis - (samplestat->getReadLength() * 5) - samplestat->getSDSampleStat();
    int32_t lastEnd = firstEndDis + (samplestat->getReadLength() * 2) + samplestat->getSDSampleStat();
    int32_t avgEnd = (firstEnd + lastEnd) / 2;
    evidence->setEnd(avgEnd);
    evidence->setCiEndLeft(firstEnd - avgEnd);
    evidence->setCiEndRight(lastEnd - avgEnd);
}

bool SpecifyingEvidenceInversion::filterEvidence(Evidence *evidence)
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

void SpecifyingEvidenceInversion::checkProveEvidence()
{
    for (int i = 0; i < preCollectSV.size(); i++)
    {
        proveEvidence(i);
    }
}

void SpecifyingEvidenceInversion::done()
{
    checkProveEvidence();
    writeFinalEvidenceAndClear();
}