#include "specifyevidenceinversion.h"

SpecifyEvidenceInversion::SpecifyEvidenceInversion()
{
    svtype = "INV";
}

void SpecifyEvidenceInversion::updateRead()
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

    int32_t diff = (readparser.getMatePos() + readparser.getLengthSequence()) - readparser.getPos();

    if (diff < 0)
    {
        return;
    }

    //Limit SVLEN
    if (diff > 1000000)
    {
        return;
    }

    if (readparser.isReverse() && readparser.isMateReverse())
    {
        checkRange();
        return;
    }

    if (!readparser.isReverse() && !readparser.isMateReverse())
    {
        checkRange();
        return;
    }

    return;
}

bool SpecifyEvidenceInversion::incrementSVFreq(int32_t overlappedpos, int32_t overlappedsvlength, int32_t pos, int32_t mpos)
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

void SpecifyEvidenceInversion::checkRange()
{
    bool added = false;
    int32_t merge = int32_t(samplestat->getAverageSampleStat()) + int32_t(samplestat->getSDSampleStat()) + (samplestat->getReadLength());
    added = incrementSVFreq(merge, merge, currentPos, currentMPos);

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

void SpecifyEvidenceInversion::proveEvidence(int index)
{
    int32_t plus = samplestat->getAverageSampleStat() + (samplestat->getSDSampleStat()) + samplestat->getReadLength();
    if (currentPos - plus > preCollectSV.at(index).getPosDiscordantRead())
    {
        if (filterEvidence(&preCollectSV.at(index)))
        {
            calculateVCF(&preCollectSV.at(index));

            // if (preCollectSV.at(index).getSvLength() > 1 && preCollectSV.at(index).getSvLength() < 1000000)
            // {
            finalEvidence.push_back(preCollectSV.at(index));
            // }
            preCollectSV.erase(preCollectSV.begin() + index);
            writeBufferEvidenceFile();
        }
        else
        {
            preCollectSV.erase(preCollectSV.begin() + index);
        }
    }
}

void SpecifyEvidenceInversion::calculateVCF(Evidence *evidence)
{
    int32_t firstPos = 0;
    int32_t lastPos = 0;
    int32_t avgPos = 0;
    int32_t firstEndDis = 0;
    int32_t firstEnd = 0;
    int32_t lastEnd = 0;
    int32_t avgEnd = 0;

    int32_t difflengthPos = (samplestat->getAverageSampleStat()) + (samplestat->getSDSampleStat() * 3) + (samplestat->getReadLength());
    int32_t difflengthEnd = (samplestat->getAverageSampleStat()) + (samplestat->getSDSampleStat() * 3) + (samplestat->getReadLength());

    firstPos = evidence->getLastPosDiscordantRead();
    lastPos = evidence->getLastPosDiscordantRead();
    firstEnd = evidence->getEndDiscordantRead();
    lastEnd = evidence->getEndDiscordantRead();

    if (lastPos<firstEnd) {
        evidence->setPos(lastPos);
        evidence->setEnd(firstEnd);
    } else {
        evidence->setPos(firstEnd);
        evidence->setEnd(lastPos);
    }

    evidence->setCiPosLeft(-difflengthPos);
    evidence->setCiPosRight(difflengthPos);
   
    evidence->setCiEndLeft(-difflengthEnd);
    evidence->setCiEndRight(difflengthEnd);
}

bool SpecifyEvidenceInversion::filterEvidence(Evidence *evidence)
{
    int32_t svLength = evidence->getEndDiscordantRead() - evidence->getPosDiscordantRead() - samplestat->getAverageSampleStat();

    if (svLength > 1000000)
    {
        return false;
    }

    if (evidence->getFrequency() <= 1)
    {
        return false;
    }
 
    return true;
}

void SpecifyEvidenceInversion::checkProveEvidence()
{
    for (int i = 0; i < preCollectSV.size(); i++)
    {
        proveEvidence(i);
    }
}

void SpecifyEvidenceInversion::done()
{
    checkProveEvidence();
    writeFinalEvidenceAndClear();
}