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

    if (readparser.isReverse())
    {
        return;
    }

    if (!readparser.isMateReverse())
    {
        return;
    }

    int32_t insertSizeFirstRead = (readparser.getMatePos() + samplestat->getReadLength()) - readparser.getPos();

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


    if (insertSizeFirstRead > samplestat->getAverageSampleStat() + (2 * samplestat->getSDSampleStat()))
    {
        checkRange();
        return;
    }
}

int32_t SpecifyingEvidenceDeletion::getSVLength()
{
    return currentMPos - currentPos - (samplestat->getAverageSampleStat() + samplestat->getSDSampleStat());
}

void SpecifyingEvidenceDeletion::checkRange()
{
    bool added = false;

    int32_t merge = int32_t(samplestat->getAverageSampleStat()) + int32_t(2 * samplestat->getSDSampleStat()) + (samplestat->getReadLength());

    added = incrementSVFreq(merge, merge, currentPos, currentMPos);

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
    bool added = false;
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
    int32_t plus = samplestat->getAverageSampleStat() + (samplestat->getSDSampleStat() * 2) + samplestat->getReadLength();
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

    int32_t merge = 0;
    lastPos = evidence->getLastPosDiscordantRead();
    firstEnd = evidence->getEndDiscordantRead();

    int32_t insertsize_max = (samplestat->getAverageSampleStat()) + (samplestat->getSDSampleStat() * 3);

    int32_t difflengthPos = (evidence->getPosDiscordantRead() + insertsize_max) - evidence->getLastPosDiscordantRead();
    int32_t difflengthEnd = insertsize_max - (evidence->getLastEndDiscordantRead() - evidence->getEndDiscordantRead());

    evidence->setPos(lastPos);
    int32_t setCiPosLeft = (evidence->getLastPosDiscordantRead() - evidence->getPosDiscordantRead());
    if (setCiPosLeft > samplestat->getReadLength())
    {
        evidence->setCiPosLeft(-(setCiPosLeft));
    }
    else
    {
        evidence->setCiPosLeft(-(samplestat->getReadLength()));
    }

    evidence->setCiPosRight(difflengthPos);
    evidence->setEnd(firstEnd);

    int32_t setCiEndRight = (evidence->getLastEndDiscordantRead() - evidence->getEndDiscordantRead());
    if (setCiEndRight > samplestat->getReadLength())
    {
        evidence->setCiEndRight((setCiEndRight));
    }
    else
    {
        evidence->setCiEndRight((samplestat->getReadLength()));
    }

    evidence->setCiEndLeft(-(difflengthEnd));
}

bool SpecifyingEvidenceDeletion::filterEvidence(Evidence *evidence)
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

    // if (evidence->getMaxMapQ() < 30)
    // {
    //     return false;
    // }

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