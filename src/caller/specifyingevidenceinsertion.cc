#include "specifyingevidenceinsertion.h"

SpecifyingEvidenceInsertion::SpecifyingEvidenceInsertion()
{
    svtype = "INS";
}

void SpecifyingEvidenceInsertion::updateRead()
{

    if (readparser.isUnmapped())
    {
        return;
    }

    currentPos = read->core.pos + 1;
    currentMPos = read->core.mpos + 1;

    if (readparser.isFirstRead() && !readparser.isReverse())
    {

        if (readparser.isMateUnmapped())
        {
            // std::cout << currentPos << " = " << currentMPos << std::endl;

            checkRange();
            return;
        }
    }

    if (readparser.isSecondRead() && readparser.isReverse())
    {
        if (readparser.isMateUnmapped())
        {
            // std::cout << currentPos << " = " << currentMPos << std::endl;

            // checkRange();
            return;
        }
    }

    if (!readparser.isPairOnSameChromosome())
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

        int diff = (readparser.getMatePos() + readparser.getLengthSequence()) - readparser.getPos();
        if (diff < 50)
        {
            return;
        }

        if (diff < samplestat->getMedianSampleStat() - samplestat->getSDSampleStat())
        {
            checkRange();
        }
        return;
    }
}

bool SpecifyingEvidenceInsertion::incrementSVFreq(int32_t overlappedpos, int32_t overlappedsvlength, int32_t pos, int32_t mpos)
{
    bool added;
    for (int positionOverlapped = 0; positionOverlapped < preCollectSV.size(); positionOverlapped++)
    {

        if (checkBetween(pos, preCollectSV.at(positionOverlapped).getPosDiscordantRead(), overlappedpos))
        {
            if (preCollectSV.at(positionOverlapped).getLastPosDiscordantRead() < currentPos)
            {
                preCollectSV.at(positionOverlapped).setLastPosDiscordantRead(currentPos);
            }

            if (currentMPos != 0)
            {
                if (preCollectSV.at(positionOverlapped).getEndDiscordantRead() > currentMPos)
                {
                    preCollectSV.at(positionOverlapped).setEndDiscordantRead(currentMPos);
                }

                if (preCollectSV.at(positionOverlapped).getLastEndDiscordantRead() < currentMPos)
                {
                    preCollectSV.at(positionOverlapped).setLastEndDiscordantRead(currentMPos);
                }
            }

            preCollectSV.at(positionOverlapped).addMapQ(readparser.getMapQuality());

            preCollectSV.at(positionOverlapped).incrementFrequency();
            added = true;
        }
    }

    return added;
}

void SpecifyingEvidenceInsertion::checkRange()
{
    bool added = false;
    int32_t merge = int32_t(samplestat->getMedianSampleStat()) + int32_t(samplestat->getSDSampleStat()) + (samplestat->getReadLength());
    added = incrementSVFreq(merge, merge, currentPos, currentMPos);

    checkProveEvidence();

    if (!added)
    {
        Evidence evidence;
        evidence.setVariantType(svtype);
        evidence.setChr(readparser.getChromosomeNameString());
        evidence.setEndChr(readparser.getChromosomeNameString());

        if (readparser.isMateUnmapped())
        {
            evidence.setComment("MATEUNMAPPED");

            evidence.setPosDiscordantRead(readparser.getPos());
            evidence.setLastPosDiscordantRead(readparser.getPos());

            evidence.setEndDiscordantRead(readparser.getPos());
            evidence.setLastEndDiscordantRead(readparser.getPos());
        }
        else
        {
            if (readparser.isFirstRead())
            {
                evidence.setPosDiscordantRead(readparser.getPos());
                evidence.setLastPosDiscordantRead(readparser.getPos());
                evidence.setEndDiscordantRead(readparser.getMatePos());
                evidence.setLastEndDiscordantRead(readparser.getMatePos());
            }
            else
            {
                evidence.setPosDiscordantRead(readparser.getMatePos());
                evidence.setLastPosDiscordantRead(readparser.getMatePos());
                evidence.setEndDiscordantRead(readparser.getPos());
                evidence.setLastEndDiscordantRead(readparser.getPos());
            }
        }

        evidence.addMapQ(readparser.getMapQuality());
        evidence.incrementFrequency();

        preCollectSV.push_back(evidence);
    }
}

void SpecifyingEvidenceInsertion::proveEvidence(int index)
{
    int32_t plus = samplestat->getMedianSampleStat() + (samplestat->getSDSampleStat()) + samplestat->getReadLength();
    if (currentPos - plus > preCollectSV.at(index).getPosDiscordantRead())
    {
        if (filterEvidence(&preCollectSV.at(index)))
        {
            calculateVCF(&preCollectSV.at(index));
            finalEvidence.push_back(preCollectSV.at(index));
            preCollectSV.erase(preCollectSV.begin() + index);
            writeBufferEvidenceFile();
            preCollectSV.clear();
        }
        else
        {
            preCollectSV.erase(preCollectSV.begin() + index);
        }
    }
}

void SpecifyingEvidenceInsertion::calculateVCF(Evidence *evidence)
{

    int32_t pos = evidence->getPosDiscordantRead();

    int32_t lastpos = evidence->getLastPosDiscordantRead();
    int32_t end = evidence->getEndDiscordantRead();
    int32_t lastend = evidence->getLastEndDiscordantRead();

    if (evidence->getComment() == "MATEUNMAPPED")
    {

        end = pos + samplestat->getMedianSampleStat() + samplestat->getSDSampleStat();
        evidence->setPos(pos);
        evidence->setEnd(end);

        evidence->setCiPosLeft(-samplestat->getReadLength()-samplestat->getMedianSampleStat());
        evidence->setCiPosRight((end - pos)+ samplestat->getMedianSampleStat());

        evidence->setCiEndLeft(pos-end);
        evidence->setCiEndRight(samplestat->getReadLength());
    }
    else
    {
        end = pos + samplestat->getMedianSampleStat() + samplestat->getSDSampleStat();
        evidence->setPos(pos);
        evidence->setEnd(end);

        evidence->setCiPosLeft(-samplestat->getReadLength()-samplestat->getMedianSampleStat());
        evidence->setCiPosRight((end - pos)+ samplestat->getMedianSampleStat());

        evidence->setCiEndLeft(pos-end);
        evidence->setCiEndRight(samplestat->getReadLength());
    }
}

bool SpecifyingEvidenceInsertion::filterEvidence(Evidence *evidence)
{
 
    if (evidence->getFrequency() <= 3)
    {
        return false;
    }

    return true;
}

void SpecifyingEvidenceInsertion::checkProveEvidence()
{
    for (int i = 0; i < preCollectSV.size(); i++)
    {
        proveEvidence(i);
    }
}

void SpecifyingEvidenceInsertion::done()
{

    std::vector<Evidence> temp;
    for (auto n : finalEvidence)
    {

        if (readdepthHelper->isRangeDisorderByMorethanRD(n.getPosDiscordantRead(), n.getEndDiscordantRead(), readdepthHelper->getAvgReadDepth() + 30))
        {
            continue;
        }
        temp.push_back(n);
    }
    finalEvidence = temp;

    checkProveEvidence();
    writeFinalEvidenceAndClear();
}