#include "specifyingevidenceinsertion.h"

SpecifyingEvidenceInsertion::SpecifyingEvidenceInsertion()
{
    svtype = "INS";
}

void SpecifyingEvidenceInsertion::updateRead()
{
    currentPos = read->core.pos + 1;
    currentMPos = read->core.mpos + 1;

    if (read->core.flag & BAM_FMUNMAP)
    {
        // checkRange();
        return;
    }

    if (!readparser.isPairOnSameChromosome())
    {
        return;
    }

    if (read->core.flag & BAM_FREAD1)
    {
        int diff = (readparser.getMatePos() + readparser.getLengthSequence()) - readparser.getPos();
        if (diff < 0)
        {
            return;
        }

        if (diff < samplestat->getMedianSampleStat() - (samplestat->getSDSampleStat() / 2))
        {
            checkRange();
            return;
        }
    }
}

void SpecifyingEvidenceInsertion::checkRange()
{
    std::string chr = bam_header->target_name[read->core.tid];
    uint32_t pos = read->core.pos + 1;

    std::string mchr = bam_header->target_name[read->core.mtid];
    uint32_t mpos = read->core.mpos + 1;

    if (mpos == 0)
    {
        mpos = pos;
    }

    if (pos == 0)
    {
        pos = mpos;
    }
    bool added;

    int32_t merge = int32_t(samplestat->getMedianSampleStat()) + int32_t(samplestat->getSDSampleStat()) + (samplestat->getReadLength());
    int positionOverlapped = findOverlapped(merge, pos, mpos);

    if (positionOverlapped >= 0)
    {
        preCollectSV.at(positionOverlapped).addAssociateRead(currentPos, currentMPos);

        if (readparser.isMateUnmapped())
        {
            if (readparser.isFirstRead())
            {
                preCollectSV.at(positionOverlapped).setForwardDirection(true);
                preCollectSV.at(positionOverlapped).NumberforwardDirection++;
                if (preCollectSV.at(positionOverlapped).getLastPosDiscordantRead() < currentPos)
                {
                    preCollectSV.at(positionOverlapped).setLastPosDiscordantRead(currentPos);
                }
            }
            else
            {
                preCollectSV.at(positionOverlapped).setBackwardDirection(true);
                preCollectSV.at(positionOverlapped).NumberbackwardDirection++;
                if (preCollectSV.at(positionOverlapped).getEndDiscordantRead() == preCollectSV.at(positionOverlapped).getPosDiscordantRead())
                {
                    preCollectSV.at(positionOverlapped).setEndDiscordantRead(pos);
                }

                if (preCollectSV.at(positionOverlapped).getLastEndDiscordantRead() < currentMPos)
                {
                    preCollectSV.at(positionOverlapped).setLastEndDiscordantRead(currentMPos);
                }
            }
        }
        else
        {
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

            if (readparser.isFirstRead())
            {
                preCollectSV.at(positionOverlapped).setForwardDirection(true);
            }
            else
            {
                preCollectSV.at(positionOverlapped).setBackwardDirection(true);
            }
        }

        //        if (readparser.isFirstRead())
        //        {
        //            preCollectSV.at(positionOverlapped).setForwardDirection(true);
        //        }
        //        else
        //        {
        //            preCollectSV.at(positionOverlapped).setFoundEvidenceAtEnd(true);
        //        }

        preCollectSV.at(positionOverlapped).addMapQ(readparser.getMapQuality());
        added = true;
        preCollectSV.at(positionOverlapped).incrementFrequency();
    }

    checkProveEvidence();

    if (!added)
    {
        Evidence evidence;
        evidence.setVariantType(svtype);
        evidence.setChr(chr);
        evidence.setEndChr(chr);
        evidence.setPosDiscordantRead(pos);
        if (readparser.isMateUnmapped())
        {
            evidence.setEndDiscordantRead(mpos);
            evidence.addAssociateRead(pos, mpos);
            evidence.setLastEndDiscordantRead(mpos);
            evidence.setComment("MATEUNMAPPED");
            if (readparser.isFirstRead())
            {
                evidence.setForwardDirection(true);
            }
            else
            {
                evidence.setBackwardDirection(true);
            }
        }
        else
        {
            evidence.setEndDiscordantRead(mpos);
            evidence.addAssociateRead(pos, mpos);
            evidence.setLastEndDiscordantRead(mpos);
        }

        if (readparser.isFirstRead())
        {
            evidence.setForwardDirection(true);
        }
        else
        {
            evidence.setBackwardDirection(true);
        }
        evidence.addMapQ(readparser.getMapQuality());
        evidence.incrementFrequency();
        evidence.setLastPosDiscordantRead(pos);
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
    // int32_t firstPos = evidence->getLastPosDiscordantRead()-50;

    int32_t pos = evidence->getPosDiscordantRead();

    int32_t lastpos = evidence->getLastPosDiscordantRead();
    int32_t end = evidence->getEndDiscordantRead();
    int32_t lastend = evidence->getLastEndDiscordantRead();

    // int32_t lastEnd = evidence->getLastPosDiscordantRead()+200;
    int32_t avgEnd = (lastpos + end) / 2;
    int32_t diff = lastend - avgEnd;
    evidence->setEnd(lastpos);
    evidence->setPos(end);

    evidence->setCiPosLeft(-samplestat->getReadLength());
    evidence->setCiPosRight(end-lastpos+samplestat->getReadLength());

    evidence->setCiEndLeft(end-lastpos-samplestat->getReadLength());
    evidence->setCiEndRight(samplestat->getReadLength());
}

bool SpecifyingEvidenceInsertion::filterEvidence(Evidence *evidence)
{
    //    readdepthHelper->getVariantVcfFormat();
    //    if (readdepthHelper->isRangeDisorderByMorethanRD(evidence->getPosDiscordantRead(),evidence->getEndDiscordantRead(),samplestat->getMedianSampleStat()+10))
    //    {
    //        return false;
    //    }

    if (evidence->getComment() == "MATEUNMAPPED")
    {

        if (evidence->getBackwardDirection() == false)
        {
            return false;
        }

        if (evidence->getForwardDirection() == false)
        {
            return false;
        }

        if (evidence->getFrequency() <= 2)
        {
            // if (evidence->getMaxMapQ() >= 30)
            // {
            return false;
            // }
        }

        return true;
    }

    // if ((int) (evidence->getFrequency())<evidence->getNumberOfZeroMapQ()-4)
    // {
    //     return false;
    // }

    //  if (evidence->getAvgMapQ() < 20)
    //     {
    //         return false;
    //     }

    //     if (evidence->getMaxMapQ() < 40)
    //     {
    //         return false;
    //     }

    if (evidence->getFrequency() < 6)
    {
        return false;
    }

    // if (evidence->getFrequency() >= samplestat->getMedianSampleStat())
    // {
    //     if (evidence->getMaxMapQ() >= 60)
    //     {
    //         return true;
    //     }
    // }

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