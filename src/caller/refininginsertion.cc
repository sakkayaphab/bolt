#include "refininginsertion.h"

RefiningInsertion::RefiningInsertion()
{
    variantresult.setVariantType("INS");
}

void RefiningInsertion::execute()
{
    variantresult.setChr(evidence.getChr());
    variantresult.setEndChr(evidence.getEndChr());
    prepareBamReader();
    first();
}

void RefiningInsertion::first()
{
    std::string findRange = convertRangeToString(evidence.getChr(), evidence.getPos() + evidence.getCiPosLeft(),
                                                 evidence.getPos() + evidence.getCiPosRight());

    std::cout << findRange << std::endl;

    const char *range = findRange.c_str();
    refineStartToEnd(range);
}

void RefiningInsertion::refineStartToEnd(const char *range)
{
    hts_itr_t *iter = NULL;

    iter = sam_itr_querys(bam_index, bam_header, range);
    if (iter == NULL)
        return;
    read = bam_init1();
    readparser.setBamHeader(bam_header);
    readparser.setBamRead(read);

    while (sam_itr_next(inFile, iter, read) >= 0)
    {

        if (readparser.isUnmapped())
        {
            continue;
        }

        if (readparser.isNotPassingFilters())
        {
            continue;
        }

        if (readparser.isPCR())
        {
            continue;
        }

        if (readparser.isSupplementaryAlignment())
        {
            continue;
        }

        auto cigar = readparser.getCigar();
        if (cigar.size() <= 1)
        {
            continue;
        }

        if (cigar.at(cigar.size() - 1).getOperatorName() == 'S' && cigar.at(cigar.size() - 1).getLength() >= 4)
        {
            mapSCEnd[readparser.getEnd()].addMapQ(readparser.getMapQuality());
            mapSCEnd[readparser.getEnd()].addLongMapping(cigar.at(cigar.size() - 1).getLength());
            mapSCEnd[readparser.getEnd()].setPosition(readparser.getEnd());
        }

        if (cigar.at(0).getOperatorName() == 'S' && cigar.at(0).getLength() >= 4)
        {
            // std::cout << readparser.getPos() << std::endl;
            mapSCStart[readparser.getPos()].addMapQ(readparser.getMapQuality());
            mapSCStart[readparser.getPos()].addLongMapping(cigar.at(0).getLength());
            mapSCStart[readparser.getPos()].setPosition(readparser.getPos());
        }
    }

    RefiningInsertion::convertMapSC();
    RefiningInsertion::clearMapSC();
    RefiningInsertion::findBreakpoint();
    RefiningInsertion::filterBreakpoint();

    hts_itr_destroy(iter);

    return;
}

void RefiningInsertion::filterBreakpoint()
{
    std::sort(vectorBP.begin(), vectorBP.end());

    for (BreakpointPosition n : vectorBP)
    {
        int32_t averagePos = 0;
        if (n.pos > n.end)
        {
            averagePos = n.pos;
        }
        else
        {
            averagePos = (n.end + n.pos) / 2;
        }

        // std::cout << n.pos << " = " << n.end << std::endl;

        variantresult.setPos(averagePos);
        variantresult.setEnd(averagePos);
        variantresult.setFrequency(n.mappingqualitylist.size());
        variantresult.setRPMapQ(*evidence.getMapQVector());

        variantresult.setMapQList(n.mappingqualitylist);
        variantresult.setChr(evidence.getChr());
        variantresult.setEndChr(evidence.getEndChr());
        variantresult.setQuailtyPass(true);
        break;
    }
}

void RefiningInsertion::findBreakpoint()
{
    for (InsertionPositionDetail n : vectorSCStart)
    {
        if (n.getLongMapping() < 4)
        {
            continue;
        }

        if (n.getFrequency() < 1)
        {
            continue;
        }

        bool added;

        for (InsertionPositionDetail m : vectorSCEnd)
        {
            if (m.getLongMapping() < 4)
            {
                continue;
            }

            if (m.getFrequency() < 1)
            {
                continue;
            }

            if (checkBetween(n.getPosition(), m.getPosition(), samplestat->getReadLength()))
            {
                BreakpointPosition tempBP;
                tempBP.pos = n.getPosition();
                tempBP.end = m.getPosition();
                tempBP.frequency = n.getFrequency() + m.getFrequency();
                tempBP.score = n.getFrequency() + m.getFrequency();
                tempBP.longmapstart = n.getLongMapping();
                tempBP.longmapend = m.getLongMapping();
                for (auto x : n.getMapQList())
                {
                    tempBP.mappingqualitylist.push_back(x);
                }
                for (auto x : m.getMapQList())
                {
                    tempBP.mappingqualitylist.push_back(x);
                }

                if (tempBP.frequency <= 3)
                {
                    continue;
                }
                added = true;
                vectorBP.push_back(tempBP);
            }
        }

        if (!added && n.getFrequency() >= 4)
        {
            BreakpointPosition tempBP;
            tempBP.pos = n.getPosition();
            tempBP.end = n.getPosition();
            tempBP.frequency = n.getFrequency();
            tempBP.score = n.getFrequency();
            tempBP.longmapstart = n.getLongMapping();
            tempBP.longmapend = n.getLongMapping();
            for (auto x : n.getMapQList())
            {
                tempBP.mappingqualitylist.push_back(x);
            } 
 
            added = true;
            vectorBP.push_back(tempBP);
        }
    }
}

bool RefiningInsertion::checkBetween(int32_t pos, int32_t targetPos, int32_t overlapped)
{
    if (targetPos - overlapped > pos)
    {
        return false;
    }

    if (targetPos + overlapped < pos)
    {
        return false;
    }

    return true;
}

void RefiningInsertion::convertMapSC()
{
    convertMapSCToVector(&mapSCStart, &vectorSCStart);
    std::sort(vectorSCStart.begin(), vectorSCStart.end());

    convertMapSCToVector(&mapSCEnd, &vectorSCEnd);
    std::sort(vectorSCEnd.begin(), vectorSCEnd.end());
}

void RefiningInsertion::convertMapSCToVector(std::map<int32_t, InsertionPositionDetail> *mapSC, std::vector<InsertionPositionDetail> *vectorSC)
{
    for (auto const &x : *mapSC)
    {
        vectorSC->push_back(x.second);
    }
}

void RefiningInsertion::clearMapSC()
{
    mapSCStart.clear();
    mapSCEnd.clear();
}