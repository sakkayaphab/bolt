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

    //large insertion by unmapped read
    //    std::cout << evidence.getPosDiscordantRead() << " // " << evidence.getEndDiscordantRead() << std::endl;
    // if (evidence.getPosDiscordantRead() == evidence.getEndDiscordantRead())
    // {
        first();
    // }
    // small - medium insertion
    // else
    // {
    //     /* code */
    // }

    //    std::cout << "---+ run complete +---" << std::endl;
}

void RefiningInsertion::first()
{
    std::string findRange = convertRangeToString(evidence.getChr(), evidence.getPos() + evidence.getCiPosLeft(),
                                                 evidence.getPos() + evidence.getCiPosRight());

    const char *range = findRange.c_str();
    // const char *mChr = evidence.getChr().c_str();

    // std::cout << range << "/" << samplestat->getReadLength() << std::endl;
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

        auto cigar = readparser.getCigar();
        if (cigar.size() <= 1)
        {
            continue;
        }

       
        if (cigar.at(0).getOperatorName() == 'S' && cigar.at(0).getLength() >= 5)
        {
            

            mapSCFirst[readparser.getPos()]++;
            if (mapMapQFirst[readparser.getPos()]<readparser.getMapQuality()) {
                mapMapQFirst[readparser.getPos()] = readparser.getMapQuality();
            }
        }

        if (cigar.at(cigar.size() - 1).getOperatorName() == 'S' && cigar.at(cigar.size() - 1).getLength() >= 5)
        {
          
            mapSCLast[readparser.getEnd()]++;
            if (mapMapQLast[readparser.getEnd()]<readparser.getMapQuality()) {
                mapMapQLast[readparser.getEnd()] = readparser.getMapQuality();
            }
        }
    }

    refineVariant(range);

    hts_itr_destroy(iter);

    return;
}

void RefiningInsertion::refineVariant(const char *range) {
    //find max
    
    int32_t position_first_hit = 0;
    int32_t position_second_hit = 0;
    int hit_position_second = 0;
    int hit_position_first = 0;
    if (getPosMaxHitValue(&mapSCLast,&mapMapQLast) < getPosMaxHitValue(&mapSCFirst,&mapMapQFirst))
    {
        position_first_hit = getPosMaxHitValue(&mapSCFirst,&mapMapQFirst);
     

        for (auto const &x : mapSCLast)
        {

            if (!isBetWeen(position_first_hit,x.first,100)) {
                continue;
            }

            if (x.second <= 3)
            {
                continue;
            }

            if (hit_position_second < x.second)
            {
                position_second_hit = x.first;
                hit_position_second = x.second;
            }
        }
    }
    else
    {
        position_second_hit = getPosMaxHitValue(&mapSCLast,&mapMapQLast);

        for (auto const &x : mapSCFirst)
        {

            if (!isBetWeen(position_second_hit,x.first,100)) {
                continue;
            }

            if (x.second <= 3)
            {
                continue;
            }

            if (hit_position_first < x.second)
            {
                position_first_hit = x.first;
                hit_position_first = x.second;
            }
        }
    }

    if (position_second_hit == 0 && position_first_hit == 0)
    {
        return;
    }
    
    if (position_second_hit == 0 )
    {
        position_second_hit = position_first_hit;
    }
    if (position_first_hit == 0)
    {
        position_second_hit = position_first_hit;
    }

    int32_t averagePos = (position_first_hit + position_second_hit) / 2;
    int readdepthAtPos = getReadDepthAtPosition(range, averagePos);

    if (readdepthAtPos < 5)
    {
        return;
    }

    if (readdepthAtPos > 500)
    {
        return;
    }

    auto ratio = (hit_position_second + hit_position_first) / (float) readdepthAtPos;
    if (ratio < 0.05)
    {
        return;
    }

    variantresult.setPos(averagePos);
    variantresult.setEnd(averagePos);
    variantresult.setFrequency(hit_position_second + hit_position_first);
    variantresult.setQuailtyPass(true);
}

bool RefiningInsertion::isBetWeen(int32_t primary,int32_t secondary,int32_t range) {
    if (secondary < primary-range) {
        return false;
    }

    if (secondary > primary+range) {
        return false;
    }

    return true;
}

int RefiningInsertion::getHitByPos(std::map<int32_t, int> *map, int32_t pos)
{
    return (*map)[pos];
}

int32_t RefiningInsertion::getPosMaxHitValue(std::map<int32_t, int> *map,std::map<int32_t, uint8_t> *mapMapQ)
{
    int hit = 0;
    int32_t position = 0;

    for (auto const &x : *map)
    {
        // if ((*mapMapQ)[x.first]<25) {
        //     continue;
        // }

        if ((*mapMapQ)[x.first]==0) {
            continue;
        }

        if (x.second > hit)
        {
            hit = x.second;
            position = x.first;
        }
    }

    return position;
}
