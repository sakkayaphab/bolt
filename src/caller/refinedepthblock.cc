#include "refinedepthblock.h"
#include <dirent.h>
RefineDepthBlock::RefineDepthBlock()
{
}

int32_t RefineDepthBlock::roundNumber(int32_t number, int32_t round)
{
    if (number % round == 0)
    {
        // std::cout << number % round << " is even " << std::endl;
        return number;
    }

    int32_t tempNumber = 0;
    if (number % round >= round / 2)
    {
        int32_t numberrounddiff = round - (number % round);
        tempNumber = numberrounddiff + number;
    }
    else
    {
        int32_t numberrounddiff = (number % round);
        tempNumber = number - numberrounddiff;
    }

    return tempNumber;
}

int32_t RefineDepthBlock::nextNumber(int32_t number, int32_t round)
{
    return roundNumber(number, round) + round;
}

int32_t RefineDepthBlock::previousNumber(int32_t number, int32_t round)
{
    return roundNumber(number, round) - round;
}

void RefineDepthBlock::execute()
{
    std::vector<std::string> vcffilelist = getPathVCFFiles();

    for (auto n : vcffilelist)
    {
        auto variantlist = getEvidenceByFilepath(n);

        if (variantlist.size() == 0)
        {
            continue;
        }

        std::vector<Evidence> result;
        if (variantlist.at(0).getSVType() == "DEL")
        {
            // result = getRefineResultDeletion(&variantlist);

            // result = getResultWithOutOverlapped(&variantlist, &variantlist);
            result = getRefineResultDeletion(&variantlist);
            result = getResultRemoveOverlapped(&result, &result);

            // } else if (variantlist.at(0).getSVType()=="DUP") {
            //     // result = getRefineResultDuplication(&result);
        }
        else if (variantlist.at(0).getSVType() == "INV")
        {

            // result = getResultWithOutOverlapped(&variantlist, &variantlist);
            // result = getResultRemoveOverlapped(&variantlist, &variantlist);
            result = getRefineResultInversion(&variantlist);
            result = getResultRemoveOverlapped(&result, &result);

            // result = variantlist;
        }
        else if (variantlist.at(0).getSVType() == "DUP")
        {

            // result = getResultWithOutOverlapped(&variantlist, &variantlist);
            result = getRefineResultDuplication(&variantlist);
            result = getResultRemoveOverlapped(&result, &result);

            // result = variantlist;
        }
        else if (variantlist.at(0).getSVType() == "BND")
        {

            // result = getResultWithOutOverlapped(&variantlist, &variantlist);
            result = getResultRemoveOverlapped(&variantlist, &variantlist);
            result = getResultRemoveOverlapped(&result, &result);

            // result = getRefineResultTranslocation(&variantlist);
            // result = variantlist;
        }
        else if (variantlist.at(0).getSVType() == "INS")
        {

            // result = getResultWithOutOverlapped(&variantlist, &variantlist);
            result = getRefineResultInsertion(&variantlist);
            result = getResultRemoveOverlapped(&result, &result);

            // result = variantlist;
        }
        else
        {
            result = variantlist;
        }
        // if (variantlist.)

        // auto result = getRefineResultDeletion(&variantlist);
        // auto result = getResultWithOutOverlapped(&variantlist, &variantlist);

        writeFile(&result);
    }

    // DepthBlockFile rdf;
    // // rdf.loadDataToCache(filemanager->getReadDepthPath() + "/" + e.getChr() + ".txt")
    // for (int i = 0; i < 40; i++)
    // {
    //          rdf.loadDataToCache(filemanager->getReadDepthPath() + "/chr1.txt");
    //     auto rd = rdf.getBlock(240000000);
    //     // auto rd = rdf.findBlockWithFile(240000000, filemanager->getReadDepthPath() + "/chr1.txt", 250);
    //     std::cout << rd.pos << std::endl;
    // }

    // std::cout << roundNumber(1124, 250) << std::endl;
    // std::cout << nextNumber(1124, 250) << std::endl;
    // std::cout << previousNumber(1124, 250) << std::endl;
}

void RefineDepthBlock::setSampleStat(SampleStat *samplestat)
{
    RefineDepthBlock::samplestat = samplestat;
}

std::vector<Evidence> RefineDepthBlock::getRefineResultTranslocation(std::vector<Evidence> *master)
{
    std::vector<Evidence> cache;
    for (auto n : *master)
    {
        rdf.loadDataToCache(filemanager->getReadDepthPath() + "/" + n.getChr() + ".txt");

        auto currentPos = roundNumber(n.getPos(), roundConfig);
        auto nextPos = nextNumber(n.getPos(), roundConfig);
        auto previousPos = previousNumber(n.getPos(), roundConfig);

        auto currentRD = rdf.getBlock(currentPos);
        auto nextRD = rdf.getBlock(nextPos);
        auto previousRD = rdf.getBlock(previousPos);

        if (n.getMark() == "")
        {
            if (n.LNGMATCH < getDivider(samplestat->getReadLength(), 1, 4, 1))
            {
                continue;
            }
        }

        if (currentRD.depth > readDepthStat.getReadDepthByChr(n.getChr()) * 2)
        {
            continue;
        }

        if (currentRD.INS1 + currentRD.INS2 >= n.getFrequency())
        {
            continue;
        }

        if (currentRD.INV1 + currentRD.INV2 >= n.getFrequency())
        {
            continue;
        }

        if (currentRD.DUP1 + currentRD.DUP2 >= n.getFrequency())
        {
            continue;
        }

        if (currentRD.DEL1 + currentRD.DEL2 >= n.getFrequency())
        {
            continue;
        }

        auto currentEnd = roundNumber(n.getEnd(), roundConfig);
        auto nextEnd = nextNumber(n.getEnd(), roundConfig);
        auto previousEnd = previousNumber(n.getEnd(), roundConfig);

        auto currentEndRD = rdf.getBlock(currentEnd);
        auto nextEndRD = rdf.getBlock(nextEnd);
        auto previousEndRD = rdf.getBlock(previousEnd);

        if (currentEndRD.INV1 + currentEndRD.INV2 >= n.getFrequency())
        {
            continue;
        }

        if (currentEndRD.DUP1 + currentEndRD.DUP2 >= n.getFrequency())
        {
            continue;
        }

        if (currentEndRD.INS1 + currentEndRD.INS2 >= n.getFrequency())
        {
            continue;
        }

        if (currentEndRD.DEL1 + currentEndRD.DEL2 >= n.getFrequency())
        {
            continue;
        }

        if (n.getFrequency() <= 2)
        {
            continue;
        }

        if (n.getMaxMapQ() < 40)
        {
            continue;
        }

        cache.push_back(n);
    }

    return cache;
}

std::vector<Evidence> RefineDepthBlock::getRefineResultInsertion(std::vector<Evidence> *master)
{

    std::vector<Evidence> cache;
    for (auto n : *master)
    {
        rdf.loadDataToCache(filemanager->getReadDepthPath() + "/" + n.getChr() + ".txt");

        auto currentPos = roundNumber(n.getPos(), roundConfig);
        // // auto nextPos = nextNumber(n.getPos(), roundConfig);
        // // auto previousPos = previousNumber(n.getPos(), roundConfig);

        auto currentRD = rdf.getBlock(currentPos);
        // auto nextRD = rdf.getBlock(nextPos);
        // auto previousRD = rdf.getBlock(previousPos);

        // if (currentRD.depth > readDepthStat.getReadDepthByChr(n.getChr()) * 2)
        // {
        //     continue;
        // }

        // if (n.getFrequency() > readDepthStat.getReadDepthByChr(n.getChr()) * 2)
        // {
        //     continue;
        // }

        if ((n.getMark() == "MATEUNMAPPED"))
        {
            // if (currentRD.depth > readDepthStat.getReadDepthByChr(n.getChr()) * 2)
            // {
            //     continue;
            // }

            //  if (currentRD.DUP1 + currentRD.DUP2 >= 2)
            // {
            //     continue;
            // }

            // if (currentRD.TRA1 + currentRD.TRA2 >= 2)
            // {
            //     continue;
            // }

            // if (currentRD.INV1 + currentRD.INV2 >= 2)
            // {
            //     continue;
            // }

            if (n.getFrequency() > readDepthStat.getReadDepthByChr(n.getChr()) * 2)
            {
                continue;
            }

            if (n.getFrequency() <= getDivider(readDepthStat.getReadDepthByChr(n.getChr()), 10, 100, 4))
            {
                continue;
            }

            if (n.getNumberOfRP() <= getDivider(readDepthStat.getReadDepthByChr(n.getChr()), 10, 100, 4))
            {
                continue;
            }

            if (n.LNGMATCH <= getDivider(samplestat->getReadLength(), 30, 100, 15))
            {
                continue;
            }

            if (n.getMaxMapQ() < 20)
            {
                continue;
            }

            if (n.getMaxRPMapQ() < 40)
            {
                continue;
            }

            if (n.getMaxMapQ() >= 60 || n.getMaxRPMapQ() >= 60)
            {
            }
            else
            {
                continue;
            }

            // continue;
        }
        else if ((n.getMark() == "SINS"))
        {
            if (n.getSvLength() < 50)
            {
                continue;
            }

            // continue;
        }
        else if ((n.getMark() == "SR"))
        {

            // if (n.getMaxMapQ() < 30)
            // {
            //     continue;
            // }

            // if (n.getMaxMapQ() < 60)
            // {
            //     if (n.getNumberOfRP() <= getDivider(samplestat->getReadLength(), 1, 100, 1))
            //     {
            //         continue;
            //     }
            //     // continue;
            // }
            // else
            // {
            //     // if (n.getNumberOfRP() <= getDivider(samplestat->getReadLength(), 20, 100, 15))
            //     // {
            //     //     continue;
            //     // }
            //     // continue;
            // }

            if (n.getNumberOfRP() <= 1)
            {
                continue;
            }

            // continue;
        }
        else if ((n.getMark() == "UNMERGE"))
        {

            continue;
        }
        else
        {
            // if (currentRD.depth > readDepthStat.getReadDepthByChr(n.getChr()) * 3)
            // {
            //     continue;
            // }

            if (n.getFrequency() > readDepthStat.getReadDepthByChr(n.getChr()) * 3)
            {
                continue;
            }

            if (n.getMaxMapQ() < 10)
            {
                continue;
            }

            if (n.getFrequency() <= 2)
            {
                continue;
            }

            if (n.getNumberOfRP() <= 3)
            {
                continue;
            }

            if (n.getFrequency() <= 4)
            {
                if (n.getFrequency() <= getDivider(readDepthStat.getReadDepthByChr(n.getChr()), 5, 100, 2))
                {
                    continue;
                }
            }

            if (n.getNumberOfRP() <= 5)
            {
                if (n.getNumberOfRP() <= getDivider(readDepthStat.getReadDepthByChr(n.getChr()), 7, 100, 3))
                {
                    continue;
                }
            }

            if (n.getMaxMapQ() < 60)
            {
                if (n.LNGMATCH <= getDivider(samplestat->getReadLength(), 15, 100, 15))
                {
                    continue;
                }
                // continue;
            }
            else
            {
                if (n.LNGMATCH <= getDivider(samplestat->getReadLength(), 20, 100, 15))
                {
                    continue;
                }
                // continue;
            }

            if (n.getMaxMapQ() >= 60 || n.getMaxRPMapQ() >= 60)
            {
            }
            else
            {
                continue;
            }

            // continue;
        }

        cache.push_back(n);
    }

    return cache;
}

std::vector<Evidence> RefineDepthBlock::getRefineResultDuplication(std::vector<Evidence> *master)
{
    std::vector<Evidence> cache;
    for (auto n : *master)
    {
        rdf.loadDataToCache(filemanager->getReadDepthPath() + "/" + n.getChr() + ".txt");

        if (n.getSvLength() < 90)
        {
            continue;
        }

        // if (n.countDiffEvidencePos(samplestat->getReadLength() * 2) >= 3)
        // {
        //     continue;
        // }

        // if (n.getPos() == n.getEvidencePos())
        // {
        //     if (n.getAltSASize()-n.countEvidencePosNear(n.getEnd(), samplestat->getReadLength()*2) > 0)
        //     {
        //         continue;
        //     }
        // }
        // else
        // {
        //     if (n.getAltSASize()-n.countEvidencePosNear(n.getPos(), samplestat->getReadLength()*2) > 0)
        //     {
        //         continue;
        //     }
        // }

        if (n.getMark() == "SR")
        {
            if (n.getFrequency() <= 1)
            {
                continue;
            }

            if (n.getMaxMapQ() < 30)
            {
                continue;
            }

            if (n.getMaxRPMapQ() == 0)
            {
                continue;
            }

            cache.push_back(n);
            continue;
        }

        if (n.getMark() == "")
        {
            if (n.LNGMATCH < 15)
            {
                continue;
            }
        }

        if (n.LNGMATCH < getDivider(samplestat->getReadLength(), 1, 5, 1))
        {
            continue;
        }

        if (n.getFrequency() <= getDivider(readDepthStat.getReadDepthByChr(n.getChr()), 5, 100, 1) && n.getFrequency() <= getDivider(readDepthStat.getReadDepthByChr(n.getChr()), 5, 100, 1))
        {
            continue;
        }

        auto currentPos = roundNumber(n.getPos(), roundConfig);
        auto nextPos = nextNumber(n.getPos(), roundConfig);
        auto previousPos = previousNumber(n.getPos(), roundConfig);

        auto currentRD = rdf.getBlock(currentPos);
        auto nextRD = rdf.getBlock(nextPos);
        auto previousRD = rdf.getBlock(previousPos);

        if (currentRD.depth > readDepthStat.getReadDepthByChr(n.getChr()) * 4)
        {
            continue;
        }

        // if (n.getSvLength() > 500)
        // {
        if (currentRD.DEL1 + currentRD.DEL2 >= n.getFrequency())
        {
            continue;
        }

        if (currentRD.INS1 + currentRD.INS2 >= n.getFrequency())
        {
            continue;
        }

        if (currentRD.INV1 + currentRD.INV2 >= n.getFrequency())
        {
            continue;
        }

        if (currentRD.TRA1 + currentRD.TRA2 >= n.getFrequency())
        {
            continue;
        }
        // }

        auto currentEnd = roundNumber(n.getEnd(), roundConfig);
        auto nextEnd = nextNumber(n.getEnd(), roundConfig);
        auto previousEnd = previousNumber(n.getEnd(), roundConfig);

        auto currentEndRD = rdf.getBlock(currentEnd);
        auto nextEndRD = rdf.getBlock(nextEnd);
        auto previousEndRD = rdf.getBlock(previousEnd);

        // if (n.getSvLength() > 500)
        // {
        if (currentEndRD.DEL1 + currentEndRD.DEL2 >= n.getFrequency())
        {
            continue;
        }

        if (currentEndRD.INS1 + currentEndRD.INS2 >= n.getFrequency())
        {
            continue;
        }

        if (currentEndRD.INV1 + currentEndRD.INV2 >= n.getFrequency())
        {
            continue;
        }

        if (currentEndRD.TRA1 + currentEndRD.TRA2 >= n.getFrequency())
        {
            continue;
        }
        // }

        if (n.getFrequency() <= 1)
        {
            continue;
        }

        if (n.getMaxMapQ() < 30)
        {
            continue;
        }

        if (n.getMaxRPMapQ() == 0)
        {
            continue;
        }

        cache.push_back(n);
    }

    return cache;
}

std::vector<Evidence> RefineDepthBlock::getRefineResultInversion(std::vector<Evidence> *master)
{
    std::vector<Evidence> cache;
    for (auto n : *master)
    {
        cache.push_back(n);
        continue;
        rdf.loadDataToCache(filemanager->getReadDepthPath() + "/" + n.getChr() + ".txt");
        auto currentPos = roundNumber(n.getPos(), roundConfig);
        auto nextPos = nextNumber(n.getPos(), roundConfig);
        auto previousPos = previousNumber(n.getPos(), roundConfig);

        auto currentRD = rdf.getBlock(currentPos);
        auto nextRD = rdf.getBlock(nextPos);
        auto previousRD = rdf.getBlock(previousPos);

        if (n.getMark() == "SR")
        {
            // continue;
            if (n.getSvLength() < samplestat->getReadLength())
            {
                continue;
            }
            // continue;
            if (n.getFrequency() <= 2)
            {
                continue;
            }

            // if (n.getMaxMapQ() < 40)
            // {
            //     continue;
            // }

            if (n.getMaxRPMapQ() == 0)
            {
                continue;
            }

            cache.push_back(n);
            continue;
        }

        // continue;

        if (n.getMark() == "")
        {
            if (n.LNGMATCH < 15)
            {
                continue;
            }
        }

        if (n.LNGMATCH < getDivider(samplestat->getReadLength(), 1, 5, 1))
        {
            continue;
        }

        if (currentRD.depth > readDepthStat.getReadDepthByChr(n.getChr()) * 4)
        {
            continue;
        }

        if (currentRD.depth == 0)
        {
            continue;
        }

        if (n.getSvLength() > 500)
        {
            if (currentRD.DEL1 + currentRD.DEL2 >= n.getFrequency())
            {
                continue;
            }

            if (currentRD.INS1 + currentRD.INS2 >= n.getFrequency())
            {
                continue;
            }

            if (currentRD.TRA1 + currentRD.TRA2 >= n.getFrequency())
            {
                continue;
            }

            if (currentRD.DUP1 + currentRD.DUP2 >= n.getFrequency())
            {
                continue;
            }
        }

        auto currentEnd = roundNumber(n.getEnd(), roundConfig);
        auto nextEnd = nextNumber(n.getEnd(), roundConfig);
        auto previousEnd = previousNumber(n.getEnd(), roundConfig);

        auto currentEndRD = rdf.getBlock(currentEnd);
        auto nextEndRD = rdf.getBlock(nextEnd);
        auto previousEndRD = rdf.getBlock(previousEnd);

        if (n.getSvLength() > 500)
        {
            if (currentEndRD.DEL1 + currentEndRD.DEL2 >= n.getFrequency())
            {
                continue;
            }

            if (currentEndRD.TRA1 + currentEndRD.TRA2 >= n.getFrequency())
            {
                continue;
            }

            if (currentEndRD.DUP1 + currentEndRD.DUP2 >= n.getFrequency())
            {
                continue;
            }

            if (currentEndRD.INS1 + currentEndRD.INS2 >= n.getFrequency())
            {
                continue;
            }
        }

        if (n.getFrequency() <= 2)
        {
            continue;
        }

        if (n.getSvLength() > 1000000)
        {
            continue;
        }

        if (n.getSvLength() < 80)
        {
            continue;
        }

        if (n.getMaxMapQ() == 0)
        {
            continue;
        }

        cache.push_back(n);
    }

    return cache;
}

std::vector<Evidence> RefineDepthBlock::getRefineResultDeletion(std::vector<Evidence> *master)
{
    std::vector<Evidence> cache;
    int count = 0;
    for (auto n : *master)
    {
        rdf.loadDataToCache(filemanager->getReadDepthPath() + "/" + n.getChr() + ".txt");

        auto currentPos = roundNumber(n.getPos(), roundConfig);
        auto nextPos = nextNumber(n.getPos(), roundConfig);
        auto previousPos = previousNumber(n.getPos(), roundConfig);

        auto currentRD = rdf.getBlock(currentPos);
        auto nextRD = rdf.getBlock(nextPos);
        auto previousRD = rdf.getBlock(previousPos);

        auto currentEnd = roundNumber(n.getEnd(), roundConfig);
        auto nextEnd = nextNumber(n.getEnd(), roundConfig);
        auto previousEnd = previousNumber(n.getEnd(), roundConfig);

        auto currentEndRD = rdf.getBlock(currentEnd);
        auto nextEndRD = rdf.getBlock(nextEnd);
        auto previousEndRD = rdf.getBlock(previousEnd);

        if (n.getSvLength() > 1000000)
        {
            continue;
        }

        if (currentRD.depth > readDepthStat.getReadDepthByChr(n.getChr()) * 3 && currentEndRD.depth > readDepthStat.getReadDepthByChr(n.getChr()) * 3)
        {
            continue;
        }

        if (n.getMark() == "SR")
        {

            // continue;
            if (n.getFrequency() <= 1)
            {
                continue;
            }

            if (n.getMaxRPMapQ() < 60)
            {
                continue;
            }

            if (n.getFrequency() <= getDivider(readDepthStat.getReadDepthByChr(n.getChr()), 5, 100, 1) && n.getFrequency() <= getDivider(readDepthStat.getReadDepthByChr(n.getChr()), 5, 100, 1))
            {
                continue;
            }

            if (n.getNumberOfRP() <= getDivider(readDepthStat.getReadDepthByChr(n.getChr()), 5, 100, 1) && n.getNumberOfRP() <= getDivider(readDepthStat.getReadDepthByChr(n.getChr()), 5, 100, 1))
            {
                continue;
            }

            if (n.getSvLength() < 20)
            {
                continue;
            }

            cache.push_back(n);
            continue;
        }
        else if (n.getMark() == "SDEL")
        {

            if (currentRD.INS1 + currentRD.INS2 >= n.getFrequency())
            {
                continue;
            }

            if (currentRD.INV1 + currentRD.INV2 >= n.getFrequency())
            {
                continue;
            }

            if (currentRD.TRA1 + currentRD.TRA2 >= n.getFrequency())
            {
                continue;
            }

            if (currentRD.DUP1 + currentRD.DUP2 >= n.getFrequency())
            {
                continue;
            }

            // continue;
            if (n.getFrequency() <= 5)
            {
                continue;
            }

            if (n.getMaxMapQ() < 60)
            {
                continue;
            }

            if (n.getFrequency() <= getDivider(readDepthStat.getReadDepthByChr(n.getChr()), 5, 100, 1) && n.getFrequency() <= getDivider(readDepthStat.getReadDepthByChr(n.getChr()), 5, 100, 1))
            {
                continue;
            }

            if (n.getSvLength() < 20)
            {
                continue;
            }

            cache.push_back(n);
            continue;
        }
        else
        {
            if (n.LNGMATCH < 20)
            {
                continue;
            }
        }

        // continue;

        if (n.LNGMATCH < getDivider(samplestat->getReadLength(), 1, 5, 1))
        {
            continue;
        }

        // continue;

        if (currentRD.depth > readDepthStat.getReadDepthByChr(n.getChr()) && currentEndRD.depth > readDepthStat.getReadDepthByChr(n.getChr()))
        {
            continue;
        }

        if (n.getSvLength() > 2000)
        {
            if (currentRD.INS1 + currentRD.INS2 >= n.getFrequency())
            {
                continue;
            }

            if (currentRD.INV1 + currentRD.INV2 >= n.getFrequency())
            {
                continue;
            }

            if (currentRD.TRA1 + currentRD.TRA2 >= n.getFrequency())
            {
                continue;
            }

            if (currentRD.DUP1 + currentRD.DUP2 >= n.getFrequency())
            {
                continue;
            }
        }

        if (n.getSvLength() > 2000)
        {
            if (currentEndRD.INS1 + currentEndRD.INS2 >= n.getFrequency())
            {
                continue;
            }

            if (currentEndRD.INV1 + currentEndRD.INV2 >= n.getFrequency())
            {
                continue;
            }

            if (currentEndRD.TRA1 + currentEndRD.TRA2 >= n.getFrequency())
            {
                continue;
            }

            if (currentEndRD.DUP1 + currentEndRD.DUP2 >= n.getFrequency())
            {
                continue;
            }
        }

        if (n.getMaxMapQ() < 30)
        {
            continue;
        }

        if (n.getMaxRPMapQ() < 30)
        {
            continue;
        }

        if (n.getSvLength() > 1000000)
        {
            continue;
        }

        if (n.getFrequency() <= 1)
        {
            continue;
        }

        if (n.getFrequency() <= getDivider(readDepthStat.getReadDepthByChr(n.getChr()), 5, 100, 1) && n.getFrequency() <= getDivider(readDepthStat.getReadDepthByChr(n.getChr()), 5, 100, 1))
        {
            continue;
        }

        if (n.getSvLength() < 500)
        {
            if (n.getFrequency() <= getDivider(readDepthStat.getReadDepthByChr(n.getChr()), 15, 100, 1) && n.getFrequency() <= getDivider(readDepthStat.getReadDepthByChr(n.getChr()), 15, 100, 1))
            {
                continue;
            }
        }

        if (n.getSvLength() < 20)
        {
            continue;
        }

        // if (n.getPos() == n.getEvidencePos())
        // {
        //     if (n.getAltSASize()-n.countEvidencePosNear(n.getEnd(), samplestat->getReadLength()) > 0)
        //     {
        //         continue;
        //     }
        // }
        // else
        // {
        //     if (n.getAltSASize()-n.countEvidencePosNear(n.getPos(), samplestat->getReadLength()) > 0)
        //     {
        //         continue;
        //     }
        // }

        cache.push_back(n);
    }

    return cache;
}

void RefineDepthBlock::writeFile(std::vector<Evidence> *master)
{
    std::ofstream myfile;
    myfile.open(filemanager->getOutputPath() + "/result.vcf", std::ios_base::app);
    for (auto n : *master)
    {
        n.setID("BOLT" + std::to_string(vcfIdNumber));
        myfile << n.getResultVcfFormatString() << std::endl;
        vcfIdNumber++;
    }

    myfile.close();
}

std::vector<Evidence> RefineDepthBlock::getResultWithOutOverlapped(std::vector<Evidence> *master, std::vector<Evidence> *slave)
{
    std::vector<Evidence> cache;
    for (auto n : *master)
    {
        bool found = false;
        for (auto m : *slave)
        {
            if (n.getPos() == m.getPos() && n.getEnd() == m.getEnd())
            {
                continue;
            }

            if (m.getPos() - 100 <= n.getPos() && n.getPos() <= m.getPos() + 100)
            {
                found = true;
                break;
            }
        }

        if (!found)
        {
            cache.push_back(n);
        }
    }

    return cache;
}

std::vector<Evidence> RefineDepthBlock::getResultRemoveOverlapped(std::vector<Evidence> *master, std::vector<Evidence> *slave)
{
    std::vector<Evidence> cache;
    for (auto n : *master)
    {
        bool found = false;
        for (auto m : cache)
        {
            if (n.getPos() == m.getPos() && n.getEnd() == m.getEnd())
            {
                found = true;
                break;
            }

            if (checkBetween(n.getPos(), m.getPos(), 10) && checkBetween(n.getEnd(), m.getEnd(), 10))
            {
                // std::cout << n.getChr() << " " << n.getPos() << " " << m.getEnd() << std::endl;
                found = true;
                break;
            }
        }

        if (!found)
        {
            cache.push_back(n);
        }
    }

    return cache;
}

bool RefineDepthBlock::checkBetween(int32_t pos, int32_t targetPos, int32_t overlapped)
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

std::vector<Evidence> RefineDepthBlock::getEvidenceByFilepath(std::string filepaht)
{
    std::vector<Evidence> cache;
    std::string line;
    std::ifstream myfile(filepaht);
    if (myfile.is_open())
    {
        // std::cout << "load evidence :" << filepaht << std::endl;
        while (getline(myfile, line))
        {
            Evidence e;
            e.setEvidenceByString(line);
            if (e.getPos() != 0)
            {

                // if (e.getMark() == "SR")
                // {
                //     cache.push_back(e);
                // }
                // else
                // {
                //     if (e.LNGMATCH < getDivider(samplestat->getReadLength(), 1, 4, 1))
                //     {
                //         continue;
                //     }

                //     if (e.getFrequency() <= 1)
                //     {
                //         continue;
                //     }

                cache.push_back(e);
                // cache.push_back(n);
                // continue;
                // }

                // if (e.getSvLength()<500 && e.getAvgMapQ()<6) {

                // } else {

                // }
            }

            // std::cout
            //     << e.getChr() << "\t"
            //     << e.getPos() << "\t"
            //     << e.getEndChr() << "\t"
            //     << e.getEnd() << "\t"
            //     << e.getFrequency() << "\t"
            //     << e.LNGMATCH << "\t"
            //     << e.convertMapQlistToCommaString() << "\t"
            //     << std::endl;
        }
        myfile.close();
        // std::cout << "End + load evidence :" << filepaht << std::endl;
    }

    return cache;
}

void RefineDepthBlock::setFileManager(FileManager *filemanager)
{
    RefineDepthBlock::filemanager = filemanager;
    readDepthStat.setFilePath(filemanager);
    readDepthStat.execute();
}

std::vector<std::string> RefineDepthBlock::getPathVCFFiles()
{
    std::vector<std::string> evidenceFilePathLists;
    DIR *d;
    struct dirent *dir;

    d = opendir(filemanager->getVariantPath().c_str());
    if (d)
    {
        while (dir = readdir(d))
        {
            if (std::string(dir->d_name).size() < 4)
            {
                continue;
            }

            if (std::string(dir->d_name).substr(std::string(dir->d_name).size() - 4) == ".vcf")
            {
                std::string tempPath = filemanager->getVariantPath() + std::string(dir->d_name);
                // std::cout << "tempPath : " << tempPath << std::endl;
                evidenceFilePathLists.push_back(tempPath);
            }
        }
        closedir(d);
    }

    return evidenceFilePathLists;
}

int RefineDepthBlock::getDivider(int value, int top, int down, int minimum)
{
    auto returnvalue = (int)(float(value) * (float(top) / float(down)));

    if (returnvalue > minimum)
    {
        return returnvalue;
    }

    return minimum;
}