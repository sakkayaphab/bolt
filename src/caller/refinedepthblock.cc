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
            result = getResultWithOutOverlapped(&variantlist, &variantlist);
            result = getRefineResultDeletion(&result);
            // } else if (variantlist.at(0).getSVType()=="DUP") {
            //     // result = getRefineResultDuplication(&result);
        }
        else if (variantlist.at(0).getSVType() == "INV")
        {

            // result = getResultWithOutOverlapped(&variantlist, &variantlist);
            result = getRefineResultInversion(&variantlist);
            // result = variantlist;
        }
        else if (variantlist.at(0).getSVType() == "DUP")
        {

            // result = getResultWithOutOverlapped(&variantlist, &variantlist);
            result = getRefineResultDuplication(&variantlist);
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

std::vector<Evidence> RefineDepthBlock::getRefineResultDuplication(std::vector<Evidence> *master)
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
            if (n.LNGMATCH < ((25 / 100) * samplestat->getReadLength()))
            {
                continue;
            }
        }

        if (currentRD.depth > readDepthStat.getReadDepthByChr(n.getChr()) * 4)
        {
            continue;
        }

        if (n.getSvLength() > 500)
        {

            //     if (currentRD.DEL1 + currentRD.DEL2 >= 10)
            //     {
            //         continue;
            //     }

            //     if (previousRD.DEL1 + previousRD.DEL2 >= 10)
            //     {
            //         continue;
            //     }

            //     if (nextRD.DEL1 + previousRD.DEL2 >= 10)
            //     {
            //         continue;
            //     }

            if (currentRD.INV1 + currentRD.INV2 >= readDepthStat.getReadDepthByChr(n.getChr()))
            {
                continue;
            }

            if (previousRD.INV1 + previousRD.INV2 >= readDepthStat.getReadDepthByChr(n.getChr()))
            {
                continue;
            }

            if (nextRD.INV1 + previousRD.INV2 >= readDepthStat.getReadDepthByChr(n.getChr()))
            {
                continue;
            }

            if (currentRD.TRA1 + currentRD.TRA2 >= readDepthStat.getReadDepthByChr(n.getChr()))
            {
                continue;
            }

            if (previousRD.TRA1 + previousRD.TRA2 >= readDepthStat.getReadDepthByChr(n.getChr()))
            {
                continue;
            }

            if (nextRD.TRA1 + previousRD.TRA2 >= readDepthStat.getReadDepthByChr(n.getChr()))
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
            //     if (currentEndRD.DEL1 + currentEndRD.DEL2 >= 10)
            //     {
            //         continue;
            //     }

            //     if (previousEndRD.DEL1 + previousEndRD.DEL2 >= 10)
            //     {
            //         continue;
            //     }

            //     if (nextEndRD.DEL1 + previousEndRD.DEL2 >= 10)
            //     {
            //         continue;
            //     }

            if (currentEndRD.INV1 + currentEndRD.INV2 >= readDepthStat.getReadDepthByChr(n.getChr()))
            {
                continue;
            }

            if (previousEndRD.INV1 + previousEndRD.INV2 >= readDepthStat.getReadDepthByChr(n.getChr()))
            {
                continue;
            }

            if (nextEndRD.INV1 + nextEndRD.INV2 >= readDepthStat.getReadDepthByChr(n.getChr()))
            {
                continue;
            }

            if (currentEndRD.TRA1 + currentEndRD.TRA2 >= readDepthStat.getReadDepthByChr(n.getChr()))
            {
                continue;
            }

            if (previousEndRD.TRA1 + previousEndRD.TRA2 >= readDepthStat.getReadDepthByChr(n.getChr()))
            {
                continue;
            }

            if (nextEndRD.TRA1 + nextEndRD.TRA2 >= readDepthStat.getReadDepthByChr(n.getChr()))
            {
                continue;
            }
        }

        if (n.getSvLength() < 100)
        {
            continue;
        }

        if (n.getFrequency() <= 1)
        {
            continue;
        }

        // // if (n.getFrequency() > 5) {
        // //     continue;
        // // }

        if (n.getMaxMapQ() < 50)
        {
            continue;
        }

        //  if (n.getMaxMapQ() < 15) {
        //      continue;
        //  }

        // if (n.getMaxMapQ() < 60)
        // {
        //     if (n.getFrequency() <= 2)
        //     {
        //         continue;
        //     }
        // }

        // if (n.getFrequency()>=30) {
        //     continue;
        // }

        // if (n.getMaxMapQ()==0) {
        //     continue;
        // }

        // if (n.getAvgMapQ()<2) {
        //     continue;
        // }

        // if (n.getSvLength()>10000000) {
        //     continue;
        // }

        cache.push_back(n);
    }

    return cache;
}

std::vector<Evidence> RefineDepthBlock::getRefineResultInversion(std::vector<Evidence> *master)
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

      

        if (n.getSvLength() > 2000)
        {
            // if (currentRD.DUP1 + currentRD.DUP2 >= 20)
            // {
            //     continue;
            // }

            // if (previousRD.DUP1 + previousRD.DUP2 >= 20)
            // {
            //     continue;
            // }

            // if (nextRD.DUP1 + previousRD.DUP2 >= 20)
            // {
            //     continue;
            // }

            if (currentRD.TRA1 + currentRD.TRA2 >= readDepthStat.getReadDepthByChr(n.getChr()))
            {
                continue;
            }

            // if (previousRD.TRA1 + previousRD.TRA2 >= readDepthStat.getReadDepthByChr(n.getChr()))
            // {
            //     continue;
            // }

            // if (nextRD.TRA1 + previousRD.TRA2 >= readDepthStat.getReadDepthByChr(n.getChr()))
            // {
            //     continue;
            // }
        }

        auto currentEnd = roundNumber(n.getEnd(), roundConfig);
        auto nextEnd = nextNumber(n.getEnd(), roundConfig);
        auto previousEnd = previousNumber(n.getEnd(), roundConfig);

        auto currentEndRD = rdf.getBlock(currentEnd);
        auto nextEndRD = rdf.getBlock(nextEnd);
        auto previousEndRD = rdf.getBlock(previousEnd);

        

        if (n.getSvLength() > 2000)
        {
            // if (currentEndRD.DUP1 + currentEndRD.DUP2 >= 20)
            // {
            //     continue;
            // }

            // if (previousEndRD.DUP1 + previousEndRD.DUP2 >= 20)
            // {
            //     continue;
            // }

            // if (nextEndRD.DUP1 + previousEndRD.DUP2 >= 20)
            // {
            //     continue;
            // }

            if (currentEndRD.TRA1 + currentEndRD.TRA2 >= readDepthStat.getReadDepthByChr(n.getChr()))
            {
                continue;
            }

            // if (previousEndRD.TRA1 + previousEndRD.TRA2 >= readDepthStat.getReadDepthByChr(n.getChr()))
            // {
            //     continue;
            // }

            // if (nextEndRD.TRA1 + nextEndRD.TRA2 >= readDepthStat.getReadDepthByChr(n.getChr()))
            // {
            //     continue;
            // }
        }

        // if (n.getSvLength() < 20)
        // {
        //     continue;
        // }

        if (n.getFrequency() <= 1)
        {
            continue;
        }

        if (n.getFrequency() >= readDepthStat.getReadDepthByChr(n.getChr()))
        {
            continue;
        }

        // if (n.getMaxMapQ() < 50)
        // {
        //     continue;
        // }

        // if (n.getAvgMapQ()<2) {
        //     continue;
        // }

        if (n.getSvLength() > 1000000)
        {
            continue;
        }

        if (n.getSvLength() < 100)
        {
            continue;
        }

        // if (n.getFrequency()<=2) {
        //         continue;
        // }

        // if (n.getFrequency()==1) {
        //     if (n.getMaxMapQ()<25) {
        //         continue;
        //     }
        // }

        // if (n.getSvLength()<300) {
        //      if (n.getMaxMapQ()<10) {
        //         continue;
        //     }
        // }

        // if (n.getSvLength() < 1000)
        // {
        //     // if (n.getFrequency()<=3) {
        //     //     continue;
        //     // }
        //     // if (currentRD.depth>10 || currentEndRD.depth>10) {
        //     //     continue;
        //     // }

        //     // if (currentEndRD.depth>100) {
        //     //     continue;
        //     // }

        //     // if (n.getMaxMapQ()<=25) {
        //     //     continue;
        //     // }
        // }

        cache.push_back(n);
    }

    return cache;
}

std::vector<Evidence> RefineDepthBlock::getRefineResultDeletion(std::vector<Evidence> *master)
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

        // if (currentRD.depth==0) {
        //     continue;
        // }

        

        

        std::cout << getDivider(10,2,10,1) << " = "<< 10*(2/10) << std::endl;

        // if (previousRD.TRA1+previousRD.TRA2 >= getDivider(readDepthStat.getReadDepthByChr(n.getChr()),1,minimumdivide,1))
        // {
        //     continue;
        // }

        // if (previousRD.INV1+previousRD.INV2 >= getDivider(readDepthStat.getReadDepthByChr(n.getChr()),1,minimumdivide,1))
        // {
        //     continue;
        // }

        if (currentRD.INS1+currentRD.INS2 >= getDivider(readDepthStat.getReadDepthByChr(n.getChr()),1,minimumdivide,1))
        {
            continue;
        }

        if (currentRD.INV1+currentRD.INV2 >= getDivider(readDepthStat.getReadDepthByChr(n.getChr()),1,minimumdivide,1))
        {
            continue;
        }

        if (currentRD.TRA1+currentRD.TRA2 >= getDivider(readDepthStat.getReadDepthByChr(n.getChr()),1,minimumdivide,1))
        {
           continue;
        }


        if (currentRD.INV1+currentRD.INV2 >= n.getFrequency())
        {
            continue;
        }

        if (currentRD.TRA1+currentRD.TRA2 >= n.getFrequency())
        {
           continue;
        }

        if (nextRD.INV1+nextRD.INV2 >= n.getFrequency())
        {
            continue;
        }

        if (nextRD.TRA1+nextRD.TRA2 >= n.getFrequency())
        {
           continue;
        }


        if (n.getSvLength() > 2000)
        {
            if (currentRD.DUP1 + currentRD.DUP2 >= 10)
            {
                continue;
            }

            if (currentRD.INV1 + currentRD.INV2 >= 10)
            {
                continue;
            }

            if (previousRD.DUP1 + previousRD.DUP2 >= 10)
            {
                continue;
            }

            if (previousRD.INV1 + previousRD.INV2 >= 10)
            {
                continue;
            }

            if (nextRD.DUP1 + previousRD.DUP2 >= 10)
            {
                continue;
            }

            if (nextRD.INV1 + previousRD.INV2 >= 10)
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

        // if (currentEndRD.depth==0) {
        //     continue;
        // }

        if (currentEndRD.INS1+currentEndRD.INS2 >= getDivider(readDepthStat.getReadDepthByChr(n.getChr()),1,minimumdivide,1))
        {
            continue;
        }

        if (currentEndRD.INV1+currentEndRD.INV2 >= getDivider(readDepthStat.getReadDepthByChr(n.getChr()),1,minimumdivide,1))
        {
            continue;
        }

        if (currentEndRD.TRA1+currentEndRD.TRA2 >= getDivider(readDepthStat.getReadDepthByChr(n.getChr()),1,minimumdivide,1))
        {
           continue;
        }

        if (currentEndRD.INV1+currentEndRD.INV2 >= n.getFrequency())
        {
            continue;
        }

        if (currentEndRD.TRA1+currentEndRD.TRA2 >= n.getFrequency())
        {
           continue;
        }

        if (previousEndRD.INV1+previousEndRD.INV2 >= n.getFrequency())
        {
            continue;
        }

        if (previousEndRD.TRA1+previousEndRD.TRA2 >= n.getFrequency())
        {
           continue;
        }

        if (n.getSvLength() > 2000)
        {
            if (currentEndRD.DUP1 + currentEndRD.DUP2 >= 10)
            {
                continue;
            }

            if (currentEndRD.INV1 + currentEndRD.INV2 >= 10)
            {
                continue;
            }

            if (previousEndRD.DUP1 + previousEndRD.DUP2 >= 10)
            {
                continue;
            }

            if (previousEndRD.INV1 + previousEndRD.INV2 >= 10)
            {
                continue;
            }

            if (nextEndRD.DUP1 + previousEndRD.DUP2 >= 10)
            {
                continue;
            }

            if (nextEndRD.INV1 + previousEndRD.INV2 >= 10)
            {
                continue;
            }
        }

        if (n.getMaxMapQ() < 40)
        {
            continue;
        }

        if (n.getSvLength() > 1000000)
        {
            continue;
        }

        if (n.getSvLength() < 50)
        {
            continue;
        }

        if (n.getFrequency() > 20)
        {
            continue;
        }

        if (n.getSvLength() < 300)
        {
            if (n.getMaxMapQ() < 10)
            {
                continue;
            }
        }
 
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

std::vector<Evidence> RefineDepthBlock::getEvidenceByFilepath(std::string filepaht)
{
    std::vector<Evidence> cache;
    std::string line;
    std::ifstream myfile(filepaht);
    if (myfile.is_open())
    {
        while (getline(myfile, line))
        {
            Evidence e;
            e.setEvidenceByString(line);
            if (e.getPos() != 0)
            {
                // if (e.getSvLength()<500 && e.getAvgMapQ()<6) {

                // } else {
                cache.push_back(e);
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
                std::string tempPath = filemanager->getVariantPath() + "/" + std::string(dir->d_name);
                evidenceFilePathLists.push_back(tempPath);
            }
        }
        closedir(d);
    }

    // d = opendir(filemanager->getSplitReadPath().c_str());
    // if (d)
    // {
    //     while (dir = readdir(d))
    //     {
    //         if (std::string(dir->d_name).size() < 4)
    //         {
    //             continue;
    //         }

    //         if (std::string(dir->d_name).substr(std::string(dir->d_name).size() - 4) == ".txt")
    //         {
    //             std::string tempPath = filemanager->getSplitReadPath() + "/" + std::string(dir->d_name);
    //             evidenceFilePathLists.push_back(tempPath);
    //         }
    //     }
    //     closedir(d);
    // }

    return evidenceFilePathLists;
}

int RefineDepthBlock::getDivider(int value,int top,int down,int minimum) {
    auto returnvalue =  (int) (float(value)*(float(top)/float(down)));

    if (returnvalue>minimum) {
        return returnvalue;
    }

    return minimum;
}