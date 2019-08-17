#include "refinedepthblock.h"
#include <dirent.h>
RefineDepthBlock::RefineDepthBlock()
{
}

int32_t RefineDepthBlock::roundNumber(int32_t number, int32_t round)
{
    if (number % round == 0)
    {
        std::cout << number % round << " is even " << std::endl;
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
        // auto result = getRefineResultDeletion(&variantlist);
        auto result = getResultWithOutOverlapped(&variantlist, &variantlist);

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

        if (n.getSvLength()>2000) {
            if (currentRD.DEL2>=2) {
                continue;
            }
            if (nextRD.DEL2>=2) {
                continue;
            }
            if (previousRD.DEL2>=2) {
                continue;
            }
        }

        if (nextRD.depth>100) {
            continue;
        }

        if (nextRD.depth>previousRD.depth) {
            continue;
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

            if (m.getPos() - 1000 <= n.getPos() && n.getPos() <= m.getPos() + 1000)
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
                cache.push_back(e);
            }

            std::cout
                << e.getChr() << "\t"
                << e.getPos() << "\t"
                << e.getEndChr() << "\t"
                << e.getEnd() << "\t"
                << e.getFrequency() << "\t"
                << e.convertMapQlistToCommaString() << "\t"
                << std::endl;
        }
        myfile.close();
    }

    return cache;
}

void RefineDepthBlock::setFileManager(FileManager *filemanager)
{
    RefineDepthBlock::filemanager = filemanager;
}

std::vector<std::string> RefineDepthBlock::getPathVCFFiles()
{
    std::vector<std::string> evidenceFilePathLists;
    DIR *d;
    struct dirent *dir;

    std::cout << filemanager->getVariantPath() << std::endl;
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

    return evidenceFilePathLists;
}