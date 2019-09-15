#include "readdepthanalysis.h"
#include <fstream>
#include <iostream>

ReadDepthAnalysis::ReadDepthAnalysis(FileManager *filemanager)
{
    ReadDepthAnalysis::filemanager = filemanager;
    // readDepthStat = readDepthStat();
    readDepthStat.setFilePath(filemanager);
    readDepthStat.execute();
}

void ReadDepthAnalysis::setSampleStat(SampleStat *samplestat)
{
}

int32_t ReadDepthAnalysis::getRound(int32_t x, int32_t max)
{
    return (x / max) * max;
}

std::vector<std::int32_t> ReadDepthAnalysis::getVectorRange(int32_t pos, int32_t end)
{
    std::vector<std::int32_t> temp;
    for (int32_t i = pos; i <= end; i += configRound)
    {
        temp.push_back(getRound(i, configRound));
    }
    return temp;
}

void ReadDepthAnalysis::setFocusReadDepth(int32_t pos, int32_t end, std::vector<ReadDepthHelper::ReadDepthVector> *focusReadDepth)
{
    ReadDepthHelper::ReadDepthVector temprdvector;

    std::vector<std::int32_t> listrange = getVectorRange(pos, end);

    // std::vector<ReadDepthHelper::ReadDepthVector> vFocus;
    for (auto n : listrange)
    {
        // std::cout << n << std::endl;
        ReadDepthHelper::ReadDepthVector datamodel;
        auto data = mapReadDepthLineSegment[n];
        datamodel.pos = data.pos;
        datamodel.depth = data.depth;
        datamodel.DEL1 = data.DEL1;
        datamodel.DUP1 = data.DUP1;
        datamodel.INS1 = data.INS1;
        datamodel.INV1 = data.INV1;
        datamodel.TRA1 = data.TRA1;

        datamodel.DEL2 = data.DEL2;
        datamodel.DUP2 = data.DUP2;
        datamodel.INS2 = data.INS2;
        datamodel.INV2 = data.INV2;
        datamodel.TRA2 = data.TRA2;

        datamodel.SCF = data.SCF;
        datamodel.SCL = data.SCL;

        focusReadDepth->push_back(datamodel);
    }
}

int ReadDepthAnalysis::getAvgReadDepth()
{
    return avgReadDepth;
}

bool ReadDepthAnalysis::filterDeletion(Evidence e)
{
    int sumStartDELStart = 0;
    int sumStartDELEnd = 0;
    int sumStartDUP = 0;
    int sumStartINV = 0;
    int sumStartTRA = 0;
    int sumStartINS = 0;
    int sumStartSCF = 0;
    int sumStartSCL = 0;
    for (auto n : startFocusReadDepth)
    {
        sumStartDELStart += n.DEL1;
        sumStartDELEnd += n.DEL2;

        sumStartDUP += n.DUP1;
        sumStartINV += n.INV1;
        sumStartTRA += n.TRA1;
        sumStartDUP += n.DUP2;
        sumStartINV += n.INV2;
        sumStartTRA += n.TRA2;
        sumStartINS += n.INS1;
        sumStartINS += n.INS2;
        sumStartSCF += n.SCF;
        sumStartSCL += n.SCL;
    }

    int sumEndDELStart = 0;
    int sumEndDELEnd = 0;
    int sumEndDUP = 0;
    int sumEndINV = 0;
    int sumEndTRA = 0;
    int sumEndINS = 0;
    int sumEndSCF = 0;
    int sumEndSCL = 0;
    for (auto n : startFocusReadDepth)
    {
        sumEndDELStart += n.DEL1;
        sumEndDELEnd += n.DEL2;

        sumEndDUP += n.DUP1;
        sumEndINV += n.INV1;
        sumEndTRA += n.TRA1;
        sumEndDUP += n.DUP2;
        sumEndINV += n.INV2;
        sumEndTRA += n.TRA2;
        sumEndINS += n.INS1;
        sumEndINS += n.INS2;
        sumEndSCF += n.SCF;
        sumEndSCL += n.SCL;
    }

    if (e.getMaxMapQ() < 10)
    {
        return false;
    }

    if (readDepthStat.getReadDepthByChr(e.getChr()) < 15)
    {

    }
    else if (e.getFrequency() < getDivider(readDepthStat.getReadDepthByChr(e.getChr()),1,10,1))
    {
        return false;
    }

    if (sumStartSCL <= 1 && sumStartSCF <= 1)
    {
        return false;
    }

    // if (sumStartTRA >= getDivider(e.getFrequency(), 1, minimumdivide, 1))
    // {
    //     return false;
    // }

    // if (sumEndTRA >= getDivider(e.getFrequency(), 1, minimumdivide, 1))
    // {
    //     return false;
    // }

    // if (sumStartINV >= getDivider(e.getFrequency(), 1, minimumdivide, 1))
    // {
    //     return false;
    // }

    // if (sumEndINV >= getDivider(e.getFrequency(), 1, minimumdivide, 1))
    // {
    //     return false;
    // }

    return true;
}

bool ReadDepthAnalysis::analyzeByEvidence(Evidence e)
{
    if (cachechr != e.getChr())
    {
        loadDataToCache(filemanager->getReadDepthPath() + "/" + e.getChr() + ".txt");
        std::cout << "✓ : " << e.getChr() << std::endl;
        cachechr = e.getChr();
    }

    startFocusReadDepth.clear();
    endFocusReadDepth.clear();

    // return true;
    // if (e.getPos() != 18185538)
    // {
    //     return false;
    // }
    setFocusReadDepth(e.getPos() + e.getCiPosLeft() - configRound, e.getPos() + e.getCiPosRight() + configRound, &startFocusReadDepth);
    setFocusReadDepth(e.getEnd() + e.getCiEndLeft() - configRound, e.getEnd() + e.getCiEndRight() + configRound, &endFocusReadDepth);
    // setFocusReadDepth(10300, 14300);

    // std::cout << "--------" << std::endl;
    // std::cout << "pos : "
    // << e.getPos() + e.getCiPosLeft() - configRound
    // << " end :"
    // << e.getPos() + e.getCiPosRight() + configRound
    // << std::endl;
    // for (auto n : focusReadDepth)
    // {
    //     std::cout << n.pos << " rd:" << n.depth << std::endl;
    // }

    // return false;
    //  return true;
    if (getReadDepthAverageFocusArea(&startFocusReadDepth) > 800)
    {
        return false;
    }

    // if (avgReadDepthFocus>10000) {
    //     return false;
    // }

    if (e.getVariantType() == "DEL")
    {
        if (getReadDepthAverageFocusArea(&startFocusReadDepth) > (readDepthStat.getReadDepthByChr(e.getChr()) * 2))
        {
            return false;
        }

        if (getReadDepthAverageFocusArea(&endFocusReadDepth) > (readDepthStat.getReadDepthByChr(e.getChr()) * 2))
        {
            return false;
        }
       
        return filterDeletion(e);
    }

    if (e.getVariantType() == "INS")
    {
        if (getReadDepthAverageFocusArea(&startFocusReadDepth) > (readDepthStat.getReadDepthByChr(e.getChr()) * 2))
        {
            return false;
        }

        if (getReadDepthAverageFocusArea(&endFocusReadDepth) > (readDepthStat.getReadDepthByChr(e.getChr()) * 2))
        {
            return false;
        }

        if (e.getComment() == "MATEUNMAPPED")
        {
            if (e.getMaxMapQ() == 0)
            {
                return false;
            }

            if (getReadDepthAverageFocusArea(&startFocusReadDepth) > 500)
            {
                return false;
            }

            if (getReadDepthAverageFocusArea(&startFocusReadDepth) < 10)
            {
                return false;
            }

            if (e.getFrequency() < 4)
            {
                return false;
            }

            return true;
        }
        // if (e.getMaxMapQ() == 0)
        // {
        //     return false;
        // }

        // if (e.getAvgMapQ() < 10)
        // {
        //     return false;
        // }

        // int mapqUnzero = e.getFrequency()-e.getNumberOfZeroMapQ();

        // if (e.getMaxMapQ() < 60)
        // {
        //     return false;
        // }

        // if (focusReadDepth.at(0).depth > getAvgReadDepth() * 2)
        // {
        //     return false;
        // }

        // if (focusReadDepth.at(focusReadDepth.size() - 1).depth > getAvgReadDepth() * 2)
        // {
        //     return false;
        // }

        if (getReadDepthAverageFocusArea(&startFocusReadDepth) > 500)
        {
            return false;
        }

        if (getReadDepthAverageFocusArea(&startFocusReadDepth) < 10)
        {
            return false;
        }

        // if (getReadDepthAverageFocusArea() < (getAvgReadDepth() - 10))
        // {
        //     return false;
        // }

        // if (getReadDepthAverageFocusArea() > (getAvgReadDepth() * 2))
        // {
        //     return false;
        // }

        return true;
    }

    if (e.getVariantType() == "DUP")
    {

        if (getReadDepthAverageFocusArea(&startFocusReadDepth) > (readDepthStat.getReadDepthByChr(e.getChr()) * 8))
        {
            return false;
        }

        if (getReadDepthAverageFocusArea(&endFocusReadDepth) > (readDepthStat.getReadDepthByChr(e.getChr()) * 2))
        {
            return false;
        }


        if (e.getMaxMapQ() < 15)
        {
            return false;
        }

        return true;
    }

    if (e.getVariantType() == "INV")
    {
        if (getReadDepthAverageFocusArea(&startFocusReadDepth) > (readDepthStat.getReadDepthByChr(e.getChr()) * 2))
        {
            return false;
        }

        if (getReadDepthAverageFocusArea(&endFocusReadDepth) > (readDepthStat.getReadDepthByChr(e.getChr()) * 2))
        {
            return false;
        }

        if (e.getMaxMapQ() < 15)
        {
            return false;
        }

        if (e.getFrequency() <= 1)
        {
            return false;
        }

        // if (e.getSvLength() > 20000)
        // {
        //     return false;
        // }

        if (getReadDepthAverageFocusArea(&startFocusReadDepth) > 1000)
        {
            return false;
        }

        // if (focusReadDepth.at(0).depth < (getAvgReadDepth() - 15))
        // {
        //     return false;
        // }

        return true;
    }

    return true;
}

int ReadDepthAnalysis::getSCFFocusArea(std::vector<ReadDepthHelper::ReadDepthVector> *focusReadDepth)
{
    int count = 0;
    for (auto n : *focusReadDepth)
    {
        count += n.SCF;
    }

    return count;
}

int ReadDepthAnalysis::getSCLFocusArea(std::vector<ReadDepthHelper::ReadDepthVector> *focusReadDepth)
{
    int count = 0;
    for (auto n : *focusReadDepth)
    {
        count += n.SCL;
    }

    return count;
}

int ReadDepthAnalysis::getReadDepthAverageFocusArea(std::vector<ReadDepthHelper::ReadDepthVector> *focusReadDepth)
{
    int sumRD = 0;
    int count = 0;
    for (auto n : *focusReadDepth)
    {
        sumRD += n.depth;
        count++;
    }
    int avgReadDepthFocus = 0;
    if (count == 0)
    {
        avgReadDepthFocus = 0;
    }
    else
    {
        avgReadDepthFocus = int(sumRD / count);
    }
    return avgReadDepthFocus;
}

void ReadDepthAnalysis::loadDataToCache(std::string filepath)
{
    mapReadDepthLineSegment.clear();

    std::string line;
    int sumRD = 0;
    int count = 0;
    std::ifstream myfile(filepath);
    if (myfile.is_open())
    {
        while (getline(myfile, line))
        {

            ReadDepthHelper::ReadDepthVector temp;
            std::vector<std::string> token = split(line, '\t');
            temp.pos = std::stol(token.at(0), nullptr, 0);
            temp.depth = std::stoi(token.at(1));
            temp.DEL1 = std::stoi(token.at(2));
            temp.DUP1 = std::stoi(token.at(3));
            temp.INS1 = std::stoi(token.at(4));
            temp.INV1 = std::stoi(token.at(5));
            temp.TRA1 = std::stoi(token.at(6));

            temp.DEL2 = std::stoi(token.at(8));
            temp.DUP2 = std::stoi(token.at(9));
            temp.INS2 = std::stoi(token.at(10));
            temp.INV2 = std::stoi(token.at(11));
            temp.TRA2 = std::stoi(token.at(12));

            temp.SCF = std::stoi(token.at(13));
            temp.SCL = std::stoi(token.at(14));

            sumRD += temp.depth;
            count++;
            mapReadDepthLineSegment[temp.pos] = temp;

            // cacheReadDepthFile.push_back(temp);
        }
        myfile.close();
    }
    else
        std::cout << "Unable to open file";

    avgReadDepthFocus = int(sumRD / count);
    // std::cout << "avgReadDepthFocus : " << avgReadDepthFocus << std::endl;
}

// void ReadDepthAnalysis::loadReadDepthStat() {
//     std::string filepath = filemanager->getReadDepthStatPath()+"/readdepthstat.txt";
//     std::string line;
//     std::ifstream myfile(filepath);
//     if (myfile.is_open())
//     {
//         while (getline(myfile, line))
//         {
//             ReadDepthHelper::ReadDepthVector temp;
//             std::vector<std::string> token = split(line, '=');

//             readdepthlist[std::stoi(token.at(1))]++;
//         }
//         myfile.close();
//     }
//     else
//         std::cout << "Unable to open file";

// }

std::vector<std::string> ReadDepthAnalysis::split(const std::string &s, char delimiter)
{
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter))
    {
        tokens.push_back(token);
    }
    return tokens;
}

int ReadDepthAnalysis::getDivider(int value, int top, int down, int minimum)
{
    auto returnvalue = (int)(float(value) * (float(top) / float(down)));

    if (returnvalue > minimum)
    {
        return returnvalue;
    }

    return minimum;
}