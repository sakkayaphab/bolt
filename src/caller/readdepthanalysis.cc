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

        // std::cout << data.SCF << " " << data.SCL << " > " << data.pos << " " << data.depth << std::endl;

        focusReadDepth->push_back(datamodel);
    }
}

int ReadDepthAnalysis::getAvgReadDepth()
{
    return avgReadDepth;
}

bool ReadDepthAnalysis::filterInversion(Evidence e)
{

    if (sumStartDEL >= getDivider(readDepthStat.getReadDepthByChr(e.getChr()), 3, 10, 1))
    {
        return false;
    }

    if (sumEndDEL >= getDivider(readDepthStat.getReadDepthByChr(e.getChr()), 3, 10, 1))
    {
        return false;
    }

    if (readDepthStat.getReadDepthByChr(e.getChr()) < 15)
    {
    }
    else if (e.getFrequency() < getDivider(readDepthStat.getReadDepthByChr(e.getChr()), 1, 20, 1))
    {
        return false;
    }

    // if (sumStartINV - e.getFrequency()>=getDivider(readDepthStat.getReadDepthByChr(e.getChr()),3,10,1)) {
    //     return false;
    // }

    // if (sumEndINV - e.getFrequency()>=getDivider(readDepthStat.getReadDepthByChr(e.getChr()),3,10,1)) {
    //     return false;
    // }

    if (e.getFrequency() <= sumStartINV - e.getFrequency())
    {
        return false;
    }

    if (e.getFrequency() <= sumEndINV - e.getFrequency())
    {
        return false;
    }

    if (getReadDepthAverageFocusArea(&startFocusReadDepth) > (readDepthStat.getReadDepthByChr(e.getChr()) * 2))
    {
        return false;
    }

    if (getReadDepthAverageFocusArea(&endFocusReadDepth) > (readDepthStat.getReadDepthByChr(e.getChr()) * 2))
    {
        return false;
    }

    if (e.getMaxMapQ() < 40)
    {
        return false;
    }

    if (e.getFrequency() <= 1)
    {
        return false;
    }

    if (getReadDepthAverageFocusArea(&startFocusReadDepth) > 1000)
    {
        return false;
    }

    if (sumStartTRA >= e.getFrequency())
    {
        return false;
    }

    if (sumEndTRA >= e.getFrequency())
    {
        return false;
    }

    if (sumStartDEL >= e.getFrequency())
    {
        return false;
    }

    if (sumEndDEL >= e.getFrequency())
    {
        return false;
    }

    return true;
}

bool ReadDepthAnalysis::filterDeletion(Evidence e)
{

    if (e.getFrequency() <= 1)
    {
        return false;
    }

    if (readDepthStat.getReadDepthByChr(e.getChr()) < 15)
    {
    }
    else if (e.getFrequency() < getDivider(readDepthStat.getReadDepthByChr(e.getChr()), 2, 100, 1))
    {
        return false;
    }

    if (sumStartSCL <= 1 && sumEndSCF <= 1)
    {
        return false;
    }

    return true;
}

bool ReadDepthAnalysis::filterInsertion(Evidence e)
{

    // std::cout << "getReadDepthAverageFocusArea " << getReadDepthAverageFocusArea(&startFocusReadDepth) << " " << (readDepthStat.getReadDepthByChr(e.getChr()) * 2) << std::endl;

    // if (getReadDepthAverageFocusArea(&startFocusReadDepth) > (readDepthStat.getReadDepthByChr(e.getChr()) * 2))
    // {
    //     return false;
    // }

    if ((e.getMark() == "MATEUNMAPPED"))
    {
        if (e.getFrequency() <= 4)
        {
            return false;
        }
    }
    else if ((e.getMark() == "SINS"))
    {
    }
    else if ((e.getMark() == "SR"))
    {
    }
    else
    {
        if (e.getFrequency() <= 3)
        {
            return false;
        }
        
    }

    if (sumStartSCL <= 3 && sumStartSCF <= 3)
    {
        return false;
    }

    return true;
}

bool ReadDepthAnalysis::filterDuplication(Evidence e)
{
    if (e.getMark() == "SR")
    {
        return true;
    }

    if (sumStartSCL <= 1 && sumStartSCF <= 1)
    {
        return false;
    }

    if (e.getFrequency() <= 1)
    {
        return false;
    }

    if (e.getMaxMapQ() < 20)
    {
        return false;
    }

    // if (readDepthStat.getReadDepthByChr(e.getChr()) < 15)
    // {

    // }
    if (e.getFrequency() <= getDivider(readDepthStat.getReadDepthByChr(e.getChr()), 5, 100, 1))
    {
        return false;
    }

    if (getReadDepthAverageFocusArea(&startFocusReadDepth) > (readDepthStat.getReadDepthByChr(e.getChr()) * 8))
    {
        return false;
    }

    if (getReadDepthAverageFocusArea(&endFocusReadDepth) > (readDepthStat.getReadDepthByChr(e.getChr()) * 8))
    {
        return false;
    }

    return true;
}

bool ReadDepthAnalysis::analyzeByEvidence(Evidence e)
{
    if (e.getMark() == "SR")
    {
        return true;
    }

    if (e.getMark() == "SDEL")
    {
        return true;
    }

    if (cachechr != e.getChr())
    {
        std::cout << e.convertToVcfString() << std::endl;

        loadDataToCache(filemanager->getReadDepthPath() + "/" + e.getChr() + ".txt");
        std::cout << "✓ : " << e.getChr() << std::endl;
        cachechr = e.getChr();
    }

    startFocusReadDepth.clear();
    endFocusReadDepth.clear();

    if (e.getVariantType() == "DEL")
    {
        setFocusReadDepth(e.getPos() + e.getCiPosLeft() - configRound, e.getPos() + e.getCiPosRight() + configRound, &startFocusReadDepth);
        setFocusReadDepth(e.getEnd() + e.getCiEndLeft() - configRound, e.getEnd() + e.getCiEndRight() + configRound, &endFocusReadDepth);
        collectNewData();

        if (getReadDepthAverageFocusArea(&startFocusReadDepth) > 500 && getReadDepthAverageFocusArea(&endFocusReadDepth) > 500)
        {
            return false;
        }

        return filterDeletion(e);
    }

    if (e.getVariantType() == "INS")
    {
        setFocusReadDepth(e.getPos() + e.getCiPosLeft() - configRound, e.getPos() + e.getCiPosRight() + configRound, &startFocusReadDepth);
        collectNewData();

        if (getReadDepthAverageFocusArea(&startFocusReadDepth) > 400)
        {
            return false;
        }

        return filterInsertion(e);
    }

    if (e.getVariantType() == "DUP")
    {
        setFocusReadDepth(e.getPos() + e.getCiPosLeft() - configRound, e.getPos() + e.getCiPosRight() + configRound, &startFocusReadDepth);
        setFocusReadDepth(e.getEnd() + e.getCiEndLeft() - configRound, e.getEnd() + e.getCiEndRight() + configRound, &endFocusReadDepth);
        collectNewData();

        if (getReadDepthAverageFocusArea(&startFocusReadDepth) > 800 && getReadDepthAverageFocusArea(&endFocusReadDepth) > 800)
        {
            return false;
        }

        return filterDuplication(e);
    }

    if (e.getVariantType() == "INV")
    {
        setFocusReadDepth(e.getPos() + e.getCiPosLeft() - configRound, e.getPos() + e.getCiPosRight() + configRound, &startFocusReadDepth);
        setFocusReadDepth(e.getEnd() + e.getCiEndLeft() - configRound, e.getEnd() + e.getCiEndRight() + configRound, &endFocusReadDepth);
        collectNewData();

        return filterInversion(e);
    }

    if (e.getVariantType() == "BND")
    {
        setFocusReadDepth(e.getPos() + e.getCiPosLeft() - configRound, e.getPos() + e.getCiPosRight() + configRound, &startFocusReadDepth);
        setFocusReadDepth(e.getEnd() + e.getCiEndLeft() - configRound, e.getEnd() + e.getCiEndRight() + configRound, &endFocusReadDepth);
        collectNewData();

        return filterTranslocation(e);
    }

    return false;
}

bool ReadDepthAnalysis::filterTranslocation(Evidence e)
{
    if (e.getFrequency() < getDivider(readDepthStat.getReadDepthByChr(e.getChr()), 4, 10, 1))
    {
        return false;
    }

    if (sumStartDEL >= getDivider(readDepthStat.getReadDepthByChr(e.getChr()), 3, 10, 1))
    {
        return false;
    }

    // if (sumEndDEL>=getDivider(readDepthStat.getReadDepthByChr(e.getChr()),3,10,1)) {
    //     return false;
    // }

    if (e.getFrequency() <= sumStartINV - e.getFrequency())
    {
        return false;
    }

    // if (e.getFrequency() <= sumEndINV - e.getFrequency())
    // {
    //     return false;
    // }

    if (getReadDepthAverageFocusArea(&startFocusReadDepth) > (readDepthStat.getReadDepthByChr(e.getChr()) * 2))
    {
        return false;
    }

    // if (getReadDepthAverageFocusArea(&endFocusReadDepth) > (readDepthStat.getReadDepthByChr(e.getChr()) * 2))
    // {
    //     return false;
    // }
    if (e.getMinMapQ() == 0)
    {
        if (e.getFrequency() < getDivider(readDepthStat.getReadDepthByChr(e.getChr()), 6, 10, 1))
        {
            return false;
        }
        return false;
    }

    if (e.getMaxMapQ() < 60)
    {
        return false;
    }

    if (e.getFrequency() <= 1)
    {
        return false;
    }

    if (getReadDepthAverageFocusArea(&startFocusReadDepth) > 1000)
    {
        return false;
    }

    if (sumStartINV >= e.getFrequency())
    {
        return false;
    }

    // if (sumEndINV >= e.getFrequency())
    // {
    //     return false;
    // }

    if (sumStartDEL >= e.getFrequency())
    {
        return false;
    }

    // if (sumEndDEL >= e.getFrequency())
    // {
    //     return false;
    // }

    return true;
}

void ReadDepthAnalysis::collectNewData()
{
    sumStartDEL = 0;
    sumStartDUP = 0;
    sumStartINV = 0;
    sumStartTRA = 0;
    sumStartINS = 0;
    sumStartSCF = 0;
    sumStartSCL = 0;

    sumEndDEL = 0;
    sumEndDUP = 0;
    sumEndINV = 0;
    sumEndTRA = 0;
    sumEndINS = 0;
    sumEndSCF = 0;
    sumEndSCL = 0;

    sumStartR1_MUN = 0;
    sumStartR2_MUN = 0;

    sumEndR1_MUN = 0;
    sumEndR2_MUN = 0;

    for (auto n : startFocusReadDepth)
    {
        sumStartDEL += n.DEL1;
        sumStartDEL += n.DEL2;

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

        sumStartR1_MUN += n.R1_MUN;
        sumStartR2_MUN += n.R1_MUN;
    }

    for (auto n : endFocusReadDepth)
    {
        sumEndDEL += n.DEL1;
        sumEndDEL += n.DEL2;

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

        sumEndR1_MUN += n.R1_MUN;
        sumEndR2_MUN += n.R1_MUN;
    }
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

            temp.R1_MUN = std::stoi(token.at(15));
            temp.R2_MUN = std::stoi(token.at(16));

            sumRD += temp.depth;
            count++;
            mapReadDepthLineSegment[temp.pos] = temp;

            // cacheReadDepthFile.push_back(temp);
        }
        myfile.close();
    }
    else
        std::cout << "Unable to open file";

    if (sumRD == 0)
    {
        avgReadDepthFocus = 0;
    }
    else
    {
        avgReadDepthFocus = int(sumRD / count);
    }

    // std::cout << "avgReadDepthFocus : " << avgReadDepthFocus << std::endl;
}

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