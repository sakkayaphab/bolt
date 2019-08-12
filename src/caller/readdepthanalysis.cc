#include "readdepthanalysis.h"
#include <fstream>
#include <iostream>

ReadDepthAnalysis::ReadDepthAnalysis(FileManager *filemanager)
{
    ReadDepthAnalysis::filemanager = filemanager;
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
    int sumDELStart = 0;
    int sumDELEnd = 0;
    int sumDUP = 0;
    int sumINV = 0;
    int sumTRA = 0;
    int sumINS = 0;
    for (auto n : startFocusReadDepth)
    {
        sumDELStart += n.DEL1;
        sumDELEnd += n.DEL2;

        sumDUP += n.DUP1;
        sumINV += n.INV1;
        sumTRA += n.TRA1;
        sumDUP += n.DUP2;
        sumINV += n.INV2;
        sumTRA += n.TRA2;
        sumINS += n.INS1;
        sumINS += n.INS2;

        // std::cout << " ----- " << "\n"
        // << "p : " << n.pos << " = "
        // << sumTRA << " // " << sumDELStart << std::endl;
    }

    int32_t svlength = e.getEndDiscordantRead() - e.getPosDiscordantRead();

    // if (sumTRA >= e.getFrequency())
    // {
    //     return false;
    // }

    // if (sumINV >= e.getFrequency())
    // {
    //     return false;
    // }

    // if (n.DEL1!=0) {
    // std::cout << e.getPos() << "\t" << n.DEL1 << "\t" << n.DEL2 << std::endl;
    // }

    if (getReadDepthAverageFocusArea(&startFocusReadDepth) > getAvgReadDepth() * 1.5)
    {
        return false;
    }

    if (e.getMaxMapQ() == 0)
    {
        return false;
    }

    if (getAvgReadDepth() > 100)
    {
        if (e.getFrequency() <= 2)
        {

            return false;
        }
    }

    // if (e.getFrequency()<=1) {
    //     return false;
    // }

    // if (e.getSvLength() < 500)
    // {
    //     if (sumINV > 50)
    //     {
    //         return false;
    //     }

    //     if (sumINS > 50)
    //     {
    //         return false;
    //     }

    //     if (sumTRA > 50)
    //     {
    //         return false;
    //     }
    // }

    if (e.getSvLength() > 2000)
    {
        // if (e.getFrequency()<=1) {
        //     if (e.getMaxMapQ()==0) {
        //         return false;
        //     }
        // }

        // if ( e.getMaxMapQ()<10) {
        //     return false;
        // }

        // if (e.getFrequency() <= 5)
        // {
        // return false;
        // }
        // return false;
    }
    else if (e.getSvLength() > 1000 && e.getSvLength() <= 2000)
    {
        //    if (e.getFrequency()<=1) {
        //         if (e.getMaxMapQ()==0) {
        //             return false;
        //         }
        //     }
        // if ( e.getMaxMapQ()<35) {
        //     return false;
        // }

        // if ( e.getAvgMapQ()<20) {
        //     return false;
        // }

        // if (e.getFrequency() <= 3)
        // {
        //     return false;
        // }

        // if (sumDEL<5) {
        //     return false;
        // }

        // if (getSCLFocusArea(&startFocusReadDepth) >= 2 || getSCFFocusArea(&endFocusReadDepth) >= 2)
        // {
        // }
        // else
        // {
        //     return false;
        // }

        // if (sumINV+sumDUP+sumTRA>10) {
        //     return false;
        // }
        // return false;
        // return false;
    }
    else if (e.getSvLength() > 500 && e.getSvLength() <= 1000)
    {
        // if (e.getFrequency()<=1) {
        //     if (e.getMaxMapQ()==0) {
        //         return false;
        //     }
        // }

        // if (e.getFrequency() <= 2 && e.getAvgMapQ() == 60)
        // {
        //     return false;
        // }

        // if (getSCLFocusArea(&startFocusReadDepth) >= 2 || getSCFFocusArea(&endFocusReadDepth) >= 2)
        // {
        // }
        // else
        // {
        //     return false;
        // }

        // if (sumDUP>5) {
        //     return false;
        // }

        // return false;
        // return false;
    }
    else
    {
        if (e.getFrequency() <= 3)
        {
            return false;
        }

        // if (getSCLFocusArea(&startFocusReadDepth) >= 10 || getSCFFocusArea(&endFocusReadDepth) >= 10)
        // {
        //     std::cout << e.getSvLength() << std::endl;
        // }
        // else
        // {
        //     return false;
        // }

        // if (e.getSvLength())
        // if (e.getMaxMapQ() < 40)
        // {
        //     return false;
        // }

        // if (e.getFrequency() <= 3)
        // {
        //     return false;
        // }

        // if (sumINV+sumDUP+sumTRA>5) {
        //     return false;
        // }

        // return false;
    }

    // if (getReadDepthAverageFocusArea() > 500)
    // {
    //     return false;
    // }

    // if (e.getMaxMapQ() < 30)
    // {
    //     return false;
    // }

    // if (getReadDepthAverageFocusArea()>avgReadDepth*2) {
    //     return false;
    // }

    // if (e.getAvgMapQ()>20) {
    //     return false;
    // }

    // if (e.getMaxMapQ() < 40)
    // {
    //     return false;
    // }

    // if (getReadDepthAverageFocusArea() > (getAvgReadDepth() * 4))
    // {
    //     return false;
    // }

    // if (e.getEndDiscordantRead() - e.getLastPosDiscordantRead() < 10)
    // {
    //     return false;
    // }

    // if (focusReadDepth.size() > 3)
    // {
    //     if ((focusReadDepth.at(0).depth) > focusReadDepth.at(1).depth)
    //     {
    //         return false;
    //     }

    //     if ((focusReadDepth.at(focusReadDepth.size() - 2).depth) > focusReadDepth.at(focusReadDepth.size() - 1).depth)
    //     {
    //         return false;
    //     }
    // }
    // else
    // {
    //     if (focusReadDepth.at(0).depth > getAvgReadDepth())
    //     {
    //         return false;
    //     }
    // }

    return true;
}

bool ReadDepthAnalysis::analyzeByEvidence(Evidence e)
{
    if (cachechr != e.getChr())
    {
        loadDataToCache(filemanager->getReadDepthPath() + "/" + e.getChr() + ".txt");
        std::cout << "e.getChr() : " << e.getChr() << std::endl;
        cachechr = e.getChr();
    }

    startFocusReadDepth.clear();
    endFocusReadDepth.clear();
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

    if (getReadDepthAverageFocusArea(&startFocusReadDepth) > 800)
    {
        return false;
    }

    // if (avgReadDepthFocus>10000) {
    //     return false;
    // }

    if (e.getVariantType() == "DEL")
    {
        return filterDeletion(e);
    }

    if (e.getVariantType() == "INS")
    {
        // return filterIns

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

        // if (e.getSvLength() < 1500)
        // {
        //     if (e.getMaxMapQ() < 30)
        //     {
        //         return false;
        //     }

        //     if (e.getFrequency() < 2)
        //     {
        //         return false;
        //     }

        //     return true;
        // }

        //    if (focusReadDepth.size() > 3)
        //     {
        //         if ((focusReadDepth.at(0).depth) < avgReadDepth)
        //         {
        //             return false;
        //         }

        //         if ((focusReadDepth.at(focusReadDepth.size() - 1).depth) < avgReadDepth)
        //         {
        //             return false;
        //         }
        //     }

        // if (e.getMaxMapQ() < 20)
        // {
        //     return false;
        // }

        // if (focusReadDepth.at(0).depth < (getAvgReadDepth() - 10))
        // {
        //     return false;
        // }

        // if (focusReadDepth.at(focusReadDepth.size() - 1).depth < (getAvgReadDepth() - 10))
        // {
        //     return false;
        // }

        // if (getReadDepthAverageFocusArea() > (getAvgReadDepth() * 3))
        // {
        //     returnfalse;
        // }

        // if (getReadDepthAverageFocusArea() < avgReadDepth-20)
        // {
        //     return false;
        // }

        // if (getReadDepthAverageFocusArea(&startFocusReadDepth) > 500)
        // {
        //     return false;
        // }

        return true;
    }

    if (e.getVariantType() == "INV")
    {
        // if (e.getMaxMapQ() < 10)
        // {
        //     return false;
        // }

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
    std::cout << "avgReadDepthFocus : " << avgReadDepthFocus << std::endl;
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

void ReadDepthAnalysis::loadAvgReadDepthStat()
{
    std::string filepath = filemanager->getReadDepthStatPath() + "/readdepthstat.txt";
    std::string line;
    std::ifstream myfile(filepath);
    int sum = 0;
    int count = 0;
    if (myfile.is_open())
    {
        while (getline(myfile, line))
        {
            ReadDepthHelper::ReadDepthVector temp;
            std::vector<std::string> token = split(line, '=');

            if (std::stoi(token.at(1)) < 10)
            {
                continue;
            }

            if (std::stoi(token.at(1)) > 800)
            {
                continue;
            }

            // std::cout <<  token.at(1) << std::endl;
            sum += std::stoi(token.at(1));
            count++;
        }
        myfile.close();
    }
    else
        std::cout << "Unable to open file";

    avgReadDepth = sum / count;
    // std::cout << sum << " " << count << std::endl;
    // std::cout << avgReadDepth << std::endl;
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
