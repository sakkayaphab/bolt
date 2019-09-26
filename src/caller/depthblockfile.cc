#include "depthblockfile.h"
#include <map>

DepthBlockFile::DepthBlockFile()
{
}

void DepthBlockFile::execute()
{
}

ReadDepthHelper::ReadDepthVector DepthBlockFile::getBlock(int32_t number)
{
    return mapReadDepthLineSegment[number];
}

ReadDepthHelper::ReadDepthVector DepthBlockFile::findBlockWithFile(int32_t number, std::string filepath, int32_t scope)
{

    std::string line;
    std::ifstream myfile(filepath);
    int count = 0;
    ReadDepthHelper::ReadDepthVector temp;
    if (myfile.is_open())
    {
        while (getline(myfile, line))
        {
            count += scope;
            if (count != number)
            {
                continue;
            }

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
            break;
        }
        myfile.close();
    }

    return temp;
}

std::vector<std::string> DepthBlockFile::split(const std::string &s, char delimiter)
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

std::string DepthBlockFile::getFilePath()
{
    return filepath;
}

void DepthBlockFile::loadDataToCache(std::string filepath)
{
    if (DepthBlockFile::filepath == filepath)
    {
        return;
    }
    DepthBlockFile::filepath = filepath;

    mapReadDepthLineSegment.clear();
    std::string line;
    // int32_t sumRD = 0;
    // int count = 0;

    std::cout << "loadDataToCache : " << filepath << std::endl; 
    std::ifstream myfile(filepath);
    if (myfile.is_open())
    {
        std::vector<std::string> lineBuffer;
        while (getline(myfile, line))
        {
            lineBuffer.push_back(line);

            if (line.size() > 10000)
            {
                addToMapReadDepthLineSegment(&lineBuffer);
                lineBuffer.clear();
            }

            // sumRD += temp.depth;
            // count++;
        }

        if (lineBuffer.size() != 0)
        {
            addToMapReadDepthLineSegment(&lineBuffer);
            lineBuffer.clear();
        }

        myfile.close();
    }
    else
        std::cout << "Unable to open file";

    std::cout << "End + loadDataToCache : " << filepath << std::endl;
    // avgReadDepthFocus = int(sumRD / count);
    // std::cout << "avgReadDepthFocus : " << avgReadDepthFocus << std::endl;
}

void DepthBlockFile::addToMapReadDepthLineSegment(std::vector<std::string> *lineBuffer)
{
    for (auto line : *lineBuffer)
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

        mapReadDepthLineSegment[temp.pos] = temp;
    }
}