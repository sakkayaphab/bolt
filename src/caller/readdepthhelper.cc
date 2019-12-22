#include "readdepthhelper.h"
#include "evidencereaddepthfilter.h"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>

ReadDepthHelper::ReadDepthHelper()
{
}

ReadDepthHelper::ReadDepthHelper(FileManager *filepath)
{
    ReadDepthHelper::filepath = filepath;
}

void ReadDepthHelper::setReadDepthMap(std::map<int32_t, EvidenceFinder::ReadDepthDetail> ReadDepthLineSegment)
{
    for (std::map<int32_t, EvidenceFinder::ReadDepthDetail>::iterator it = ReadDepthLineSegment.begin(); it != ReadDepthLineSegment.end(); ++it)
    {
        ReadDepthVector rdv;
        rdv.pos = it->first;
        rdv.depth = it->second.RD;
        rdv.DEL1 = it->second.DEL1;
        rdv.INS1 = it->second.INS1;
        rdv.DUP1 = it->second.DUP1;
        rdv.INV1 = it->second.INV1;
        rdv.TRA1 = it->second.TRA1;

        rdv.DEL2 = it->second.DEL2;
        rdv.INS2 = it->second.INS2;
        rdv.DUP2 = it->second.DUP2;
        rdv.INV2 = it->second.INV2;
        rdv.TRA2 = it->second.TRA2;

        rdv.SCF = it->second.SCF;
        rdv.SCL = it->second.SCL;

        rdv.R1_MUN = it->second.R1_MUN;
        rdv.R2_MUN = it->second.R2_MUN;

        vecRDLine.push_back(rdv);
    }

    std::sort(vecRDLine.begin(), vecRDLine.end());
}

void ReadDepthHelper::writeReadDepthLineFile(std::string path)
{
    std::ofstream myfile;
    myfile.open(path, std::ios_base::app);
    int32_t currentPos = 0;
    bool first = true;
    for (auto n : vecRDLine)
    {
        

        if (currentPos+range < n.pos)
        {
            if (n.pos==0) {
                continue;
            }

            for (; currentPos < n.pos;)
            {
                currentPos += range;
                myfile << currentPos << "\t"
                       << 0 << "\t"
                       << 0 << "\t"
                       << 0 << "\t"
                       << 0 << "\t"
                       << 0 << "\t"
                       << 0 << "\t"
                       << "="
                       << "\t"
                       << 0 << "\t"
                       << 0 << "\t"
                       << 0 << "\t"
                       << 0 << "\t"
                       << 0 << "\t"
                       << 0 << "\t"
                       << 0 << "\t"
                       << 0 << "\t"
                       << 0 << std::endl;
            }
        }
        else
        {
            myfile << n.pos << "\t"
                   << n.depth << "\t"
                   << n.DEL1 << "\t"
                   << n.DUP1 << "\t"
                   << n.INS1 << "\t"
                   << n.INV1 << "\t"
                   << n.TRA1 << "\t"
                   << "="
                   << "\t"
                   << n.DEL2 << "\t"
                   << n.DUP2 << "\t"
                   << n.INS2 << "\t"
                   << n.INV2 << "\t"
                   << n.TRA2 << "\t"
                   << n.SCF << "\t"
                   << n.SCL << "\t"
                   << n.R1_MUN << "\t"
                   << n.R2_MUN << std::endl;
            currentPos = n.pos;
        }
    }
    myfile.close();
}

std::vector<std::string> ReadDepthHelper::split(const std::string &s, char delimiter)
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

void ReadDepthHelper::loadReadDepthFile(std::string path)
{
    std::string line;
    std::ifstream myfile(path);
    if (myfile.is_open())
    {
        while (getline(myfile, line))
        {
            auto list = split(line, '\t');

            ReadDepthVector rdv;
            rdv.pos = std::stol(list.at(0));
            rdv.depth = std::stol(list.at(1));
            vecRDLine.push_back(rdv);
            // std::cout << rdv.posDiscordantRead << " / " << rdv.depth << '\n';
        }
        myfile.close();
    }
    else
    {
        std::cout << "Unable to open file";
    }
}

void ReadDepthHelper::calculateAvgReaddepth()
{
    uint64_t sumRD = 0;
    for (auto n : vecRDLine)
    {
        if (n.depth < 1000 && n.depth >= 0)
        {
            sumRD += n.depth;
        }
    }

    uint64_t size = vecRDLine.size();

    if (size != 0)
    {
        int avg = int(sumRD / size);

        setAvgReadDepth(avg);
        std::cout << target_chromosome << " size : " << size << " rd : " << sumRD << " AVG : " << avg << std::endl;
    }
    else
    {
        setAvgReadDepth(0);
    }

}

int ReadDepthHelper::getAvgReadDepth()
{
    return avgReadDepth;
}

void ReadDepthHelper::setAvgReadDepth(int avgReadDepth)
{
    ReadDepthHelper::avgReadDepth = avgReadDepth;
}

void ReadDepthHelper::findEvidence()
{
    EvidenceReadDepthFilter evidenceRDF(filepath);
    int32_t pCurrentPos;
    int32_t pCurrentEnd;
    bool start = false;
    int count = 0;
    for (auto n : vecRDLine)
    {
        count++;
        if (n.depth < 5)
        {
            if (!start)
            {
                if (count - 2 >= 0)
                {
                    pCurrentPos = vecRDLine.at(count - 2).pos;
                }
                else
                {
                    pCurrentPos = vecRDLine.at(count - 1).pos;
                }

                start = true;
            }

            continue;
        }
        else
        {
            if (start)
            {
                pCurrentEnd = n.pos;
                int32_t diff = pCurrentEnd - pCurrentPos;
                if (diff < 50000)
                {
                    if (evidenceRDF.haveNoneSymbol(target_chromosome, pCurrentPos, pCurrentEnd))
                    {
                    }
                    else
                    {
                        VariantRangeRD vrrd;
                        vrrd.pos = pCurrentPos;
                        vrrd.end = pCurrentEnd;
                        LowRD.push_back(vrrd);
                    }
                }

                start = false;
            }

            continue;
        }
    }
}

void ReadDepthHelper::setRange(int32_t mrange)
{
    range = mrange;
}

void ReadDepthHelper::setTargetChromosome(std::string chr)
{
    target_chromosome = chr;
}

void ReadDepthHelper::setOutputPath(std::string path)
{
    outputpath = path;
}

std::vector<std::string> ReadDepthHelper::getVariantVcfFormat()
{
    std::vector<std::string> variantvcf;

    for (auto n : LowRD)
    {
        std::string buf;
        //CHROM
        buf.append(target_chromosome + "\t");
        // POS
        buf.append(std::to_string(n.pos) + "\t");
        // ID
        buf.append(".\t");
        // REF
        buf.append(".\t");
        // ALT
        buf.append(".\t");
        // QUAL
        buf.append(".\t");
        // FILTER
        buf.append(".\t");
        // INFO
        buf.append("END=" + std::to_string(n.end) + ";");
        buf.append("SVTYPE=DEL;");
        buf.append("SVLEN=" + std::to_string(n.end - n.pos) + ";");
        variantvcf.push_back(buf);
    }

    return variantvcf;
}

void ReadDepthHelper::writeVcf(std::string path)
{
    auto vvfLists = getVariantVcfFormat();
    std::ofstream myfile;
    myfile.open(path, std::ios_base::app);
    for (auto n : vvfLists)
    {
        myfile << n << std::endl;
    }
    myfile.close();
}

bool ReadDepthHelper::isRangeDisorderByMorethanRD(int32_t pos, int32_t end, int readdepth)
{
    std::vector<ReadDepthVector> collect;
    bool firstSector;

    for (int i = 0; i < vecRDLine.size(); i++)
    {
        if (!firstSector)
        {
            if (vecRDLine.at(i).pos >= pos)
            {
                collect.push_back(vecRDLine.at(i));
                firstSector = true;
            }
        }
        else
        {
            if (vecRDLine.at(i).pos > end)
            {
                break;
            }
            collect.push_back(vecRDLine.at(i));
        }
    }

    for (auto n : collect)
    {
        if (n.depth > readdepth)
        {
            return true;
        }
    }

    return false;
}

void ReadDepthHelper::writeReadDepthStat(std::string path)
{
    filepath->mutexLockReadDepthStat();
    std::ofstream myfile;
    myfile.open(path, std::ios_base::app);
    myfile << target_chromosome << "=" << std::to_string(getAvgReadDepth()) << std::endl;
    myfile.close();
    filepath->mutexUnlockReadDepthStat();
}
