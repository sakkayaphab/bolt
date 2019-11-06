#include "faidx.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>

Faidx::Faidx()
{
}

void Faidx::initialize()
{
    std::ifstream file(indexfilepath);
    std::string line;

    while (std::getline(file, line))
    {
        std::istringstream iss(line);
        std::string token;

        FileFormat fftoken;
        int i = 0;
        while (std::getline(iss, token, '\t'))
        {
            if (i == 0)
            {
                fftoken.NAME = token;
            }
            else if (i == 1)
            {
                fftoken.LENGTH = strtoll(token.c_str(), NULL, 10);
            }
            else if (i == 2)
            {
                fftoken.OFFSET = strtoll(token.c_str(), NULL, 10);
            }
            else if (i == 3)
            {
                fftoken.LINEBASES = strtoll(token.c_str(), NULL, 10);
            }
            else if (i == 4)
            {
                fftoken.LINEWIDTH = strtoll(token.c_str(), NULL, 10);
                faidxlist.push_back(fftoken);
            }
            else
            {
                std::cerr << "fasta index error because it has tab more than normal" << std::endl;
                exit(1);
            }
            i++;
        }
    }
}

bool Faidx::hasChromosomeName(std::string chr)
{
    for (int i = 0; i < faidxlist.size(); i++)
    {
        if (faidxlist.at(i).NAME == chr)
        {
            return true;
        }
    }

    return false;
}

uint64_t Faidx::getLengthbyChromosome(std::string chr)
{
    for (int i = 0; i < faidxlist.size(); i++)
    {
        if (faidxlist.at(i).NAME == chr)
        {
            return faidxlist.at(i).LENGTH;
        }
    }
    return 0;
}

uint64_t Faidx::getOffsetbyChromosome(std::string chr)
{
    for (int i = 0; i < faidxlist.size(); i++)
    {
        if (faidxlist.at(i).NAME == chr)
        {
            return faidxlist.at(i).OFFSET;
        }
    }
    return 0;
}

uint64_t Faidx::getLinebasesbyChromosome(std::string chr)
{
    for (int i = 0; i < faidxlist.size(); i++)
    {
        if (faidxlist.at(i).NAME == chr)
        {
            return faidxlist.at(i).LINEBASES;
        }
    }
    return 0;
}

uint64_t Faidx::getLinewidthbyChromosome(std::string chr)
{
    for (int i = 0; i < faidxlist.size(); i++)
    {
        if (faidxlist.at(i).NAME == chr)
        {
            return faidxlist.at(i).LINEWIDTH;
        }
    }
    return 0;
}

uint64_t Faidx::getQualoffsetbyChromosome(std::string chr)
{
    for (int i = 0; i < faidxlist.size(); i++)
    {
        if (faidxlist.at(i).NAME == chr)
        {
            return faidxlist.at(i).QUALOFFSET;
        }
    }
    return 0;
}

void Faidx::setIndexFilePath(std::string t_indexfilepath)
{
    indexfilepath = t_indexfilepath;
}

std::string Faidx::getIndexFilePath()
{
    return indexfilepath;
}

uint64_t Faidx::getApproximateLineStart(std::string chromosome, uint64_t start)
{
    return floor(start / getLinebasesbyChromosome(chromosome));
}

uint64_t Faidx::getApproximateLineEnd(std::string chromosome, uint64_t end)
{
    return ceil(end / getLinebasesbyChromosome(chromosome)) + 1;
}

uint64_t Faidx::getApproximateOffsetStart(std::string chromosome, uint64_t start)
{

    uint64_t apppointstart = getApproximateLineStart(chromosome, start);
    uint64_t diff = start - apppointstart * getLinebasesbyChromosome(chromosome);

    return apppointstart * getLinewidthbyChromosome(chromosome) + diff;
}

uint64_t Faidx::getApproximateOffsetEnd(std::string chromosome, uint64_t end)
{

    uint64_t apppointend = getApproximateLineEnd(chromosome, end);
    uint64_t diff = apppointend * getLinebasesbyChromosome(chromosome) - end;

    return (apppointend * getLinewidthbyChromosome(chromosome) - 1) - diff;
}

uint64_t Faidx::getOffsetEndByPosition(std::string chromosome, uint64_t end)
{
    if (!hasChromosomeName(chromosome))
    {
        std::cout << "not found this chromosome : " << chromosome << std::endl;
    }

    uint64_t startoffset = getOffsetbyChromosome(chromosome);

    return startoffset + getApproximateOffsetEnd(chromosome, end);
}

uint64_t Faidx::getOffsetStartByPosition(std::string chromosome, uint64_t start)
{
    if (!hasChromosomeName(chromosome))
    {
        std::cout << "not found this chromosome : " << chromosome << std::endl;
    }

    uint64_t startoffset = getOffsetbyChromosome(chromosome);

    return startoffset + getApproximateOffsetStart(chromosome, start);
}