#include "fasta.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>

Fasta::Fasta()
{
}

void Fasta::setFilePath(std::string t_filepath)
{
    filepath = t_filepath;
}

std::string Fasta::getFilePath()
{
    return filepath;
}

std::string Fasta::getSeqByOffset(uint64_t start, uint64_t end)
{
    std::ifstream file(filepath);
    file.seekg(start);

    uint64_t sumOffset = start;
    uint64_t currentPosition = file.tellg();

    std::string output;
    std::string line;
    while (std::getline(file, line))
    {
        sumOffset += line.size() + 1;

        if (sumOffset == end)
        {
            output += line.c_str();
            break;
        }
        else if (sumOffset > end)
        {
            uint64_t diff = sumOffset - end;
            output += line.substr(0, line.size() - diff + 1);
            break;
        }

        output += line.c_str();
    }

    return output;
}