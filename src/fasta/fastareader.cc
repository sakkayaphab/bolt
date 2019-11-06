#include "fastareader.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

FastaReader::FastaReader()
{
}

void FastaReader::initialize()
{
    faidx.initialize();
}

void FastaReader::setIndexFilePath(std::string t_indexfilepath)
{
    faidx.setIndexFilePath(t_indexfilepath);
}

std::string FastaReader::getIndexFilePath()
{
    return faidx.getIndexFilePath();
}

void FastaReader::setFilePath(std::string filepath)
{
    fasta.setFilePath(filepath);
}

std::string FastaReader::getFilePath()
{
    return fasta.getFilePath();
}

std::string FastaReader::getSeqbyPosition(std::string chromosome,uint64_t start,uint64_t  end) {
    uint64_t offsetstart = faidx.getOffsetStartByPosition(chromosome,start);
    uint64_t offsetend  = faidx.getOffsetEndByPosition(chromosome,end);

    return fasta.getSeqByOffset(offsetstart,offsetend);
}