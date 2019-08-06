#include "readparser.h"
#include <string>
#include <iostream>
#include <algorithm>
#include <cassert>
#include <vector>
#include <cstring>
#include <stdlib.h>
#include <ctype.h>

void ReadParser::setBamRead(bam1_t *bam_read)
{
    source_bamread = bam_read;
}

void ReadParser::setBamHeader(bam_hdr_t *t_bam_header)
{
    bam_header = t_bam_header;
}

ReadParser::ReadParser()
{
}

ReadParser::~ReadParser()
{
}

std::vector<ReadParser::Cigar> ReadParser::getCigar()
{
    const auto cigar = bam_get_cigar(source_bamread);

    std::vector<ReadParser::Cigar> cigarlist;

    for (int k = 0; k < getLengthCigar(); k++)
    {

        const int op = bam_cigar_op(cigar[k]);
        const int ol = bam_cigar_oplen(cigar[k]);

        Cigar cigarL;
        cigarL.setOperatorName(op);
        cigarL.setLength(ol);
        cigarlist.push_back(cigarL);
    }

    return cigarlist;
}

int32_t ReadParser::getPos()
{
    return source_bamread->core.pos + 1;
}

int32_t ReadParser::getEnd()
{
    std::vector<ReadParser::Cigar> cigars = getCigar();
    int count = 0;
    int position = 0;
    int sizeCigar = cigars.size();
    for (Cigar c : cigars)
    {
        if (c.getOperatorName() == 'S' && position == 0)
        {
            count -= c.getLength();
        }
        else if (c.getOperatorName() == 'I')
        {
        }
        else if (c.getOperatorName() == 'H')
        {
        }
        else if (c.getOperatorName() == 'S' && position != sizeCigar)
        {
        }
        else
        {
            count += c.getLength();
        }
        position++;
    }

    return getPos() + count + 1;
}

int32_t ReadParser::getMatePos()
{
    return source_bamread->core.mpos + 1;
}

int32_t ReadParser::getMateEnd()
{
    return source_bamread->core.mpos + source_bamread->core.l_qseq + 1;
}

int32_t ReadParser::getInsertSize()
{
    if (getPos() < getMatePos())
    {
        return getMateEnd() - getPos();
    }

    return getEnd() - getMatePos();
}

int32_t ReadParser::getPosOfSeq()
{
    std::vector<ReadParser::Cigar> cigars = getCigar();

    int32_t shift = 0;
    if (cigars.at(0).getOperatorName() == 'S')
    {
        shift = cigars.at(0).getLength();
    }

    return getPos() - shift;
}

int32_t ReadParser::getEndOfSeq()
{
    return getPosOfSeq() + source_bamread->core.l_qseq;
}

bool ReadParser::isMateUnmapped()
{
    if ((source_bamread->core.flag & BAM_FMUNMAP))
    {
        return true;
    }

    return false;
}

bool ReadParser::isUnmapped()
{
    if ((source_bamread->core.flag & BAM_FUNMAP))
    {
        return true;
    }

    return false;
}

std::string ReadParser::getChromosomeNameString()
{
    std::string chr = bam_header->target_name[source_bamread->core.tid];
    return chr;
}

std::string ReadParser::getMateChromosomeNameString()
{
    std::string mchr = bam_header->target_name[source_bamread->core.mtid];
    return mchr;
}

bool ReadParser::isPairOnSameChromosome()
{
    if (bam_header->target_name[source_bamread->core.tid] == bam_header->target_name[source_bamread->core.mtid])
    {
        return true;
    }

    return false;
}

uint32_t ReadParser::getLengthCigar()
{
    return source_bamread->core.n_cigar;
}

void ReadParser::printCigar()
{
    const auto cigar = bam_get_cigar(source_bamread);

    for (int k = 0; k < getLengthCigar(); k++)
    {
        const int op = bam_cigar_op(cigar[k]);
        const int ol = bam_cigar_oplen(cigar[k]);

        // std::cout << op << "=" << ol << " ";
    }

    // std::cout << std::endl;

    // uint8_t *q = bam_get_seq(source_bamread);

    // int32_t lseq = source_bamread->core.l_qseq;

    // for (int i = 0; i < source_bamread->core.l_qseq; ++i)
    // {
    //     printf("%c", seq_nt16_str[bam_seqi(bam_get_seq(source_bamread), i)]);
    // }
    // printf("\n");

    // for (int i = 0; i < lseq; ++i)
    // {
    //     //   printf("%c", seq_nt16_str[bam_seqi(bam_get_seq(source_bamread), i)]);
    //     printf("%c", seq_nt16_str[bam_seqi(q, i)]);
    // }

    // for (int i = 0; i < lseq; ++i)
    // {
    // printf("qual:");
    // for (int i = 0; i < source_bamread->core.l_qseq; ++i)
    // {
    //     printf("%c", bam_get_qual(source_bamread)[i]);
    // }
    // printf("\n");
    // }

    // std::cout << std::endl
    //           << " = " << q << " / " << source_bamread->core.l_qseq << std::endl;
}

bool ReadParser::hasLastCigarSoftclipped()
{
    const auto cigar = bam_get_cigar(source_bamread);
    int last = getLengthCigar() - 1;
    const int op = bam_cigar_op(cigar[last]);
    const int ol = bam_cigar_oplen(cigar[last]);

    if (op == BAM_CSOFT_CLIP)
    {
        return true;
    }

    return false;
}

bool ReadParser::hasFirstCigarSoftclipped()
{
    const auto cigar = bam_get_cigar(source_bamread);
    const int op = bam_cigar_op(cigar[0]);
    const int ol = bam_cigar_oplen(cigar[0]);

    if (op == BAM_CSOFT_CLIP)
    {
        return true;
    }

    return false;
}

bool ReadParser::hasLastCigarHardclipped()
{
    const auto cigar = bam_get_cigar(source_bamread);
    int last = getLengthCigar() - 1;
    const int op = bam_cigar_op(cigar[last]);
    const int ol = bam_cigar_oplen(cigar[last]);

    if (op == BAM_CHARD_CLIP)
    {
        return true;
    }

    return false;
}

bool ReadParser::hasFirstCigarHardclipped()
{
    const auto cigar = bam_get_cigar(source_bamread);
    const int op = bam_cigar_op(cigar[0]);
    const int ol = bam_cigar_oplen(cigar[0]);

    if (op == BAM_CHARD_CLIP)
    {
        return true;
    }

    return false;
}

std::string ReadParser::getSequence()
{
    std::string buf;
    for (int i = 0; i < source_bamread->core.l_qseq; ++i)
    {
        buf += seq_nt16_str[bam_seqi(bam_get_seq(source_bamread), i)];
    }
    return buf;
}

int *ReadParser::getBaseQuality()
{
    int *s = new int[source_bamread->core.l_qseq];
    int qual;
    for (int i = 0; i < source_bamread->core.l_qseq; ++i)
    {
        qual = bam_get_qual(source_bamread)[i];
        s[i] = qual;
        // std::cout << s[i] << std::endl;
    }

    return s;
}

uint8_t ReadParser::getMapQuality()
{
    return source_bamread->core.qual;
}

char complement(char n)
{
    switch (n)
    {
    case 'A':
        return 'T';
    case 'T':
        return 'A';
    case 'G':
        return 'C';
    case 'C':
        return 'G';
    }

    return ' ';
}

std::string ReadParser::getReverseComplement(std::string nucs)
{
    transform(
        begin(nucs),
        end(nucs),
        begin(nucs),
        complement);

    std::reverse(nucs.begin(), nucs.end());

    return nucs;
}

void ReadParser::replaceToReverseComplement(std::string *nucs)
{
    transform(
        begin(*nucs),
        end(*nucs),
        begin(*nucs),
        complement);

    std::reverse(nucs->begin(), nucs->end());
}

int32_t ReadParser::getLengthSequence()
{
    return source_bamread->core.l_qseq;
}

bool ReadParser::isReverse()
{
    if ((source_bamread->core.flag & BAM_FREVERSE))
    {
        return true;
    }

    return false;
}

bool ReadParser::isMateReverse()
{
    if ((source_bamread->core.flag & BAM_FMREVERSE))
    {
        return true;
    }

    return false;
}

bool ReadParser::isFirstRead()
{
    if ((source_bamread->core.flag & BAM_FREAD1))
    {
        return true;
    }

    return false;
}

bool ReadParser::isSecondRead()
{
    if ((source_bamread->core.flag & BAM_FREAD2))
    {
        return true;
    }

    return false;
}

std::string ReadParser::getEdgeSeqFromStartSeq(uint32_t number)
{
    return getSequence().substr(0, number);
}

std::string ReadParser::getEdgeSeqFromEndSeq(uint32_t number)
{
    return getSequence().substr(source_bamread->core.l_qseq - number, source_bamread->core.l_qseq);
}

bool ReadParser::isNotPassingFilters()
{
    if ((source_bamread->core.flag & BAM_FQCFAIL))
    {
        return true;
    }

    return false;
}

bool ReadParser::isPCR()
{
    if ((source_bamread->core.flag & BAM_FDUP))
    {
        return true;
    }

    return false;
}
bool ReadParser::isSupplementaryAlignment()
{
    if ((source_bamread->core.flag & BAM_FSUPPLEMENTARY))
    {
        return true;
    }

    return false;
}

int ReadParser::getLastToStartMissMatchPosMD()
{
    std::vector<ReadParser::AlignMD> md = getAlignMD();

    int sizeAcc = 0;
    // std::cout << "#md.size() :" << md.size() << std::endl;
    for (int i = md.size() - 1; i >= 0; i--)
    {
        AlignMD lastmd = md.at(i);
        sizeAcc += lastmd.size;
        // std::cout << "#" << lastmd.operate << " "<< lastmd.size << std::endl;
        if (lastmd.operate != 'M')
        {
            break;
        }
    }

    return sizeAcc;
}

int ReadParser::getStartToEndMissMatchPosMD()
{
    std::vector<ReadParser::AlignMD> md = getAlignMD();

    int sizeAcc = 0;
    for (int i = 0; i < md.size(); i++)
    {
        AlignMD lastmd = md.at(i);
        sizeAcc += lastmd.size;
        // std::cout << "#" << lastmd.operate << " = " << lastmd.size << std::endl;
        if (lastmd.operate != 'M')
        {
            break;
        }
    }

    return sizeAcc;
}

std::vector<ReadParser::AlignMD> ReadParser::getAlignMD()
{
    std::vector<ReadParser::AlignMD> alignMDs;
    const char *mdtagchar = "MD";
    auto aux = bam_aux_get(source_bamread, mdtagchar);
    auto auxChar = bam_aux2Z(aux);
    // std::cout << auxChar << std::endl;

    std::string collectNumber;
    std::string operate = "M";
    bool foundDel = false;
    int delCount = 0;
    for (int i = 0; i < strlen(auxChar); i++)
    {
        const char symbol = toupper(auxChar[i]);

        if (isdigit(symbol))
        {
            collectNumber += symbol;

            if (foundDel)
            {
                ReadParser::AlignMD Aalignmd;
                Aalignmd.operate = 'D';
                Aalignmd.size = delCount;
                alignMDs.push_back(Aalignmd);
                collectNumber = "";
                delCount = 0;
                foundDel = false;
            }
        }
        else
        {
            if (foundDel)
            {
                delCount++;
                continue;
            }

            if (collectNumber != "")
            {
                ReadParser::AlignMD Aalignmd;
                Aalignmd.operate = 'M';
                Aalignmd.size = std::stoi(collectNumber);
                alignMDs.push_back(Aalignmd);
                collectNumber = "";
            }

            if (symbol == '^')
            {
                foundDel = true;
                delCount = 0;
            }

            else if (symbol == 'A')
            {
                ReadParser::AlignMD Aalignmd;
                Aalignmd.operate = 'A';
                Aalignmd.size = 1;
                alignMDs.push_back(Aalignmd);
            }
            else if (symbol == 'G')
            {
                ReadParser::AlignMD Aalignmd;
                Aalignmd.operate = 'G';
                Aalignmd.size = 1;
                alignMDs.push_back(Aalignmd);
            }
            else if (symbol == 'T')
            {
                ReadParser::AlignMD Aalignmd;
                Aalignmd.operate = 'T';
                Aalignmd.size = 1;
                alignMDs.push_back(Aalignmd);
            }
            else if (symbol == 'C')
            {
                ReadParser::AlignMD Aalignmd;
                Aalignmd.operate = 'C';
                Aalignmd.size = 1;
                alignMDs.push_back(Aalignmd);
            }
        }
    }

    if (collectNumber != "")
    {
        ReadParser::AlignMD Aalignmd;
        Aalignmd.operate = 'M';
        Aalignmd.size = std::stoi(collectNumber);
        alignMDs.push_back(Aalignmd);
        collectNumber = "";
    }

    // for (auto n : alignMDs)
    // {
    //     std::cout << n.operate << " // " << n.size << std::endl;
    // }

    return alignMDs;
}