#include "splitread.h"

SplitRead::SplitRead(ReadParser *readparser, SampleStat *samplestate)
{
    SplitRead::readparser = readparser;
    SplitRead::samplestate = samplestate;
}

void SplitRead::updateRead()
{
    satag = readparser->getSATag();
    if (satag.size() == 0)
    {
        return;
    }

    findTandemDuplication();

    // if (readparser->isFirstRead()) {
    //     // if (readparser->is)
    // }

    // if (satag.size() != 0)
    // {
    //     for (auto n : satag)
    //     {
    //         std::cout << n.chrname << " " << n.pos << " " << n.strand << std::endl;
    //     }
    // std::cout << "-----------" << std::endl;
    // }
}

void SplitRead::findTandemDuplication()
{

    for (ReadParser::SATag sa : satag)
    {
        if (sa.cigar.size() != 2)
        {
            continue;
        }

        if (readparser->getChromosomeNameString() != sa.chrname)
        {
            continue;
        }
        

        if (readparser->getPos() < sa.pos)
        {
            if (sa.pos - readparser->getPos() > (samplestate->getReadLength()*1.5))
            {
                continue;
            }

            if (readparser->isReverse() && sa.strand == "-")
            {
                goto findTandemDupFirst;
            }

            if (!readparser->isReverse() && sa.strand == "+")
            {
                goto findTandemDupFirst;
            }

            continue;

        findTandemDupFirst:

            if (!(readparser->hasFirstCigarSoftclipped()))
            {
                continue;
            }

            if (sa.cigar.at(0).getOperatorName() == 'M' && sa.cigar.at(1).getOperatorName() == 'S')
            {
                mapDUP[std::make_pair(readparser->getPos(), sa.pos + sa.cigar.at(0).getLength())].NumberOfMatchRead++;
            }
        } else {

            if (readparser->getPos()-sa.pos > (samplestate->getReadLength()*1.5))
            {
                continue;
            }


            if (readparser->isReverse() && sa.strand == "-")
            {
                goto findTandemDupSecond;
            }

            if (!readparser->isReverse() && sa.strand == "+")
            {
                goto findTandemDupSecond;
            }

            continue;

        findTandemDupSecond:

            if (!(readparser->hasLastCigarSoftclipped()))
            {
                continue;
            }

            if (sa.cigar.at(0).getOperatorName() == 'S' && sa.cigar.at(sa.cigar.size()-1).getOperatorName() == 'M')
            {
                mapDUP[std::make_pair( sa.pos,readparser->getEnd())].NumberOfMatchRead++;
            }
        }
    }
}

void SplitRead::printResult() {
    for (auto const &x : mapDUP)
    {
        if (x.second.NumberOfMatchRead>=2) {
            std::cout 
            << " pos : " << x.first.first 
            << " end : " << x.first.second
            << " NumberOfMatchRead : " << x.second.NumberOfMatchRead
            << std::endl;
        }
    }
}