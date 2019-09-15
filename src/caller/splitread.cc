#include "splitread.h"
#include <iostream>
#include <fstream>
SplitRead::SplitRead(std::string chrname, ReadParser *readparser, SampleStat *samplestate, FileManager *filepath)
{
    SplitRead::readparser = readparser;
    SplitRead::samplestate = samplestate;
    SplitRead::chrname = chrname;
    SplitRead::filepath = filepath;
}

void SplitRead::updateRead()
{
    satag = readparser->getSATag();
    if (satag.size() == 0)
    {
        return;
    }

    // findInversion();
    findDeletion();
    // findTandemDuplication();

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

void SplitRead::findInversion()
{
    // std::cout << "-----------" << std::endl;
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

            if (readparser->isReverse() && sa.strand == "+")
            {
                goto findInvFirst;
            }

            if (!readparser->isReverse() && sa.strand == "-")
            {
                goto findInvFirst;
            }

            continue;

        findInvFirst:

            if (readparser->hasFirstCigarSoftclipped())
            {
                if (sa.cigar.at(0).getOperatorName() == 'S' && sa.cigar.at(sa.cigar.size() - 1).getOperatorName() == 'M')
                {
                    mapINV[std::make_pair(readparser->getPos(), sa.pos)].NumberOfMatchRead++;
                    mapINV[std::make_pair(readparser->getEnd(), sa.pos)].MatchLists.push_back(sa.cigar.at(sa.cigar.size() - 1).getLength());
                    mapINV[std::make_pair(readparser->getEnd(), sa.pos)].MapQLists.push_back(readparser->getMapQuality());
                }
            }

            if (readparser->hasLastCigarSoftclipped())
            {
                if (sa.cigar.at(0).getOperatorName() == 'M' && sa.cigar.at(sa.cigar.size() - 1).getOperatorName() == 'S')
                {
                    mapINV[std::make_pair(readparser->getEnd(), sa.pos+sa.cigar.at(0).getLength())].NumberOfMatchRead++;
                    mapINV[std::make_pair(readparser->getEnd(), sa.pos+sa.cigar.at(0).getLength())].MatchLists.push_back(sa.cigar.at(0).getLength());
                    mapINV[std::make_pair(readparser->getEnd(), sa.pos+sa.cigar.at(0).getLength())].MapQLists.push_back(readparser->getMapQuality());
                }
            }
        }
        else
        {

            if (readparser->isReverse() && sa.strand == "+")
            {
                goto findInvSecond;
            }

            if (!readparser->isReverse() && sa.strand == "-")
            {
                goto findInvSecond;
            }

            continue;

        findInvSecond:

            if (readparser->hasFirstCigarSoftclipped())
            {
                if (sa.cigar.at(0).getOperatorName() == 'S' && sa.cigar.at(sa.cigar.size() - 1).getOperatorName() == 'M')
                {
                    mapINV[std::make_pair(sa.pos, readparser->getPos())].NumberOfMatchRead++;
                    mapINV[std::make_pair(sa.pos, readparser->getPos())].MatchLists.push_back(sa.cigar.at(sa.cigar.size() - 1).getLength());
                    mapINV[std::make_pair(sa.pos, readparser->getPos())].MapQLists.push_back(readparser->getMapQuality());
                }
            }

            if (readparser->hasLastCigarSoftclipped())
            {
                if (sa.cigar.at(0).getOperatorName() == 'M' && sa.cigar.at(sa.cigar.size() - 1).getOperatorName() == 'S')
                {
                    mapINV[std::make_pair(sa.pos+sa.cigar.at(0).getLength(), readparser->getEnd())].NumberOfMatchRead++;
                    mapINV[std::make_pair(sa.pos+sa.cigar.at(0).getLength(), readparser->getEnd())].MatchLists.push_back(sa.cigar.at(0).getLength());
                    mapINV[std::make_pair(sa.pos+sa.cigar.at(0).getLength(), readparser->getEnd())].MapQLists.push_back(readparser->getMapQuality());
                }
            }
        }
    }
}

void SplitRead::findDeletion()
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
            // if (sa.pos - readparser->getPos() > (samplestate->getReadLength()*1.5))
            // {
            //     continue;
            // }

            if (readparser->isReverse() && sa.strand == "-")
            {
                goto findDelFirst;
            }

            if (!readparser->isReverse() && sa.strand == "+")
            {
                goto findDelFirst;
            }

            continue;

        findDelFirst:

            if (!(readparser->hasLastCigarSoftclipped()))
            {
                continue;
            }

            if (sa.cigar.at(0).getOperatorName() == 'S' && sa.cigar.at(sa.cigar.size() - 1).getOperatorName() == 'M')
            {
                mapDEL[std::make_pair(readparser->getEnd(), sa.pos)].NumberOfMatchRead++;
                mapDEL[std::make_pair(readparser->getEnd(), sa.pos)].MatchLists.push_back(sa.cigar.at(sa.cigar.size() - 1).getLength());
                mapDEL[std::make_pair(readparser->getEnd(), sa.pos)].MapQLists.push_back(readparser->getMapQuality());
            }
        }
        else
        {

            if (readparser->isReverse() && sa.strand == "-")
            {
                goto findDelSecond;
            }

            if (!readparser->isReverse() && sa.strand == "+")
            {
                goto findDelSecond;
            }

            continue;

        findDelSecond:

            if (!(readparser->hasFirstCigarSoftclipped()))
            {
                continue;
            }

            if (sa.cigar.at(0).getOperatorName() == 'M' && sa.cigar.at(sa.cigar.size() - 1).getOperatorName() == 'S')
            {
                mapDEL[std::make_pair(sa.pos, readparser->getPos())].NumberOfMatchRead++;
                mapDEL[std::make_pair(sa.pos, readparser->getPos())].MatchLists.push_back(sa.cigar.at(0).getLength());
                mapDEL[std::make_pair(sa.pos, readparser->getPos())].MapQLists.push_back(readparser->getMapQuality());
            }
        }
    }
}

void SplitRead::findTandemDuplication()
{

    // if (satag.size() != 1)
    // {
    //     return;
    // }

    for (ReadParser::SATag sa : satag)
    {
        // if (sa.mapQ < 30)
        // {
        //     continue;
        // }

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
            if (sa.pos - readparser->getPos() > (samplestate->getReadLength() * 1.5))
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
                mapDUP[std::make_pair(readparser->getPos(), sa.pos + sa.cigar.at(0).getLength())].MatchLists.push_back(sa.cigar.at(0).getLength());
                mapDUP[std::make_pair(readparser->getPos(), sa.pos + sa.cigar.at(0).getLength())].MapQLists.push_back(readparser->getMapQuality());
            }
        }
        else
        {

            if (readparser->getPos() - sa.pos > (samplestate->getReadLength() * 1.5))
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

            if (sa.cigar.at(0).getOperatorName() == 'S' && sa.cigar.at(sa.cigar.size() - 1).getOperatorName() == 'M')
            {
                mapDUP[std::make_pair(sa.pos, readparser->getEnd())].NumberOfMatchRead++;
                mapDUP[std::make_pair(sa.pos, readparser->getEnd())].MatchLists.push_back(sa.cigar.at(sa.cigar.size() - 1).getLength());
                mapDUP[std::make_pair(sa.pos, readparser->getEnd())].MapQLists.push_back(readparser->getMapQuality());
            }
        }
    }
}

void SplitRead::printResult()
{
    // printDeletion();
    // printDuplication();
    printInversion();
}

void SplitRead::printDeletion()
{
    std::cout << "Deletion" << std::endl;
    for (auto const &x : mapDEL)
    {
        if (x.second.NumberOfMatchRead >= 2)
        {
            Evidence evidence;
            evidence.setPos(x.first.first);
            evidence.setEnd(x.first.second);
            evidence.setChr(chrname);
            evidence.setEndChr(chrname);
            evidence.setFrequency(x.second.NumberOfMatchRead);
            evidence.setVariantType("DEL");
            evidence.setMark("SR");
            evidence.setMapQList(x.second.MapQLists);

            if (evidence.getEnd() - evidence.getPos() < 100)
            {
                continue;
            }

            if (evidence.getMaxMapQ() < 30)
            {
                continue;
            }

            if (evidence.getMinMapQ() == 0)
            {
                continue;
            }

            vecDEL.push_back(evidence);

            // std::cout
            //     << " pos : " << x.first.first
            //     << " end : " << x.first.second
            //     << " NumberOfMatchRead : " << x.second.NumberOfMatchRead
            //     << std::endl;
        }
    }

    removeDuplicateResult(&vecDEL);
    for (auto x : vecDEL)
    {
        writeFile(x);
    }
}

void SplitRead::printDuplication()
{
    //   std::cout << "duplication" << std::endl;
    for (auto const &x : mapDUP)
    {
        if (x.second.NumberOfMatchRead >= 2)
        {
            Evidence evidence;
            evidence.setPos(x.first.first);
            evidence.setEnd(x.first.second);
            evidence.setChr(chrname);
            evidence.setEndChr(chrname);
            evidence.setFrequency(x.second.NumberOfMatchRead);
            evidence.setVariantType("DUP");
            evidence.setMark("SR");
            // evidence.LNGMATCH = x.second.
            evidence.setMapQList(x.second.MapQLists);

            if (evidence.getEnd() - evidence.getPos() < 100)
            {
                continue;
            }

            if (evidence.getMaxMapQ() < 50)
            {
                continue;
            }

            if (evidence.getMinMapQ() == 0)
            {
                continue;
            }

            vecDUP.push_back(evidence);

            // std::cout
            //     << " pos : " << x.first.first
            //     << " end : " << x.first.second
            //     << " NumberOfMatchRead : " << x.second.NumberOfMatchRead
            //     << std::endl;
        }
    }

    // removeDuplicateResult(&vecDUP);
    for (auto x : vecDUP)
    {
        writeFile(x);
    }
}


void SplitRead::printInversion()
{
    for (auto const &x : mapINV)
    {
        if (x.second.NumberOfMatchRead >= 1)
        {
            Evidence evidence;
            evidence.setPos(x.first.first);
            evidence.setEnd(x.first.second);
            evidence.setChr(chrname);
            evidence.setEndChr(chrname);
            evidence.setFrequency(x.second.NumberOfMatchRead);
            evidence.setVariantType("INV");
            evidence.setMark("SR");
            // evidence.LNGMATCH = x.second.
            evidence.setMapQList(x.second.MapQLists);

            if (evidence.getEnd() - evidence.getPos() < 100)
            {
                continue;
            }

            if (evidence.getMaxMapQ() < 50)
            {
                continue;
            }

            if (evidence.getMinMapQ() == 0)
            {
                continue;
            }

            vecINV.push_back(evidence);

            // std::cout
            //     << " pos : " << x.first.first
            //     << " end : " << x.first.second
            //     << " NumberOfMatchRead : " << x.second.NumberOfMatchRead
            //     << std::endl;
        }
    }

    // removeDuplicateResult(&vecDUP);
    for (auto x : vecINV)
    {
        writeFile(x);
    }
}

void SplitRead::removeDuplicateResult(std::vector<Evidence> *vec)
{
    std::vector<Evidence> temp;
    bool added = false;
    for (auto n : *vec)
    {
        added = false;
        for (auto m : temp)
        {
            if (n.getPos() == m.getPos() && n.getEnd() == m.getEnd())
            {
                continue;
            }

            if (checkBetween(n.getPos(), m.getPos(), 200) && checkBetween(n.getEnd(), m.getEnd(), 200))
            {
                added = true;
                break;
            }
        }

        if (!added)
        {
            temp.push_back(n);
        }
    }

    *vec = temp;
}



bool SplitRead::checkBetween(int32_t pos, int32_t targetPos, int32_t overlapped)
{
    if (targetPos - overlapped > pos)
    {
        return false;
    }

    if (targetPos + overlapped < pos)
    {
        return false;
    }

    return true;
}

int SplitRead::writeFile(Evidence vr)
{
    std::ofstream myfile;
    myfile.open(filepath->getOutputPath() + "/analysis/splitread/" + vr.getChr() + "." + vr.getVariantType() + ".txt", std::ios_base::app);
    // std::cout << filepath->getOutputPath() + "/analysis/splitread/" + vr.getChr() + "." + vr.getVariantType() + ".txt" << std::endl;
    vr.setID("BOLT" + std::to_string(vcfIdNumber));
    myfile << vr.getResultVcfFormatString() << std::endl;
    vcfIdNumber++;
    myfile.close();
    return 0;
}