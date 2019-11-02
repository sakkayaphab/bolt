#include "splitread.h"
#include <iostream>
#include <fstream>
SplitRead::SplitRead(std::string chrname, ReadParser *readparser, SampleStat *samplestat, FileManager *filepath)
{
    SplitRead::readparser = readparser;
    SplitRead::samplestat = samplestat;
    SplitRead::chrname = chrname;
    SplitRead::filepath = filepath;
}

void SplitRead::updateRead()
{
    findDeletionInRead();
    findInsertionInRead();
    satag = readparser->getSATag();
    if (satag.size() == 0)
    {
        return;
    }

    findInversion();
    findDeletion();
    findInsertion();
    findTandemDuplication();
}

void SplitRead::findInversion()
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
                    mapINV[std::make_pair(readparser->getEnd(), sa.pos + sa.cigar.at(0).getLength())].NumberOfMatchRead++;
                    mapINV[std::make_pair(readparser->getEnd(), sa.pos + sa.cigar.at(0).getLength())].MatchLists.push_back(sa.cigar.at(0).getLength());
                    mapINV[std::make_pair(readparser->getEnd(), sa.pos + sa.cigar.at(0).getLength())].MapQLists.push_back(readparser->getMapQuality());
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
                    mapINV[std::make_pair(sa.pos + sa.cigar.at(0).getLength(), readparser->getEnd())].NumberOfMatchRead++;
                    mapINV[std::make_pair(sa.pos + sa.cigar.at(0).getLength(), readparser->getEnd())].MatchLists.push_back(sa.cigar.at(0).getLength());
                    mapINV[std::make_pair(sa.pos + sa.cigar.at(0).getLength(), readparser->getEnd())].MapQLists.push_back(readparser->getMapQuality());
                }
            }
        }
    }
}

void SplitRead::findInsertionInRead()
{
    std::vector<ReadParser::Cigar> cigar = readparser->getCigar();

    if (cigar.size() <= 2)
    {
        return;
    }

    // CIGAR = 20S71M6D43M54D117M = 20 + 71(M) + 43(M) + 117(M)
    int32_t incrementPos = 0;

    for (int i = 0; i < cigar.size(); i++)
    {
        if (cigar.at(i).getOperatorName() == 'M')
        {
            incrementPos += cigar.at(i).getLength();
        }

        if (cigar.at(i).getOperatorName() == 'I')
        {

            if (cigar.at(i).getLength() > 10)
            {

                mapSmallINS[std::make_pair(readparser->getPos() + incrementPos, readparser->getPos() + incrementPos + cigar.at(i).getLength())].NumberOfMatchRead++;
                mapSmallINS[std::make_pair(readparser->getPos() + incrementPos, readparser->getPos() + incrementPos + cigar.at(i).getLength())].MatchLists.push_back(cigar.at(i).getLength());
                mapSmallINS[std::make_pair(readparser->getPos() + incrementPos, readparser->getPos() + incrementPos + cigar.at(i).getLength())].MapQLists.push_back(readparser->getMapQuality());

                // std::cout << "D" << cigar.at(i).getLength() << " "
                // << readparser->getPos() + incrementPos << " = " << readparser->getPos() + incrementPos + cigar.at(i).getLength()
                // << " " << mapSmallDEL[std::make_pair(readparser->getPos() + incrementPos, readparser->getPos() + incrementPos + cigar.at(i).getLength())].NumberOfMatchRead
                // << std::endl;
            }

            incrementPos += cigar.at(i).getLength();
        }
    }
}

void SplitRead::findDeletionInRead()
{
    std::vector<ReadParser::Cigar> cigar = readparser->getCigar();

    if (cigar.size() <= 2)
    {
        return;
    }

    // CIGAR = 20S71M6D43M54D117M = 20 + 71(M) + 43(M) + 117(M)
    int32_t incrementPos = 0;

    for (int i = 0; i < cigar.size(); i++)
    {
        if (cigar.at(i).getOperatorName() == 'M')
        {
            incrementPos += cigar.at(i).getLength();
        }

        if (cigar.at(i).getOperatorName() == 'D')
        {

            if (cigar.at(i).getLength() > 10)
            {

                mapSmallDEL[std::make_pair(readparser->getPos() + incrementPos, readparser->getPos() + incrementPos + cigar.at(i).getLength())].NumberOfMatchRead++;
                mapSmallDEL[std::make_pair(readparser->getPos() + incrementPos, readparser->getPos() + incrementPos + cigar.at(i).getLength())].MatchLists.push_back(cigar.at(i).getLength());
                mapSmallDEL[std::make_pair(readparser->getPos() + incrementPos, readparser->getPos() + incrementPos + cigar.at(i).getLength())].MapQLists.push_back(readparser->getMapQuality());

                // std::cout << "D" << cigar.at(i).getLength() << " "
                // << readparser->getPos() + incrementPos << " = " << readparser->getPos() + incrementPos + cigar.at(i).getLength()
                // << " " << mapSmallDEL[std::make_pair(readparser->getPos() + incrementPos, readparser->getPos() + incrementPos + cigar.at(i).getLength())].NumberOfMatchRead
                // << std::endl;
            }

            incrementPos += cigar.at(i).getLength();
        }
    }
}

bool SplitRead::haveSmallDeletion()
{
    // std::vector<ReadParser::Cigar> cigar = readparser.getCigar();
    // if (cigar.size() <= 1)
    // {
    //     return false;
    // }

    // int32_t posDel = 0;
    // int32_t indelpos = 0;
    // int32_t indelend = 0;

    // for (auto c : cigar)
    // {
    //     if (c.getOperatorName() == 'M')
    //     {
    //         posDel += c.getLength();
    //     }

    //     if (c.getOperatorName() == 'D')
    //     {
    //         indelpos = posDel + readparser->getPos();
    //         indelend = posDel + readparser->getPos() + c.getLength();
    //     }
    // }

    // if (indelend-indelpos<50) {
    //     continue;
    // }

    // if (indelend-indelpos>500) {
    //     continue;
    // }

    // if (indelpos != 0 && indelend != 0)
    // {
    //     mapSmallDEL[std::make_pair(indelpos, indelend)].NumberOfMatchRead++;
    //     mapSmallDEL[std::make_pair(indelpos, indelend)].MatchLists.push_back(indelend-indelpos);
    //     mapSmallDEL[std::make_pair(indelpos, indelend)].MapQLists.push_back(readparser->getMapQuality());
    // }
}

void SplitRead::findInsertion()
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

        if (readparser->isReverse() && sa.strand == "-")
        {
            goto findInsertion;
        }

        if (!readparser->isReverse() && sa.strand == "+")
        {
            goto findInsertion;
        }

    findInsertion:

        // if (readparser->hasFirstCigarSoftclipped() && sa.cigar.at(0).getOperatorName() == 'S')
        // {

        //     if (checkBetween(readparser->getPos(), sa.pos, getDivider(samplestat->getReadLength(), 1, 10, 1)))
        //     {
        //         mapINS[std::make_pair(readparser->getPos(), sa.pos)].NumberOfMatchRead++;
        //         mapINS[std::make_pair(readparser->getPos(), sa.pos)].MatchLists.push_back(sa.cigar.at(sa.cigar.size() - 1).getLength());
        //         mapINS[std::make_pair(readparser->getPos(), sa.pos)].MapQLists.push_back(readparser->getMapQuality());
        //     }
        // }

        if (readparser->hasLastCigarSoftclipped() && sa.cigar.at(0).getOperatorName() == 'S')
        {
            // std::cout << "findInsertion SA" << std::endl;
            int softclip = readparser->getSoftClippedSequenceEnd().size();
            int match = 0;
            if (sa.cigar.at(sa.cigar.size() - 1).getOperatorName() == 'M')
            {
                match = sa.cigar.at(sa.cigar.size() - 1).getLength();
            }

            if (match == 0)
            {
                continue;
            }

            int lsize = softclip - match;

            if (lsize<50)
            {
                continue;
            }

            if (checkBetween(readparser->getEnd(), sa.pos, getDivider(samplestat->getReadLength(), 1, 10, 1)))
            {
                std::cout << readparser->getEnd() << " " << sa.pos << std::endl;
                mapINS[std::make_pair(readparser->getEnd(), sa.pos)].NumberOfMatchRead++;
                mapINS[std::make_pair(readparser->getEnd(), sa.pos)].MatchLists.push_back(lsize);
                mapINS[std::make_pair(readparser->getEnd(), sa.pos)].MapQLists.push_back(readparser->getMapQuality());
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

                // std::cout << readparser->getEnd() << " = " << sa.pos << std::endl;
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
            // if (sa.pos - readparser->getPos() > (samplestat->getReadLength() * 1.5))
            // {
            //     continue;
            // }

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

            // if (readparser->getPos() - sa.pos > (samplestat->getReadLength() * 1.5))
            // {
            //     continue;
            // }

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
    printDeletion();
    printSmallDeletion();
    printSmallInsertion();
    printDuplication();
    printInversion();
    printInsertion();
}

std::vector<Evidence> SplitRead::convertMapToEvidenceList(std::map<std::pair<int32_t, int32_t>, RefiningSV::MatchRead> *mapSV, std::string svtype, std::string mark)
{
    std::vector<Evidence> vecTemp;
    for (auto const &x : *mapSV)
    {

        if (x.second.MapQLists.size() >= 1)
        {
            Evidence evidence;
            evidence.setPos(x.first.first);
            evidence.setEnd(x.first.second);
            evidence.setChr(chrname);
            evidence.setEndChr(chrname);
            evidence.setFrequency(x.second.MapQLists.size());
            evidence.setVariantType(svtype);
            evidence.setMark(mark);
            evidence.setMapQList(x.second.MapQLists);

            vecTemp.push_back(evidence);
        }
    }

    return vecTemp;
}

void SplitRead::mergeEvidence(std::vector<Evidence> *vecTemp)
{
    std::vector<Evidence> newEvidenceTempList;
    std::sort(vecTemp->begin(), vecTemp->end());

    for (auto m : *vecTemp)
    {
        bool incremented = false;
        for (int32_t i = 0; i < newEvidenceTempList.size(); i++)
        {
            if (checkBetween(m.getPos(), newEvidenceTempList.at(i).getPos(), 2) && checkBetween(m.getEnd(), newEvidenceTempList.at(i).getEnd(), 2))
            {
                incremented = true;
                for (auto x : *m.getMapQVector())
                {
                    newEvidenceTempList.at(i).addMapQ(x);
                    newEvidenceTempList.at(i).incrementFrequency();
                }
                break;
            }
        }

        if (!incremented)
        {
            newEvidenceTempList.push_back(m);
        }
    }

    *vecTemp = newEvidenceTempList;
}

void SplitRead::setAllCIEvidence(std::vector<Evidence> *elist, int32_t rangePos)
{
    for (int32_t i = 0; i < elist->size(); i++)
    {
        elist->at(i).setCiPosLeft(-rangePos);
        elist->at(i).setCiPosRight(rangePos);
        elist->at(i).setCiEndLeft(-rangePos);
        elist->at(i).setCiEndRight(rangePos);
    }
}

void SplitRead::filterFrequencyLowerThan(int number, std::vector<Evidence> *elist)
{
    std::vector<Evidence> newEvidenceTempList;

    for (auto n : *elist)
    {
        if (n.getFrequency() <= number)
        {
            continue;
        }

        newEvidenceTempList.push_back(n);
    }

    *elist = newEvidenceTempList;
}

void SplitRead::filterMapQLowerThan(uint8_t mapq, std::vector<Evidence> *elist)
{
    std::vector<Evidence> newEvidenceTempList;

    for (auto n : *elist)
    {
        if (n.getMaxMapQ() < mapq)
        {
            continue;
        }

        newEvidenceTempList.push_back(n);
    }

    *elist = newEvidenceTempList;
}

void SplitRead::filterEvidenceList(std::vector<Evidence> *elist)
{
    std::vector<Evidence> newEvidenceTempList;

    for (auto n : *elist)
    {
        if (n.getMaxMapQ() < 60)
        {
            continue;
        }

        if (n.getMinMapQ() == 0)
        {
            continue;
        }

        newEvidenceTempList.push_back(n);
    }

    *elist = newEvidenceTempList;
}

void SplitRead::filterLengthMinEvidenceList(std::vector<Evidence> *elist, int32_t min)
{
    std::vector<Evidence> newEvidenceTempList;

    for (auto n : *elist)
    {
        if (n.getSvLength() < min)
        {
            continue;
        }

        newEvidenceTempList.push_back(n);
    }

    *elist = newEvidenceTempList;
}

void SplitRead::filterLengthMaxEvidenceList(std::vector<Evidence> *elist, int32_t max)
{
    std::vector<Evidence> newEvidenceTempList;

    for (auto n : *elist)
    {
        if (n.getSvLength() > max)
        {
            continue;
        }

        newEvidenceTempList.push_back(n);
    }

    *elist = newEvidenceTempList;
}

void SplitRead::printDuplication()
{
    auto vecTemp = convertMapToEvidenceList(&mapDUP, "DUP", "SR");
    mergeEvidence(&vecTemp);
    setAllCIEvidence(&vecTemp, samplestat->getReadLength());
    filterEvidenceList(&vecTemp);
    // filterLengthMinEvidenceList(&vecTemp, 100);
    filterLengthMaxEvidenceList(&vecTemp, samplestat->getReadLength() * 2);

    for (auto x : vecTemp)
    {
        writeFile(x);
    }
}

void SplitRead::printInversion()
{
    auto vecTemp = convertMapToEvidenceList(&mapINV, "INV", "SR");
    mergeEvidence(&vecTemp);
    setAllCIEvidence(&vecTemp, samplestat->getReadLength() * 2);
    filterEvidenceList(&vecTemp);
    // filterLengthMinEvidenceList(&vecTemp, 50);
    filterLengthMaxEvidenceList(&vecTemp, 1000000);
    filterFrequencyLowerThan(1, &vecTemp);
    // filterMapQLowerThan(60, &vecTemp);

    for (auto x : vecTemp)
    {
        writeFile(x);
    }
}

void SplitRead::printInsertion()
{
    auto vecTemp = convertMapToEvidenceList(&mapINS, "INS", "SR");
    // filterLengthMinEvidenceList(&vecTemp, 49);
    // setAllCIEvidence(&vecTemp, samplestat->getReadLength() * 2);
    filterFrequencyLowerThan(0, &vecTemp);
    for (auto x : vecTemp)
    {
        writeFile(x);
    }
}

void SplitRead::printSmallInsertion()
{
    auto vecTemp = convertMapToEvidenceList(&mapSmallINS, "INS", "SINS");
    // mergeEvidence(&vecTemp);

    filterLengthMinEvidenceList(&vecTemp, 49);
    filterFrequencyLowerThan(1, &vecTemp);

    for (auto x : vecTemp)
    {
        writeFile(x);
    }
}

void SplitRead::printSmallDeletion()
{
    auto vecTemp = convertMapToEvidenceList(&mapSmallDEL, "DEL", "SDEL");
    // mergeEvidence(&vecTemp);
    filterFrequencyLowerThan(1, &vecTemp);

    for (auto x : vecTemp)
    {
        writeFile(x);
    }
}

void SplitRead::printDeletion()
{
    auto vecTemp = convertMapToEvidenceList(&mapDEL, "DEL", "SR");
    mergeEvidence(&vecTemp);
    setAllCIEvidence(&vecTemp, samplestat->getReadLength() * 2);
    filterEvidenceList(&vecTemp);
    // filterLengthMinEvidenceList(&vecTemp, 50);
    filterLengthMaxEvidenceList(&vecTemp, 1000000);
    filterFrequencyLowerThan(1, &vecTemp);

    for (auto x : vecTemp)
    {
        // std::cout << x.getResultVcfFormatString() << std::endl;

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

int SplitRead::getDivider(int value, int top, int down, int minimum)
{
    auto returnvalue = (int)(float(value) * (float(top) / float(down)));

    if (returnvalue > minimum)
    {
        return returnvalue;
    }

    return minimum;
}