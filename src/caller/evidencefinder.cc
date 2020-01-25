#include "evidencefinder.h"
#include <stdio.h>
#include <string>
#include <iostream>
#include <cstring>
#include <vector>
#include <iterator>
#include <algorithm>
#include "specifyingevidencedeletion.h"
#include "specifyingevidenceinsertion.h"
#include "specifyingevidenceinversion.h"
#include "specifyingevidencetandemduplication.h"
#include "specifyingevidencetranslocation.h"
#include <unistd.h>
#include "readdepthhelper.h"
#include "splitread.h"

EvidenceFinder::EvidenceFinder(SampleStat *samplestat_T, FileManager *filepath_T, std::string *target_chromosome_T)
{
    samplestat = samplestat_T;
    filepath = filepath_T;
    target_chromosome = target_chromosome_T;
}

void getRound(uint32_t *x, uint32_t *max, uint32_t *round)
{
    *round = (*x / *max) * *max;
}

void EvidenceFinder::setHtsIndex(hts_idx_t *index)
{
    bam_index = index;
}

void EvidenceFinder::findEvidence()
{
    hts_itr_t *iterT = NULL;
    samFile *inT = NULL;
    // bam_hdr_t *bam_header = NULL;
    inT = sam_open(filepath->getSamplePath().c_str(), "r");
    if (inT == NULL)
        return;

    // target_chromosome std::string to array char
    int n = target_chromosome->length();
    char char_array_chrTarget[n + 1];
    strcpy(char_array_chrTarget, target_chromosome->c_str());
    iterT = sam_itr_querys(bam_index, bam_header, char_array_chrTarget);

    read = bam_init1();
    setupReadParser(read);
    int countT = 0;

    ReadDepthDetail readdepthdetail;
    uint32_t coverage = 0;
    uint32_t configRound = 250;
    uint32_t roundedPos = 0;
    uint32_t roundedCurrentPos = 0;

    int svtype;

    SpecifyingEvidenceDeletion seDeletion;
    seDeletion.setSampleStat(samplestat);
    seDeletion.setRead(read, bam_header);
    seDeletion.setOutputPath(filepath->getTempEvidencePath() + "/" + *target_chromosome + ".DEL.txt");

    SpecifyingEvidenceInsertion seInsertion;
    seInsertion.setSampleStat(samplestat);
    seInsertion.setRead(read, bam_header);
    seInsertion.setOutputPath(filepath->getTempEvidencePath() + "/" + *target_chromosome + ".INS.txt");

    SpecifyingEvidenceInversion seInversion;
    seInversion.setSampleStat(samplestat);
    seInversion.setRead(read, bam_header);
    seInversion.setOutputPath(filepath->getTempEvidencePath() + "/" + *target_chromosome + ".INV.txt");

    SpecifyingEvidenceTandemDuplication seTandemDuplication;
    seTandemDuplication.setSampleStat(samplestat);
    seTandemDuplication.setRead(read, bam_header);
    seTandemDuplication.setOutputPath(filepath->getTempEvidencePath() + "/" + *target_chromosome + ".DUP.txt");

    SpecifyingEvidenceTranslocation seTranslocation;
    seTranslocation.setSampleStat(samplestat);
    seTranslocation.setRead(read, bam_header);
    seTranslocation.setOutputPath(filepath->getTempEvidencePath() + "/" + *target_chromosome + ".TRA.txt");

    std::vector<ReadParser::Cigar> cigar;
    SplitRead splitread(*target_chromosome, &readparser, samplestat, filepath);
    int32_t markduppos = 0;
    int32_t markdupend = 0;
    int32_t markdupmatepos = 0;
    int32_t markdupmateend = 0;
    uint8_t markdupmapq = 0;

    while (sam_itr_next(inT, iterT, read) >= 0)
    {

        if (readparser.isUnmapped())
        {
            continue;
        }

        if (readparser.isNotPassingFilters())
        {
            continue;
        }

        if (readparser.isPCR())
        {
            continue;
        }

        if (readparser.isSupplementaryAlignment())
        {
            continue;
        }

        if (markduppos==readparser.getPos() && markdupend==readparser.getEnd()) {
            if (markdupmatepos==readparser.getMatePos() && markdupmateend==readparser.getMateEnd()) {
                if (markdupmapq == readparser.getMapQuality()) {
                    continue;
                }
            }

        }

        markduppos = readparser.getPos();
        markdupend = readparser.getEnd();
        markdupmatepos = readparser.getMatePos();
        markdupmateend = readparser.getMateEnd();
        markdupmapq = readparser.getMapQuality();


        splitread.updateRead();

        currentPos = read->core.pos + 1;
        currentMPos = read->core.mpos + 1;

        cigar = readparser.getCigar();

        seDeletion.updateRead();
        seInsertion.updateRead();
        seInversion.updateRead();
        seTandemDuplication.updateRead();
        seTranslocation.updateRead();

        //Read depth
        getRound(&currentPos, &configRound, &roundedPos);
        if (roundedCurrentPos == roundedPos)
        {
            readdepthdetail.RD++;

            if (cigar.at(0).getOperatorName() == 'S' && cigar.at(0).getLength() >= 8)
            {
                readdepthdetail.SCF++;
            }

            if (cigar.at(cigar.size() - 1).getOperatorName() == 'S' && cigar.at(cigar.size() - 1).getLength() >= 8)
            {
                readdepthdetail.SCL++;
            }

            // count sv
            updateReadDepthSV(&readdepthdetail);
        }
        else
        {
            coverage = readdepthdetail.RD * samplestat->getReadLength() / configRound;
            ReadDepthLineSegment[roundedCurrentPos].RD += coverage;
            ReadDepthLineSegment[roundedCurrentPos].DEL1 += readdepthdetail.DEL1;
            ReadDepthLineSegment[roundedCurrentPos].DUP1 += readdepthdetail.DUP1;
            ReadDepthLineSegment[roundedCurrentPos].INV1 += readdepthdetail.INV1;
            ReadDepthLineSegment[roundedCurrentPos].INS1 += readdepthdetail.INS1;
            ReadDepthLineSegment[roundedCurrentPos].TRA1 += readdepthdetail.TRA1;

            ReadDepthLineSegment[roundedCurrentPos].DEL2 += readdepthdetail.DEL2;
            ReadDepthLineSegment[roundedCurrentPos].DUP2 += readdepthdetail.DUP2;
            ReadDepthLineSegment[roundedCurrentPos].INV2 += readdepthdetail.INV2;
            ReadDepthLineSegment[roundedCurrentPos].INS2 += readdepthdetail.INS2;
            ReadDepthLineSegment[roundedCurrentPos].TRA2 += readdepthdetail.TRA2;

            ReadDepthLineSegment[roundedCurrentPos].SCF += readdepthdetail.SCF;
            ReadDepthLineSegment[roundedCurrentPos].SCL += readdepthdetail.SCL;

            ReadDepthLineSegment[roundedCurrentPos].R1_MUN += readdepthdetail.R1_MUN;
            ReadDepthLineSegment[roundedCurrentPos].R2_MUN += readdepthdetail.R2_MUN;

            readdepthdetail.RD = 0;
            coverage = 0;
            readdepthdetail.DEL1 = 0;
            readdepthdetail.DUP1 = 0;
            readdepthdetail.INS1 = 0;
            readdepthdetail.INV1 = 0;
            readdepthdetail.TRA1 = 0;
            readdepthdetail.DEL2 = 0;
            readdepthdetail.DUP2 = 0;
            readdepthdetail.INS2 = 0;
            readdepthdetail.INV2 = 0;
            readdepthdetail.TRA2 = 0;

            readdepthdetail.SCF = 0;
            readdepthdetail.SCL = 0;

            readdepthdetail.R1_MUN = 0;
            readdepthdetail.R2_MUN = 0;
            roundedCurrentPos = roundedPos;
        }
    }

    splitread.printResult();

    std::cout << "✓ " << *target_chromosome << std::endl;


    hts_itr_destroy(iterT);
    bam_destroy1(read);
    sam_close(inT);

    ReadDepthHelper rdh(filepath);
    rdh.setRange(250);
    rdh.setReadDepthMap(ReadDepthLineSegment);
    rdh.setTargetChromosome(*target_chromosome);
    rdh.writeReadDepthLineFile(filepath->getReadDepthPath() + "/" + *target_chromosome + ".txt");
    rdh.findEvidence();
    rdh.calculateAvgReaddepth();
    //    rdh.writeVcf(filepath->getTempEvidencePath() + "/" + *target_chromosome + ".RDDEL.txt");
    rdh.writeReadDepthStat(filepath->getReadDepthStatPath());
    //    if (rdh.getAvgReadDepth()>2000)
    //    {
    //        return;
    //    }

    seInsertion.setReadDepthHelper(&rdh);
    seInsertion.done();

    seDeletion.done();

    seInversion.done();

    seTandemDuplication.done();

    seTranslocation.done();
}

void EvidenceFinder::checkNormalRead(ReadDepthDetail *rdd)
{
    if (readparser.isFirstRead())
    {
        if (readparser.isMateUnmapped() && !readparser.isReverse())
        {
            rdd->R1_MUN++;
            return;
        }

        // if (readparser.isReverse()) {
        //     // rdd->ABN_READ++;
        //     return;
        // }

        // if (!readparser.isMateReverse()) {
        //     // rdd->ABN_READ++;
        //     return;
        // }

        // insertSizeFirstRead = (readparser.getMatePos() + readparser.getLengthSequence()) - readparser.getPos();
        // if (insertSizeFirstRead>samplestat->getAverageSampleStat() + samplestat->getSDSampleStat() + samplestat->getReadLength() + 200) {
        //     // rdd->ABN_READ++;
        //     return;
        // }

        // if (insertSizeFirstRead < samplestat->getAverageSampleStat() - samplestat->getSDSampleStat() - samplestat->getReadLength() - 200)
        // {
        //     // rdd->ABN_READ++;
        //     return;
        // }

        // if (!(readparser.getPos()<readparser.getMatePos())) {
        //     //  rdd->ABN_READ++;
        //     return;
        // }

        // rdd->N_READ++;
        return;
    }
    else
    {
        if (readparser.isMateUnmapped() && readparser.isReverse())
        {
            rdd->R2_MUN++;
            return;
        }

        // if (!readparser.isReverse()) {
        //     rdd->ABN_READ++;
        //     return;
        // }

        // if (readparser.isMateReverse()) {
        //     rdd->ABN_READ++;
        //     return;
        // }

        // insertSizeSecondRead = (readparser.getPos() + readparser.getLengthSequence()) - readparser.getMatePos();
        // if (insertSizeSecondRead < samplestat->getAverageSampleStat() - samplestat->getSDSampleStat() - samplestat->getReadLength() - 200) {
        //     // rdd->ABN_READ++;
        //     return;
        // }

        // if (insertSizeSecondRead > samplestat->getAverageSampleStat() + samplestat->getSDSampleStat() + samplestat->getReadLength() + 200) {
        //     // rdd->ABN_READ++;
        //     return;
        // }

        // if (!(readparser.getMatePos()<readparser.getPos())) {
        //     //  rdd->ABN_READ++;
        //     return;
        // }

        // rdd->N_READ++;
        return;
    }
}

void EvidenceFinder::updateReadDepthSV(ReadDepthDetail *rdd)
{
    //checkNormalRead
    checkNormalRead(rdd);

    if (readparser.isFirstRead())
    {
        if (isDeletion())
        {
            rdd->DEL1++;
        }

        if (isTandemDuplication())
        {
            rdd->DUP1++;
        }

        if (isInsertion())
        {
            rdd->INS1++;
        }

        if (isInversion())
        {
            rdd->INV1++;
        }

        if (isTranslocation())
        {
            rdd->TRA1++;
        }
    }

    if (readparser.isSecondRead())
    {
        if (isDeletion())
        {
            rdd->DEL2++;
        }

        if (isTandemDuplication())
        {
            rdd->DUP2++;
        }

        if (isInsertion())
        {
            rdd->INS2++;
        }

        if (isInversion())
        {
            rdd->INV2++;
        }

        if (isTranslocation())
        {
            rdd->TRA2++;
        }
    }

}

void EvidenceFinder::execute()
{
    findEvidence();
}

void EvidenceFinder::setSampleStat(SampleStat *samplestat_T)
{
    samplestat = samplestat_T;
}

void EvidenceFinder::setFilePath(FileManager *filepath_T)
{
    filepath = filepath_T;
}

void EvidenceFinder::setTargetChromosome(std::string *target_chromosome_T)
{
    target_chromosome = target_chromosome_T;
}

void EvidenceFinder::setBamHeader(bam_hdr_t *t_bam_header)
{
    bam_header = t_bam_header;
}

void EvidenceFinder::setupReadParser(bam1_t *alnT)
{
    read = alnT;
    readparser.setBamRead(alnT);
    readparser.setBamHeader(bam_header);
}

bool EvidenceFinder::isDefaultOrientation()
{
    if (readparser.isFirstRead())
    {
        if (readparser.isReverse())
        {
            return false;
        }

        if (!readparser.isMateReverse())
        {
            return false;
        }
    }

    if (readparser.isSecondRead())
    {
        if (!readparser.isReverse())
        {
            return false;
        }

        if (readparser.isMateReverse())
        {
            return false;
        }
    }

    return true;
}

bool EvidenceFinder::isDeletion()
{
    if (!readparser.isPairOnSameChromosome())
    {
        return false;
    }

    if (readparser.isUnmapped()) {
        return false;
    }

    if (readparser.isMateUnmapped()) {
        return false;
    }

    if (isTranslocation())
    {
        return false;
    }

    if (!isDefaultOrientation())
    {
        return false;
    }

    if (readparser.isFirstRead())
    {
        insertSizeFirstRead = (readparser.getMatePos() + readparser.getLengthSequence()) - readparser.getPos();

        if (insertSizeFirstRead > samplestat->getAverageSampleStat() + (samplestat->getSDSampleStat()*3))
        {
            return true;
        }
    }

    if (readparser.isSecondRead())
    {
        insertSizeSecondRead = (readparser.getPos() + readparser.getLengthSequence()) - readparser.getMatePos();

        if (insertSizeSecondRead > samplestat->getAverageSampleStat() + (samplestat->getSDSampleStat()*3))
        {
            return true;
        }
    }

    return false;
}

bool EvidenceFinder::isTandemDuplication()
{
    if (readparser.isUnmapped()) {
        return false;
    }

    if (readparser.isMateUnmapped()) {
        return false;
    }

    if (!readparser.isPairOnSameChromosome())
    {
        return false;
    }

    if (isTranslocation())
    {
        return false;
    }

    if (!isDefaultOrientation())
    {
        return false;
    }

    if (readparser.isFirstRead())
    {
        if (readparser.getPos() > readparser.getMatePos())
        {
            return true;
        }
    }

    if (readparser.isSecondRead())
    {
        if (readparser.getPos() < readparser.getMatePos())
        {
            return true;
        }
    }

    return false;
}

bool EvidenceFinder::isInversion()
{
    if (!readparser.isPairOnSameChromosome())
    {
        return false;
    }

    if (readparser.isUnmapped()) {
        return false;
    }

    if (readparser.isMateUnmapped()) {
        return false;
    }

    if (isTranslocation())
    {
        return false;
    }

    if (readparser.isReverse() && readparser.isMateReverse())
    {
        return true;
    }

    if (!readparser.isReverse() && !readparser.isMateReverse())
    {
        return true;
    }

    return false;
}

bool EvidenceFinder::isInsertion()
{
    
    if (readparser.isUnmapped()) {
        return false;
    }

    if (isTranslocation())
    {
        return false;
    }

    if (readparser.isFirstRead())
    {
        if (readparser.isReverse())
        {
            return false;
        }

        if (readparser.isMateUnmapped())
        {
            return true;
        }
    }

    if (readparser.isSecondRead())
    {
        if (!readparser.isReverse())
        {
            return false;
        }

        if (readparser.isMateUnmapped())
        {
            return true;
        }
    }

    if (!readparser.isPairOnSameChromosome())
    {
        return false;
    }

    if (!isDefaultOrientation())
    {
        return false;
    }

    if (readparser.isFirstRead())
    {
        insertSizeFirstRead = (readparser.getMatePos() + readparser.getLengthSequence()) - readparser.getPos();

        if (insertSizeFirstRead<0) {
            return false;
        }

        if (insertSizeFirstRead < samplestat->getAverageSampleStat() - (samplestat->getSDSampleStat()*2))
        {
            return true;
        }
    }

    if (readparser.isSecondRead())
    {
        insertSizeSecondRead = (readparser.getPos() + readparser.getLengthSequence()) - readparser.getMatePos();

        if (insertSizeSecondRead<0) {
            return false;
        }

        if (insertSizeSecondRead < samplestat->getAverageSampleStat() - (samplestat->getSDSampleStat()*2))
        {
            return true;
        }
        
    }

    
    return false;
}

bool EvidenceFinder::isTranslocation()
{
    if (readparser.isUnmapped()) {
        return false;
    }

    if (readparser.isMateUnmapped()) {
        return false;
    }

    if (!isDefaultOrientation())
    {
        return false;
    }

    if (readparser.isPairOnSameChromosome())
    {
        return false;
    }

    return true;
}