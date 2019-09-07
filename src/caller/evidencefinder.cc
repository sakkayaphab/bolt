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
    while (sam_itr_next(inT, iterT, read) >= 0)
    {
        if ((read->core.flag & BAM_FUNMAP))
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
        
        std::vector<ReadParser::SATag> satag = readparser.getSATag();
        if (satag.size()!=0) {
            for (auto n:satag) {
                std::cout << 
                n.chrname << 
                " " <<
                n.pos <<
                " " <<
                n.strand << std::endl;
            }
        }

        currentPos = read->core.pos + 1;
        currentMPos = read->core.mpos + 1;

        cigar = readparser.getCigar();
        // if (readparser.getMapQuality() > 0)
        // {
            seDeletion.updateRead();
            seInsertion.updateRead();
            seInversion.updateRead();
            seTandemDuplication.updateRead();
            seTranslocation.updateRead();
        // }

        // std::cout << readparser.getPosDiscordantRead() << ",preCollectDEL : " << preCollectDEL.size() << ",collectDeletionInfoLists : " << collectDeletionInfoLists.size() << std::endl;

        //Read depth
        getRound(&currentPos, &configRound, &roundedPos);
        if (roundedCurrentPos == roundedPos)
        {
            readdepthdetail.RD++;

            if (cigar.at(0).getOperatorName() == 'S' && cigar.at(0).getLength() >= 2)
            {
                readdepthdetail.SCF++;
            }

            if (cigar.at(cigar.size() - 1).getOperatorName() == 'S' && cigar.at(0).getLength() >= 2)
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

    // std::cout << "---------------------------------------" << std::endl;
    // std::cout << "✓ " << *target_chromosome << std::endl;
    // std::cout << "number of read : " << countT << std::endl;
    // std::cout << "number of region deletion : " << countDEL << std::endl;
    // std::cout << "size of read depth line segment : " << ReadDepthLineSegment.size() << std::endl;
    // std::cout << "---------------------------------------" << std::endl;

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
    rdh.writeReadDepthStat(filepath->getReadDepthStatPath() + "/" + "readdepthstat.txt");
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

void EvidenceFinder::checkNormalRead(ReadDepthDetail *rdd) {
    if (readparser.isFirstRead()) {
        if (readparser.isMateUnmapped()) {
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
        // if (insertSizeFirstRead>samplestat->getMedianSampleStat() + samplestat->getSDSampleStat() + samplestat->getReadLength() + 200) {
        //     // rdd->ABN_READ++;
        //     return;
        // }

        // if (insertSizeFirstRead < samplestat->getMedianSampleStat() - samplestat->getSDSampleStat() - samplestat->getReadLength() - 200)
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
    } else {
         if (readparser.isMateUnmapped()) {
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
        // if (insertSizeSecondRead < samplestat->getMedianSampleStat() - samplestat->getSDSampleStat() - samplestat->getReadLength() - 200) {
        //     // rdd->ABN_READ++;
        //     return;
        // }

        // if (insertSizeSecondRead > samplestat->getMedianSampleStat() + samplestat->getSDSampleStat() + samplestat->getReadLength() + 200) {
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

    if (readparser.isSecondRead())
    {
        insertSizeSecondRead = (readparser.getPos() + readparser.getLengthSequence()) - readparser.getMatePos();

        if (readparser.isMateUnmapped())
        {
            rdd->INS2++;
            return;
        }

        if (!readparser.isPairOnSameChromosome())
        {
            rdd->TRA2++;
            return;
        }

        if (insertSizeSecondRead < 0)
        {
            return;
        }

        if (readparser.isReverse() && readparser.isMateReverse())
        {
            rdd->INV2++;
            return;
        }

        if (!readparser.isReverse() && !readparser.isMateReverse())
        {
            rdd->INV2++;
            return;
        }

        if (insertSizeSecondRead > samplestat->getMedianSampleStat() + samplestat->getSDSampleStat() + samplestat->getReadLength() + 200)
        {
            rdd->DEL2++;
            return;
        }

        if (insertSizeSecondRead < samplestat->getMedianSampleStat() - samplestat->getSDSampleStat() - samplestat->getReadLength() - 120)
        {
            rdd->INS2++;
            return;
        }

        if (!readparser.isReverse())
        {
            return;
        }

        if (readparser.isMateReverse())
        {
            return;
        }

        if (readparser.getPos() < readparser.getMatePos())
        {
            rdd->DUP2++;
            return;
        }
        return;
    }

    if (readparser.isFirstRead())
    {
        insertSizeFirstRead = (readparser.getMatePos() + readparser.getLengthSequence()) - readparser.getPos();
        if (readparser.isMateUnmapped())
        {
            rdd->INS1++;
            return;
        }

        if (!readparser.isPairOnSameChromosome())
        {
            rdd->TRA1++;
            return;
        }

        if (insertSizeFirstRead < 0)
        {
            return;
        }

        if (readparser.isReverse() && readparser.isMateReverse())
        {
            rdd->INV1++;
            return;
        }

        if (!readparser.isReverse() && !readparser.isMateReverse())
        {
            rdd->INV1++;
            return;
        }

        if (insertSizeFirstRead > samplestat->getMedianSampleStat() + samplestat->getSDSampleStat() + samplestat->getReadLength() + 200)
        {
            rdd->DEL1++;
            return;
        }

        if (insertSizeFirstRead < samplestat->getMedianSampleStat() - samplestat->getSDSampleStat() - samplestat->getReadLength() - 120)
        {
            rdd->INS1++;
            return;
        }

        if (readparser.isReverse())
        {
            return;
        }

        if (!readparser.isMateReverse())
        {
            return;
        }

        if (readparser.getPos() > readparser.getMatePos())
        {
            rdd->DUP1++;
            return;
        }
        return;
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