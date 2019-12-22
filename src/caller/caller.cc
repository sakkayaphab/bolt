#include "caller.h"
#include "samplestat.h"
#include <iostream>
#include <stdio.h>
#include <iomanip>
#include "task.h"
#include <string.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <cstdlib>
#include <iostream>
#include <unistd.h>
#include <fasta/fastareader.h>
#include "refiningtandemduplication.h"
#include "refiningtranslocation.h"
#include "refininginversion.h"
#include "refininginsertion.h"
#include "refiningdeletion.h"
#include "evidenceprovider.h"
#include <fstream>
#include <unistd.h>
#include <mutex>
#include <iostream>
#include <dirent.h>
#include "variantresultfilter.h"
#include <tbb/tbb.h>
#include "readdepthhelper.h"
#include "readdepthanalysis.h"
#include "refinedepthblock.h"
#include "depthblockfile.h"

Caller::~Caller()
{
    if (inFile != NULL)
    {
        sam_close(inFile);
    }

    if (bam_index != NULL)
    {
        hts_idx_destroy(bam_index);
    }
}

Caller::Caller(std::string samplepath_T, std::string referencepath_T, std::string outputpath_T)
{
    filepath.setOutputPath(outputpath_T);
    filepath.setReferencePath(referencepath_T);
    filepath.setSamplePath(samplepath_T);
    filepath.initialize();
    prepareHts();
    execSampleStat();
    applyBamHeader();

    // ReadDepthStat readDepthStat;
    // readDepthStat.setFilePath(&filepath);
    // readDepthStat.execute();
}

void Caller::setnumberofpair_stat(int n)
{
    numberofpair_stat = n;
}

void Caller::execSampleStat()
{
    samplestat.setSamplePath(filepath.getSamplePath());
    samplestat.execute();
}

void Caller::showinfo()
{
    std::cout << "+" << std::setfill('-') << std::setw(80) << "-" << std::endl;
    std::cout << "| Sample => " << filepath.getSamplePath() << std::endl;
    std::cout << "| Reference => " << filepath.getReferencePath() << std::endl;
    std::cout << "| Output => " << filepath.getOutputPath() << std::endl;
    std::cout << "+" << std::setfill('-') << std::setw(80) << "-" << std::endl;

    std::cout << "+" << std::setfill('-') << std::setw(80) << "-" << std::endl;
    std::cout << "| Read length => " << samplestat.getReadLength() << std::endl;
    std::cout << "| Average insert size => " << samplestat.getAverageSampleStat() << std::endl;
    std::cout << "| Standard deviation  => " << samplestat.getSDSampleStat() << std::endl;
    std::cout << "+" << std::setfill('-') << std::setw(80) << "-" << std::endl;
}

void Caller::prepareHts()
{
    inFile = sam_open(filepath.getSamplePath().c_str(), "r");
    if (inFile == NULL)
    {
        std::cout << "error : sam_open" << std::endl;
        exit(1);
        return;
    }

    bam_index = sam_index_load(inFile, filepath.getSamplePath().c_str());
    if (bam_index == NULL)
    {
        std::cout << "error : sam_index_load" << std::endl;
        return;
    }
}

void Caller::setParallel(int number)
{
    numberofparallel = number;
}

void Caller::applyBamHeader()
{
    bam_header = *getBamHeader();
}

bam_hdr_t *Caller::getBamHeader()
{
    samFile *inT = NULL;
    bam_hdr_t *headerT = NULL;

    inT = sam_open(filepath.getSamplePath().c_str(), "r");
    if (inT == NULL)
        return NULL;
    headerT = sam_hdr_read(inT);

    inT = sam_open(filepath.getSamplePath().c_str(), "r");
    if (inT == NULL)
        return NULL;
    if ((headerT = sam_hdr_read(inT)) == 0)
        return NULL;
    sam_close(inT);

    return headerT;
}

void Caller::execute()
{
    std::cout << std::endl;
    std::cout << "------------------------------" << std::endl;
    std::cout << "# Find evidence :" << std::endl;
    std::cout << "------------------------------" << std::endl;
    // find evidence by using threads

    // tbb::task_scheduler_init init(35);
    // std::mutex mxRead;
    // int32_t tasks = bam_header.n_targets;
    // //    std::cout << tasks << std::endl;
    // tbb::parallel_for(0, tasks, [&](int i) {
    //     //        mxRead.lock();
    //     i++;
    //     std::string tp(bam_header.target_name[i - 1]);
    //     Task task(samplestat, &filepath, tp);
    //     task.setHtsIndex(bam_index);
    //     task.setBamHeader(bam_header);
    //     //        mxRead.unlock();
    //     task.execute();
    //     task.showTaskDone();
    // });

    // std::cout << "start executing" << std::endl;
    // find evidence by using multiple process
    // int max_active = numberofparallel;
    int max_active = numberofparallel;
    int number_active = 0;
    bool done = false;
    int32_t tasks = bam_header.n_targets;
    int i = 0;
    for (; !done; ++number_active)
    {
        i++;
        std::string tp(bam_header.target_name[i - 1]);
        Task task(samplestat, &filepath, tp);
        task.setHtsIndex(bam_index);
        task.setBamHeader(bam_header);
        if (i >= tasks - 1)
        {
            done = true;
        }

        for (; number_active >= max_active; --number_active)
        {
            wait(NULL);
        }

        auto pid = fork();
        if (pid < 0)
        {
        }
        if (pid == 0)
        {
            // std::cout << "executing : " << tp << std::endl;
            task.execute();
            exit(0);
        }
        else
        {
            
        }
    }

    for (i = 0; i < max_active; i++)
    {
        wait(NULL);
    }
}

void Caller::mergeSplitRead()
{
    std::vector<std::string> evidenceFilePathLists;

    DIR *d;
    struct dirent *dir;
    d = opendir(filepath.getSplitReadPath().c_str());
    if (d)
    {
        while (dir = readdir(d))
        {

            if (std::string(dir->d_name).size() < 4)
            {
                continue;
            }

            if (std::string(dir->d_name).substr(std::string(dir->d_name).size() - 4) == ".txt")
            {
                std::string tempPath = filepath.getSplitReadPath() + "/" + std::string(dir->d_name);
                std::cout << tempPath << std::endl;
                evidenceFilePathLists.push_back(tempPath);
            }
        }
        closedir(d);
    }

    std::sort(evidenceFilePathLists.begin(), evidenceFilePathLists.end());

    //write file

    // int count = 0;
    // for (auto n : evidenceFilePathLists)
    // {
    //     std::vector<std::string> cache;

    //     std::string line;
    //     std::ifstream myfile(n);
    //     if (myfile.is_open())
    //     {

    //         while (getline(myfile, line))
    //         {
    //                 cache.push_back(line);
    //                 count++;

    //         }
    //         myfile.close();

    //         std::ofstream writefile;
    //         writefile.open(filepath.getOutputPath() + "/result.vcf", std::ios_base::app);
    //         for (auto a : cache)
    //         {
    //             writefile << a << std::endl;
    //         }
    //         writefile.close();

    //     }

    //     cache.clear();
    // }
}

void Caller::refineDelpthBlock()
{
    std::cout << std::endl;
    std::cout << "------------------------------" << std::endl;
    std::cout << "# Filter Breakpoint with Depth block : " << std::endl;
    std::cout << "------------------------------" << std::endl;

    RefineDepthBlock rdb;
    rdb.setFileManager(&filepath);
    rdb.setSampleStat(&samplestat);
    // rdb.filemanager = filepath;
    rdb.execute();
}

void Caller::catEvidenceFile()
{
    std::vector<std::string> evidenceFilePathLists;

    DIR *d;

    d = opendir(filepath.getTempEvidencePath().c_str());
    if (d)
    {
        struct dirent *dir;
        while (dir = readdir(d))
        {

            if (std::string(dir->d_name).size() < 4)
            {
                continue;
            }

            if (std::string(dir->d_name).substr(std::string(dir->d_name).size() - 4) == ".txt")
            {
                std::string tempPath = filepath.getTempEvidencePath() + "/" + std::string(dir->d_name);
                evidenceFilePathLists.push_back(tempPath);
            }
        }
        closedir(d);
    }

    d = opendir(filepath.getSplitReadPath().c_str());
    if (d)
    {
        struct dirent *dir;
        while (dir = readdir(d))
        {

            if (std::string(dir->d_name).size() < 4)
            {
                continue;
            }

            if (std::string(dir->d_name).substr(std::string(dir->d_name).size() - 4) == ".txt")
            {
                std::string tempPath = filepath.getSplitReadPath() + "/" + std::string(dir->d_name);
                evidenceFilePathLists.push_back(tempPath);
            }
        }
        closedir(d);
    }

    if (evidenceFilePathLists.size() == 0)
    {
        return;
    }

    std::sort(evidenceFilePathLists.begin(), evidenceFilePathLists.end());

    //write file
    std::string writeFinal = filepath.getAllEvidencePath();
    int count = 0;

    ReadDepthAnalysis rda(&filepath);

    for (auto n : evidenceFilePathLists)
    {
        std::vector<std::string> cache;

        std::string line;
        std::ifstream myfile(n);
        if (myfile.is_open())
        {
            std::string svtype = n.substr(n.size() - 7, 3);

            // if (svtype != "INS")
            // {
            //     continue;
            // }

            while (getline(myfile, line))
            {
                Evidence e;
                e.setEvidenceByString(line);
                if (rda.analyzeByEvidence(e))
                {
                    cache.push_back(line);
                    count++;
                }
            }
            myfile.close();

            std::ofstream writefile;
            writefile.open(writeFinal, std::ios_base::app);
            for (auto a : cache)
            {
                writefile << a << std::endl;
            }
            writefile.close();

            // std::cout << "=================" << std::endl;
        }

        cache.clear();
    }

}

void Caller::mergeReadDepthFile()
{
    std::vector<std::string> evidenceFilePathLists;

    DIR *d;
    struct dirent *dir;
    d = opendir(filepath.getTempEvidencePath().c_str());
    if (d)
    {
        while (dir = readdir(d))
        {

            if (std::string(dir->d_name).size() < 4)
            {
                continue;
            }

            if (std::string(dir->d_name).substr(std::string(dir->d_name).size() - 4) == ".txt")
            {
                std::string tempPath = filepath.getTempEvidencePath() + "/" + std::string(dir->d_name);
                evidenceFilePathLists.push_back(tempPath);
            }
        }
        closedir(d);
    }

    std::sort(evidenceFilePathLists.begin(), evidenceFilePathLists.end());

    //write file
    std::string writeFinal = filepath.getAllEvidencePath();
    int count = 0;
    ReadDepthAnalysis rda(&filepath);
    for (auto n : evidenceFilePathLists)
    {
        std::vector<std::string> cache;

        std::string line;
        std::ifstream myfile(n);
        if (myfile.is_open())
        {

            std::string svtype = n.substr(n.size() - 7, 3);
            // if (svtype != "INV")
            // {
            //     continue;
            // }
            std::string chr = n.substr(filepath.getTempEvidencePath().size() + 1, n.size() - filepath.getTempEvidencePath().size() - 9);
            std::cout << chr << std::endl;

            while (getline(myfile, line))
            {
                Evidence e;
                e.setEvidenceByString(line);
                if (rda.analyzeByEvidence(e))
                {
                    cache.push_back(line);
                    count++;
                }
            }
            myfile.close();

            std::ofstream writefile;
            writefile.open(writeFinal, std::ios_base::app);
            for (auto a : cache)
            {
                writefile << a << std::endl;
            }
            writefile.close();

            std::cout << "=================" << std::endl;
        }

        cache.clear();
    }

    std::cout << count << std::endl;
}

void Caller::catfile()
{
    std::cout << std::endl;
    std::cout << "------------------------------" << std::endl;
    std::cout << "# Merge and filter evidence :" << std::endl;
    std::cout << "------------------------------" << std::endl;

    // std::cout << ""
    catEvidenceFile();
}

int Caller::writeFile(Evidence vr)
{
    std::ofstream myfile;
    myfile.open(filepath.getOutputPath() + "/analysis/variant/" + vr.getChr() + "." + vr.getVariantType() + ".vcf", std::ios_base::app);
    // if (vr.getVariantType()=="DEL") {
    // myfile.open(filepath.getEvidencePath() + "/" + vr.getChromosome() + "." + vr.getVariantType() + ".vcf", std::ios_base::app);
    // std::cout << filepath.getOutputPath() + "/analysis/variant/" + vr.getChr() + "." + vr.getVariantType() + ".vcf" << std::endl;
    // }else {
    //     return 0;
    // }
    vr.setID("BOLT" + std::to_string(vcfIdNumber));
    myfile << vr.getResultVcfFormatString() << std::endl;
    //  std::cout << "mapq : " << vr.getMapQVector()->size() << vr.convertMapQlistToCommaString() << std::endl;
    vcfIdNumber++;
    myfile.close();
    return 0;
}

int Caller::findBreakPoint()
{
    std::cout << std::endl;
    std::cout << "------------------------------" << std::endl;
    std::cout << "# Refine breakpoint : " << std::endl;
    std::cout << "------------------------------" << std::endl;

    // removeResult();

    //    ReadDepthHelper readdepthHelper;
    //    readdepthHelper.loadReadDepthFile()

    FastaReader fastaReader;
    fastaReader.setFilePath((filepath.getReferencePath()));
    fastaReader.setIndexFilePath((filepath.getReferencePath()) + ".fai");
    fastaReader.initialize();

    samFile *inFile = sam_open(filepath.getSamplePath().c_str(), "r");
    if (inFile == NULL)
    {
        return 1;
    }

    hts_idx_t *bam_index = sam_index_load(inFile, filepath.getSamplePath().c_str());
    if (bam_index == NULL)
    {
        return 1;
    }

    std::mutex mxRead;
    std::mutex mxWriteFile;
    EvidenceProvider ep(&filepath);
    int sizeLoop = ep.getEvidenceSize();
    // std::cout << "sizeLoop : " << sizeLoop << std::endl;

    int countRunEvidence = 0;
    // tbb::task_scheduler_init init(1);

    ReadDepthAnalysis rda(&filepath);
    // float progress = 0.0;
    // float incrementevery = float(1)/float(sizeLoop);
    // std::cout << incrementevery << std::endl;

    tbb::parallel_for(0, sizeLoop, [&](int i) {
        

        mxRead.lock();
        // int barWidth = 50;
        // std::cout << "[";
        // int pos = barWidth * progress;
        // for (int i = 0; i < barWidth; ++i) {
        //     if (i < pos) std::cout << "=";
        //     else if (i == pos) std::cout << ">";
        //     else std::cout << " ";
        // }
        // std::cout << "] " << int(progress * 100.0) << " %\r";
        // std::cout.flush();

        // progress += incrementevery;

        if (ep.isEmpty())
        {
            std::cout << "end" << std::endl;
        }
        Evidence thisEvidence = ep.getEvidence();
        Evidence variantresult;
        // std::cout << thisEvidence.getPos() << " " << thisEvidence.getChr()
        // << " / " << thisEvidence.getEnd() << " " << thisEvidence.getEndChr()
        // << std::endl;

        mxRead.unlock();

        // if (thisEvidence.getVariantType() != "INV")
        // {
        //     goto skip;
        // }
        
        //    std::cout << thisEvidence.getPos() << " / " << thisEvidence.getEnd() << std::endl;
        if (thisEvidence.getVariantType() == "DEL")
        {
            // EvidenceFilter ef;
            // if (ef.passFilterEvidence(&thisEvidence))
            // {
            RefiningDeletion rfd;
            rfd.setHtsIndex(bam_index);
            rfd.setFilePath(&filepath);
            rfd.setEvidence(thisEvidence);
            rfd.setSampleStat(&samplestat);
            rfd.setFastaReader(fastaReader);
            rfd.execute();
            variantresult = rfd.getVariantResult();
            // std::cout << variantresult.getResultVcfFormatString() << std::endl;
            // }
        }
        else if (thisEvidence.getVariantType() == "DUP")
        {
            RefiningTandemDuplication rfd;
            rfd.setHtsIndex(bam_index);
            rfd.setFilePath(&filepath);
            rfd.setEvidence(thisEvidence);
            rfd.setSampleStat(&samplestat);
            rfd.setFastaReader(fastaReader);
            rfd.execute();
            variantresult = rfd.getVariantResult();
            // std::cout << variantresult.getResultVcfFormatString() << std::endl;
        }
        else if (thisEvidence.getVariantType() == "INS")
        {

            RefiningInsertion rfd;
            rfd.setHtsIndex(bam_index);
            rfd.setFilePath(&filepath);
            rfd.setEvidence(thisEvidence);
            rfd.setSampleStat(&samplestat);
            rfd.setFastaReader(fastaReader);
            rfd.execute();
            variantresult = rfd.getVariantResult();
            // std::cout << variantresult.getResultVcfFormatString() << std::endl;
        }
        else if (thisEvidence.getVariantType() == "INV")
        {
            RefiningInversion rfd;
            rfd.setHtsIndex(bam_index);
            rfd.setFilePath(&filepath);
            rfd.setEvidence(thisEvidence);
            rfd.setSampleStat(&samplestat);
            rfd.setFastaReader(fastaReader);
            rfd.execute();
            variantresult = rfd.getVariantResult();
            // std::cout << variantresult.getResultVcfFormatString() << std::endl;
        }
        else if (thisEvidence.getVariantType() == "BND")
        {
            RefiningTranslocation rfd;
            rfd.setHtsIndex(bam_index);
            rfd.setFilePath(&filepath);
            rfd.setEvidence(thisEvidence);
            rfd.setSampleStat(&samplestat);
            rfd.setFastaReader(fastaReader);
            rfd.execute();
            variantresult = rfd.getVariantResult();
        }

        // skip:

        VariantResultFilter vrf;
        if (vrf.passFilterSV(&variantresult))
        {
            std::cout << variantresult.getResultVcfFormatString() << std::endl;

            mxWriteFile.lock();
            // rda.analyzeByBreakPoint(variantresult);
            // if (variantresult.isQuailtyPass()) {
            writeFile(variantresult);
            // }
            mxWriteFile.unlock();
        }
        countRunEvidence++;
    });
}

void Caller::removeResult()
{
    auto thispathfile = filepath.getOutputPath() + "/result.vcf";
    if (remove(thispathfile.c_str()) != 0)
        perror("Error deleting file");
    else
        puts("File successfully deleted");
}

void Caller::debugEvidenceProvider()
{
    int count = 0;
    int loop = 110;
    tbb::parallel_for(0, loop, [&](int i) {
        count++;
        if (count == 90)
        {
            loop += 10;
        }
        std::cout << "Hello" << i << "/" << count << std::endl;
    });

    return;
}
