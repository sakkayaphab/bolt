#include "cli.h"
#include <iostream>
#include "caller/caller.h"
#include "caller/evidence.h"
#include "caller/readparser.h"
#include <string>
#include "caller/editdistance.h"

Cli::Cli(int m_argc, char **m_argv)
{
    argc = m_argc;
    argv = m_argv;

    if (argc > 1)
    {
        args_lists.assign(argv + 1, argv + argc);
    }
}

void Cli::ShowHelp()
{
    std::cout << std::endl;
    std::cout << "NAME:" << std::endl;
    std::cout << "\tBolt - a bioinformatics tool" << std::endl;
    std::cout << std::endl;
    std::cout << "USAGE:" << std::endl;
    std::cout << "\tbolt command [command options] [arguments...]" << std::endl;
    std::cout << std::endl;
    std::cout << "COMMANDS:" << std::endl;
    std::cout << "\tcall\tcall variant" << std::endl;
    std::cout << "\tversion\tShows version" << std::endl;
    std::cout << std::endl;
}

void Cli::showHelpCallSV()
{
    std::cout << std::endl;
    std::cout << "USAGE:" << std::endl;
    std::cout << "\tbolt call [command options] [arguments...]" << std::endl;
    std::cout << std::endl;
    std::cout << "COMMAND OPTIONS:" << std::endl;
    std::cout << "\t-b\tsample file path (*required)" << std::endl;
    std::cout << "\t-r\treference file path (*required)" << std::endl;
    std::cout << "\t-o\toutput path (*required)" << std::endl;
    std::cout << std::endl;
}

std::string Cli::getCommand()
{
    if (args_lists.size() > 0)
    {
        return args_lists.at(0);
    }
    return "";
}

int Cli::callSV()
{

    if (args_lists.size() < 2)
    {
        showHelpCallSV();
        return 1;
    }

    if (args_lists.at(1) == "-h" || args_lists.at(1) == "-help")
    {
        showHelpCallSV();
        return 0;
    }

    // Find BAM
    bool foundBam = false;
    std::string bamPath;
    for (auto n : args_lists)
    {
        if (foundBam)
        {
            bamPath = n;
            break;
        }

        if (n == "-b")
        {
            foundBam = true;
        }
    }

    if (bamPath == "")
    {
        std::cout << "not found bam file path" << std::endl;
        return 1;
    }

    std::cout << "bam file = " << bamPath << std::endl;

    // Find Reference
    bool foundRef = false;
    std::string refPath;
    for (auto n : args_lists)
    {
        if (foundRef)
        {
            refPath = n;
            break;
        }

        if (n == "-r")
        {
            foundRef = true;
        }
    }

    if (refPath == "")
    {
        std::cout << "not found reference file path" << std::endl;
        return 1;
    }

    std::cout << "reference file = " << refPath << std::endl;

    // Find output
    bool foundOut = false;
    std::string outPath;
    for (auto n : args_lists)
    {
        if (foundOut)
        {
            outPath = n;
            break;
        }

        if (n == "-o")
        {
            foundOut = true;
        }
    }

    if (outPath == "")
    {
        std::cout << "not found output file path" << std::endl;
        return 1;
    }

    std::cout << "output file = " << outPath << std::endl;

    return 0;
}

int Cli::debug()
{

    // std::string  sample = "/data/users/wichadak/kan/wgsim/survi.20x.bam";
    //    std::string  reference = "/data/users/wichadak/kan/reference/ucsc_hg19.fa";
    //    std::string sample = "/data/users/wichadak/kan/sample/NA12878_S1.bam";
    // std::string sample = "/home/sakkayaphab/kan/sample/G2223.remdup.uniqMap.TS.bam";
    // std::string sample = "/home/sakkayaphab/kan/wgsim/hx1f4s4full_3rdfixedv2.bam";
    // std::string sample = "/home/sakkayaphab/kan/sample/HG001.GRCh38_full_plus_hs38d1_analysis_set_minus_alts.300x.bam";
    // std::string sample = "/home/sakkayaphab/kan/sample/NA12878_S1.bam";
    // std::string sample = "/home/sakkayaphab/kan/sample/G5091.bam";

    // std::string sample = "/home/sakkayaphab/kan/sample/Sim-A_30x_bwa.bam";

    std::string sample =  "/home/sakkayaphab/kan/sample/com/sample/ERR174339/ERR174339.bam";
    // std::string sample =  "/home/sakkayaphab/kan/sample/com/sample/SRR1910366/SRR1910366.bam";

    // std::string sample = "/home/sakkayaphab/kan/sample/NA12878.hiseq.wgs.bwa.raw.bam";
    // std::string sample = "/home/sakkayaphab/kan/sra/sratoolkit.2.9.6-centos_linux64/bin/SRR390728.bam";
    // std::string sample = "/home/sakkayaphab/kan/wgsim/survi.20x.bam";
    // std::string sample = "/home/sakkayaphab/kan/wgsim/survi.50x.bam";
    // std::string sample = "/home/sakkayaphab/kan/wgsim/survi.100x.bam";
    // std::string reference = "/home/sakkayaphab/kan/reference/Homo_sapiens_assembly19.fasta";
    //  std::string reference = "/home/sakkayaphab/kan/reference/ucsc_hg19.fa";
    std::string reference = "/home/sakkayaphab/kan/reference/hs37d5.fa";

    //   std::string reference = "/home/sakkayaphab/kan/reference/GRCh38_full_plus_hs38d1_analysis_set_minus_alts/GRCh38_full_plus_hs38d1_analysis_set_minus_alts.fa";
    // std::string reference = "/home/sakkayaphab/kan/reference/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/WholeGenomeFasta/genome.fa";
    // std::string reference = "/home/sakkayaphab/kan/reference/Homo_sapiens_assembly18.fasta";
    // std::string reference = "/data/users/wichadak/kan/reference/ucsc_hg19.fa";
    // std::string output = "na12878";
    std::string output = "temp";
    // std::string output = "/data/users/wichadak/kan/bolt/temp";
    // std::string output = "/data/users/wichadak/kan/bolt/temp";
    // std::string output = "temp100";
    // std::string output = "temp50";
    // std::string name = "hello";
    // std::cout << name.substr(0,name.size()) << std::endl;
    // return 0;

    Caller caller(sample, reference, output);
    caller.showinfo();
    // caller.setParallel(40);
    // caller.execute();
    // caller.catfile();
    // caller.findBreakPoint();
    caller.refineDelpthBlock();


    return 0;
}