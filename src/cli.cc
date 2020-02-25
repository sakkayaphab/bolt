#include "cli.h"
#include <iostream>
#include "caller/caller.h"
#include "caller/evidence.h"
#include "caller/readparser.h"
#include <string>
#include "caller/editdistance.h"
#include <thread>
#include <cstdlib>
#include "graph/graph.h"

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
    std::cout << "\t-t\tnumber of threads to use" << std::endl;

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
        return EXIT_FAILURE;
    }

    if (args_lists.at(1) == "-h" || args_lists.at(1) == "-help")
    {
        showHelpCallSV();
        return EXIT_SUCCESS;
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
        return EXIT_FAILURE;
    }

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
        return EXIT_FAILURE;
    }

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
        return EXIT_FAILURE;
    }

    // Find threads
    bool foundThread = false;
    std::string outThread;
    for (auto n : args_lists)
    {
        if (foundThread)
        {
            outThread = n;
            break;
        }

        if (n == "-t")
        {
            foundThread = true;
        }
    }


    unsigned int threads = std::thread::hardware_concurrency();
    if (outThread != "")
    {
        try
        {
            std::cout << outThread << std::endl;
            threads = std::stoi(outThread);
        }
        catch (std::invalid_argument const &e)
        {
            std::cout << "Bad input: std::invalid_argument thrown" << '\n';
            return EXIT_FAILURE;
        }
        catch (std::out_of_range const &e)
        {
            std::cout << "Integer overflow: std::out_of_range thrown" << '\n';
            return EXIT_FAILURE;
        }
    }

    Caller caller(bamPath, refPath, outPath);
    caller.showinfo();
    caller.setParallel(threads);
//    caller.execute();
//    caller.catfile();
    caller.findBreakPoint();
    caller.refineDelpthBlock();

    return EXIT_SUCCESS;
}

int Cli::debug()
{
    Graph graph;
    int kmer = 3;
    graph.buildGraph("ABC",kmer);
    graph.buildGraph("BCD",kmer);
    graph.buildGraph("DEF",kmer);

//    graph.addNodeToNode("1","2");
//    graph.addNodeToNode("1","3");
//    graph.addNodeToNode("2","2");
//    graph.addNodeToNode("2","4");
//    graph.addNodeToNode("2","4");
//    graph.addNodeToNode("3","1");
//    graph.addNodeToNode("3","2");
//    graph.addNodeToNode("3","5");
//    graph.addNodeToNode("4","3");
//    graph.addNodeToNode("4","6");
//    graph.addNodeToNode("5","6");
//    graph.addNodeToNode("6","3");

//    graph.addNodeToNode("0","1");
//    graph.addNodeToNode("0","2");
//    graph.addNodeToNode("1","2");
//    graph.addNodeToNode("1","3");
//    graph.addNodeToNode("2","3");
//    graph.addNodeToNode("3","4");
//    graph.addNodeToNode("4","0");
//    graph.addNodeToNode("4","1");
//    graph.addNodeToNode("4","5");

    graph.showAllNode();

//    std::vector<std::stack<std::string>> vStackResult = graph.findDFS("AGTCATATC","GGAGACGGT");
    std::vector<std::stack<std::string>> vStackResult = graph.findDFS("AGATCAGTGAGTTCACCACAGCCGCGCAAGGCACAGGCTCTCCACGGTCCTGCAGCCTCAGTTCCTGGGAAGCTCACAAGCAGGGGGCTGGCGTGTGAACT","ACAGCCGCGCAAGGCACAGGCTCTCCACGGTCCTGCAGCCTCAGTTCCTGGGAAGCTCACAAGCAGGGGGCTGGCGTGTGAACTCCGGTCCCGCTGTGGTC");
    return EXIT_SUCCESS;
    graph.clearMaskedEdges();
    if (vStackResult.size()==0) {
        return EXIT_SUCCESS;
    }
    for (std::stack<std::string> ss:vStackResult) {
        std::vector<std::string> vResult = graph.covertStackToVector(ss);
        graph.reverseVectorString(&vResult);
        graph.showVectorString(&vResult);
    }

    return EXIT_SUCCESS;
}