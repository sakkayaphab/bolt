#include "cli.h"
#include <iostream>
#include <string>
#include <thread>
#include <cstdlib>
#include <fstream>
#include "caller/caller.h"
#include "caller/evidence.h"
#include "caller/readparser.h"
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

void Cli::ShowHelp() const
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

void Cli::showHelpCallSV() const
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

std::string Cli::getCommand() const
{
    if (args_lists.size() > 0)
    {
        return args_lists.at(0);
    }
    return "";
}

std::string Cli::getArgumentValue(const std::string& flag) const
{
    for (size_t i = 0; i < args_lists.size(); ++i)
    {
        if (args_lists[i] == flag && i + 1 < args_lists.size())
        {
            return args_lists[i + 1];
        }
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

    // Parse arguments using helper function
    std::string bamPath = getArgumentValue("-b");
    std::string refPath = getArgumentValue("-r");
    std::string outPath = getArgumentValue("-o");
    std::string threadStr = getArgumentValue("-t");

    // Validate required arguments
    if (bamPath.empty())
    {
        std::cout << "Error: BAM file path is required (-b)" << std::endl;
        return EXIT_FAILURE;
    }

    if (refPath.empty())
    {
        std::cout << "Error: Reference file path is required (-r)" << std::endl;
        return EXIT_FAILURE;
    }

    if (outPath.empty())
    {
        std::cout << "Error: Output path is required (-o)" << std::endl;
        return EXIT_FAILURE;
    }

    // Basic file existence checks
    {
        std::ifstream bamFile(bamPath);
        if (!bamFile.good())
        {
            std::cout << "Error: Cannot access BAM file: " << bamPath << std::endl;
            return EXIT_FAILURE;
        }
    }
    
    {
        std::ifstream refFile(refPath);
        if (!refFile.good())
        {
            std::cout << "Error: Cannot access reference file: " << refPath << std::endl;
            return EXIT_FAILURE;
        }
    }


    // Parse thread count
    unsigned int threads = std::thread::hardware_concurrency();
    if (!threadStr.empty())
    {
        try
        {
            int threadCount = std::stoi(threadStr);
            if (threadCount <= 0)
            {
                std::cout << "Error: Thread count must be positive" << std::endl;
                return EXIT_FAILURE;
            }
            threads = static_cast<unsigned int>(threadCount);
            std::cout << "Using " << threads << " threads" << std::endl;
        }
        catch (const std::invalid_argument& e)
        {
            std::cout << "Error: Invalid thread count format" << std::endl;
            return EXIT_FAILURE;
        }
        catch (const std::out_of_range& e)
        {
            std::cout << "Error: Thread count out of range" << std::endl;
            return EXIT_FAILURE;
        }
    }


    Caller caller(bamPath, refPath, outPath);
    caller.showinfo();
    caller.setParallel(threads);
    caller.execute();
    caller.catfile();
    caller.findBreakPoint();
    caller.refineDelpthBlock();

    return EXIT_SUCCESS;
}

int Cli::debug() const
{
    std::cout << "Hello world" << std::endl;

    return EXIT_SUCCESS;
}