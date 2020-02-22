#include <iostream>
#include <stdio.h>
#include "cli.h"

int main(int argc, char **argv)
{
    
    Cli cli(argc, argv);

    if (cli.getCommand()=="call")
    {
        return cli.callSV();
    }
    else if (cli.getCommand()=="debug")
    {
        return cli.debug();
    }
    else if (cli.getCommand()=="version")
    {
        std::cout << "VERSION:" << std::endl;
        std::cout << "\t0.5.0" << std::endl;
        return EXIT_SUCCESS;
    }
    else if (cli.getCommand()=="")
    {
        cli.ShowHelp();
        return EXIT_SUCCESS;
    }
    else
    {
        std::cout << cli.getCommand() <<" : command not found" << std::endl;
        std::cout << std::endl;
        cli.ShowHelp();
        return EXIT_FAILURE;
    }

    return 0;
}