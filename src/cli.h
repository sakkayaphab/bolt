#ifndef CLI_H
#define CLI_H
#include <string>
#include <vector>

class Cli
{
  private:
    int argc;
    char **argv;

    std::string name;
    std::string description;
    std::string version;
    std::vector<std::string> args_lists;

  public:
    Cli(int argc, char **argv);
    std::string getCommand() const;
    int callSV();
    void showHelpCallSV() const;
    int debug() const;
    void ShowHelp() const;

private:
    std::string getArgumentValue(const std::string& flag) const;

};

#endif