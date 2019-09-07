#ifndef READDEPTH_H
#define READDEPTH_H

#include <string>
#include <vector>
#include "readdepthhelper.h"

class ReadDepth {
private:
    std::string filepath;
    std::map<std::string, int> rdmap;
public:
    ReadDepth();
    void loadFile(std::string filepath);
    void execute();
    int getReadDepthByChr(std::string chr);
};

#endif