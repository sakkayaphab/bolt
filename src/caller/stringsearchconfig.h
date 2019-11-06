#ifndef STRINGSEARCHCONFIG_H
#define STRINGSEARCHCONFIG_H
#include <string>

class StringSearchConfig
{
private:
    int32_t allowmissmatch = 0;
    int32_t maxAllowAlign = 0;
    int32_t maxContinueMissMatch = 0;
public:
    StringSearchConfig();
    void setAllowMissMatch(int32_t number);
    int32_t getAllowMissMatch();

    void setMaxAllowAlign(int32_t number);
    int32_t getMaxAllowAlign();

    void setMaxContinueMissMatch(int32_t number);
    int32_t getMaxContinueMissMatch();
};

#endif