#include "stringsearchconfig.h"
#include <iostream>

StringSearchConfig::StringSearchConfig()
{
}

void StringSearchConfig::setAllowMissMatch(int32_t number)
{
    allowmissmatch = number;
}

int32_t StringSearchConfig::getAllowMissMatch()
{
    return allowmissmatch;
}

void StringSearchConfig::setMaxAllowAlign(int32_t number)
{
    maxAllowAlign = number;
}

int32_t StringSearchConfig::getMaxAllowAlign()
{
    return maxAllowAlign;
}

void StringSearchConfig::setMaxContinueMissMatch(int32_t number)
{
    maxContinueMissMatch = number;
}

int32_t StringSearchConfig::getMaxContinueMissMatch()
{
    return maxContinueMissMatch;
}