#include "dataline.h"

Dataline::Dataline()
{
}

void Dataline::setTextCHROM(std::string value)
{
    CHROM = value;
}
void Dataline::setTextPOS(std::string value)
{
    POS = value;
}

void Dataline::setTextID(std::string value)
{
    ID = value;
}

void Dataline::setTextREF(std::string value)
{
    REF = value;
}

void Dataline::setTextALT(std::string value)
{
    ALT = value;
}

void Dataline::setTextQUAL(std::string value)
{
    QUAL = value;
}

void Dataline::setTextFILTER(std::string value)
{
    FILTER = value;
}

void Dataline::setTextINFO(std::string value)
{
    INFO = value;
}

void Dataline::setTextFORMAT(std::string value)
{
    FORMAT = value;
}

void Dataline::setTextSAMPLE(std::vector<std::string> value)
{
    SAMPLE = value;
}

std::string Dataline::getTextCHROM()
{
    return CHROM;
}
std::string Dataline::getTextPOS()
{
    return POS;
}
std::string Dataline::getTextID()
{
    return ID;
}
std::string Dataline::getTextREF()
{
    return REF;
}
std::string Dataline::getTextALT()
{
    return ALT;
}
std::string Dataline::getTextQUAL()
{
    return QUAL;
}
std::string Dataline::getTextFILTER()
{
    return FILTER;
}
std::string Dataline::getTextINFO()
{
    return INFO;
}

std::string Dataline::getTextFORMAT()
{
    return FORMAT;
}
std::string Dataline::getTextSAMPLE()
{
    return "";
}