#include "metainfoline.h"

MetaInfoline::MetaInfoline()
{
    
}

std::string MetaInfoline::getDateStringNow()
{
    time_t now = time(0);
    struct tm tstruct;
    char buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "%Y%m%d", &tstruct);

    return buf;
}

void MetaInfoline::setheader_fileformat(std::string value)
{
    fileformat = value;
}
void MetaInfoline::setheader_fileDate(std::string value)
{
    fileDate = value;
}
void MetaInfoline::setheader_source(std::string value)
{
    source = value;
}
void MetaInfoline::setheader_reference(std::string value)
{
    reference = value;
}

std::string MetaInfoline::getTextFileformat()
{
    std::string text = "##fileformat="+fileformat;
    return text;
}

std::string MetaInfoline::getTextFileDate()
{
   std::string text = "##fileDate="+fileDate;
    return text;
}

std::string MetaInfoline::getTextSource()
{
    std::string text = "##source="+source;
    return text;
}

std::string MetaInfoline::getTextReference()
{
    std::string text = "##reference="+reference;
    return text;
}