#include "vcf.h"

Vcf::Vcf()
{

}

void Vcf::addDataline(Dataline t_dataline)
{
    datalines.push_back(t_dataline);
}

Dataline *Vcf::getDatalineAt(uint32_t number)
{
    return &datalines.at(number);
}

void Vcf::setMetaInfoLines(MetaInfoline meta)
{
    metainfolines = meta;
}

MetaInfoline *Vcf::getMetaInfoline()
{
    return &metainfolines;
}