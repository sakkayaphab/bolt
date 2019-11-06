#ifndef VCF_H
#define VCF_H
#include <string>
#include <vector>
#include "dataline.h"
#include "metainfoline.h"

class Vcf
{
private:


protected:
  std::vector<Dataline> datalines;
  MetaInfoline metainfolines;

public:
  Vcf();
  void addDataline(Dataline dataline);
  Dataline *getDatalineAt(uint32_t number);
  MetaInfoline *getMetaInfoline();
  void setMetaInfoLines(MetaInfoline meta);
};

#endif