#ifndef REFININGTRANSLOCATION_H
#define REFININGTRANSLOCATION_H
#include "samplestat.h"
#include <stdint.h>
#include <string>
#include <map>
#include <vector>
#include "refinesv.h"

class RefineTranslocation : public RefineSV
{
private:
  void first();
  void refineStartToEnd(const char *range);

public:
  RefineTranslocation();
  void execute();
  void calculateFinalBreakpoint(std::map<std::pair<int32_t, int32_t>, RefineSV::MatchRead> *listPosition);
};

#endif