#ifndef REFININGTANDEMDUPLICATION_H
#define REFININGTANDEMDUPLICATION_H
#include "samplestat.h"
#include <stdint.h>
#include <string>
#include <map>
#include <vector>
#include "refinesv.h"

class RefineTandemDuplication : public RefineSV
{
private:
  void first();
  void second();
  void refineStartToEnd(const char *range);
  void refineEndToStart(const char *range);

  Evidence resultFirst;
  Evidence resultSecond;

public:
  RefineTandemDuplication();
  void execute();
  Evidence calculateFinalBreakpoint(std::map<std::pair<int32_t, int32_t>, RefineSV::MatchRead> *listPosition);
  Evidence getBestResult(Evidence r1, Evidence r2);
};

#endif