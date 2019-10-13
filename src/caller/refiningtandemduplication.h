#ifndef REFININGTANDEMDUPLICATION_H
#define REFININGTANDEMDUPLICATION_H
#include "samplestat.h"
#include <stdint.h>
#include <string>
#include <map>
#include <vector>
#include "refiningsv.h"

class RefiningTandemDuplication : public RefiningSV
{
private:
  void first();
  void second();
  void refineStartToEnd(const char *range);
  void refineEndToStart(const char *range);

  Evidence resultFirst;
  Evidence resultSecond;

public:
  RefiningTandemDuplication();
  void execute();
  Evidence calculateFinalBreakpoint(std::map<std::pair<int32_t, int32_t>, RefiningSV::MatchRead> *listPosition);
  Evidence getBestResult(Evidence r1, Evidence r2);
};

#endif