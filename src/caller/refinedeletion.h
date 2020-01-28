#ifndef REFININGDELETION_H
#define REFININGDELETION_H
#include "samplestat.h"
#include <stdint.h>
#include <string>
#include <map>
#include <vector>
#include "refinesv.h"
#include "readparser.h"
#include <iostream>
#include <queue>

class RefineDeletion : public RefineSV
{
private:

  void first();
  void second();
  void refineStartToEnd(const char *range);
  void refineEndToStart(const char *range);

  Evidence resultFirst;
  Evidence resultSecond;

public:
  RefineDeletion();
  ~RefineDeletion();
  void execute();
  void approximate();
  Evidence calculateFinalBreakpoint(std::map<std::pair<int32_t, int32_t>, RefineSV::MatchRead> *listPosition);
  int getNumberMapQ(std::vector<uint8_t> mapqlist,uint8_t start,uint8_t end);
  Evidence getBestResult(Evidence r1, Evidence r2);
};

#endif