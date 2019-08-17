#ifndef REFININGDELETION_H
#define REFININGDELETION_H
#include "samplestat.h"
#include <stdint.h>
#include <string>
#include <map>
#include <vector>
#include "refiningsv.h"
#include "readparser.h"
#include <iostream>
#include <queue>

class RefiningDeletion : public RefiningSV
{
private:

  void first();
  void second();
  void refineStartToEnd(const char *range);
  void refineEndToStart(const char *range);

public:
  RefiningDeletion();
  ~RefiningDeletion();
  void execute();
  void approximate();
  void calculateFinalBreakpoint(std::map<std::pair<int32_t, int32_t>, RefiningSV::MatchRead> *listPosition);
  int getNumberMapQ(std::vector<uint8_t> mapqlist,uint8_t start,uint8_t end);
};

#endif