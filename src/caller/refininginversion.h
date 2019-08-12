#ifndef REFININGINVERSION_H
#define REFININGINVERSION_H
#include "samplestat.h"
#include <stdint.h>
#include <string>
#include <map>
#include <vector>
#include "refiningsv.h"

class RefiningInversion : public RefiningSV
{
private:
  int minimum = 2;
  int minimumMaxMatchSeq = 10;
    void first();
  void second();
  void refineStartToEnd(const char *range);
  void refineEndToStart(const char *range);

public:
  RefiningInversion();
     void execute();
     void calculateFinalBreakpoint(std::map<std::pair<int32_t, int32_t>, RefiningSV::MatchRead> *listPosition);
};

#endif