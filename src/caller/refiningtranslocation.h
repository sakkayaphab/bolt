#ifndef REFININGTRANSLOCATION_H
#define REFININGTRANSLOCATION_H
#include "samplestat.h"
#include <stdint.h>
#include <string>
#include <map>
#include <vector>
#include "refiningsv.h"

class RefiningTranslocation : public RefiningSV
{
private:
  void first();
  void refineStartToEnd(const char *range);

public:
  RefiningTranslocation();
  void execute();
  void calculateFinalBreakpoint(std::map<std::pair<int32_t, int32_t>, RefiningSV::MatchRead> *listPosition);
};

#endif