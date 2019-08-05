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

public:
  RefiningTandemDuplication();
  void execute();
};

#endif