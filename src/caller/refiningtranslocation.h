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
};

#endif