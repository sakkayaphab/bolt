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
};

#endif