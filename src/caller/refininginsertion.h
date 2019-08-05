#ifndef REFININGINSERTION_H
#define REFININGINSERTION_H
#include "samplestat.h"
#include <stdint.h>
#include <string>
#include <map>
#include <vector>
#include "refiningsv.h"

class RefiningInsertion : public RefiningSV
{
private:
  void first();
  void refineStartToEnd(const char *range);
  int32_t getPosMaxHitValue(std::map<int32_t, int> *map,std::map<int32_t, uint8_t> *mapMapQ);
  int getHitByPos(std::map<int32_t, int> *map, int32_t pos);
  bool isBetWeen(int32_t primary, int32_t secondary, int32_t range);
  void refineVariant(const char *range);

  std::map<int32_t, int> mapSCFirst;
  std::map<int32_t, uint8_t> mapMapQFirst;
  std::map<int32_t, int> mapSCLast;
  std::map<int32_t, uint8_t> mapMapQLast;

public:
  RefiningInsertion();
  void execute();
};

#endif