#ifndef REFININGINSERTION_H
#define REFININGINSERTION_H
#include "samplestat.h"
#include <stdint.h>
#include <string>
#include <map>
#include <vector>
#include "refiningsv.h"
#include "insertionpositiondetail.h"

class RefiningInsertion : public RefiningSV
{
private:
  void first();
  void refineStartToEnd(const char *range);

  std::map<int32_t, InsertionPositionDetail> mapSCStart;
  std::map<int32_t, InsertionPositionDetail> mapSCEnd;


  std::vector<InsertionPositionDetail> vectorSCStart;
  std::vector<InsertionPositionDetail> vectorSCEnd;

  struct BreakpointPosition
  {
    int32_t pos;
    int32_t end;
    int frequency = 0;
    int longmapstart =0;
    int longmapend = 0;
    int score = 0;
    std::vector<uint8_t> mappingqualitylist;

    bool operator<(const BreakpointPosition &rhs) const
    {
        return (frequency < rhs.frequency);
    }
  };

  std::vector<BreakpointPosition> vectorBP;
  
  
  void convertMapSC();
  void clearMapSC();
  void convertMapSCToVector(std::map<int32_t, InsertionPositionDetail> *mapSC,std::vector<InsertionPositionDetail> *vectorSC);
  void findBreakpoint();
  bool checkBetween(int32_t pos, int32_t targetPos, int32_t overlapped);
  void filterBreakpoint();

public:
  RefiningInsertion();
  void execute();
};

#endif