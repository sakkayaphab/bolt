#ifndef REFININGINSERTION_H
#define REFININGINSERTION_H
#include "samplestat.h"
#include <stdint.h>
#include <string>
#include <map>
#include <vector>
#include "refiningsv.h"
#include "insertionpositiondetail.h"
#include "editdistance.h"

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

  struct CountRefineSeq
  {
    std::string seq;
    int count=0;
  };
  
  
  
  void convertMapSC();
  void clearMapSC();
  std::vector<InsertionPositionDetail> convertMapSCToVector(std::map<int32_t, InsertionPositionDetail> mapSC);
  void findBreakpoint();
  bool checkBetween(int32_t pos, int32_t targetPos, int32_t overlapped);
  void filterBreakpoint();

  std::vector<CountRefineSeq> mergeString(std::vector<std::string> fragmentlist,bool fromstart);
  bool compareEditDistance(std::string s1,std::string s2,bool fromstart);
  void substringSeq(std::string *s1,std::string *s2,bool fromstart);

public:
  RefiningInsertion();
  void execute();
  
};

#endif