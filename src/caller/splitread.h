#ifndef SPLITREAD_H
#define SPLITREAD_H
#include <string>
#include "readparser.h"
#include "samplestat.h"
#include "refiningsv.h"

class SplitRead
{
private:
  ReadParser *readparser;
  SampleStat *samplestate;
  std::vector<ReadParser::SATag> satag;

  std::map<std::pair<int32_t, int32_t>, RefiningSV::MatchRead> mapDUP;

public:
  SplitRead(ReadParser *readparser,SampleStat *samplestate);
  void updateRead();
  void findTandemDuplication();
  void printResult();
};

#endif