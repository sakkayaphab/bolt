#ifndef SPLITREAD_H
#define SPLITREAD_H
#include <string>
#include "readparser.h"
#include "samplestat.h"
#include "refiningsv.h"
#include "evidence.h"
#include "filemanager.h"

class SplitRead
{
private:
  ReadParser *readparser;
  SampleStat *samplestate;
  std::string chrname;
  std::vector<ReadParser::SATag> satag;
  FileManager *filepath;

  std::map<std::pair<int32_t, int32_t>, RefiningSV::MatchRead> mapDUP;
  std::map<std::pair<int32_t, int32_t>, RefiningSV::MatchRead> mapDEL;

  std::vector<Evidence> vecDEL;

  int vcfIdNumber = 0;

public:
  SplitRead(std::string chrname,ReadParser *readparser, SampleStat *samplestate,FileManager *filepath);
  void updateRead();
  void findTandemDuplication();
  void findDeletion();
  void printResult();
  void removeDuplicateResult(std::vector<Evidence> *vec);
  bool checkBetween(int32_t pos, int32_t targetPos, int32_t overlapped);
  int writeFile(Evidence vr);
};

#endif