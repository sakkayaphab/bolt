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
  std::map<std::pair<int32_t, int32_t>, RefiningSV::MatchRead> mapINV;

  std::vector<Evidence> vecINV;

  int vcfIdNumber = 0;

public:
  SplitRead(std::string chrname,ReadParser *readparser, SampleStat *samplestate,FileManager *filepath);
  void updateRead();
  void findTandemDuplication();
  void findDeletion();
  void findInversion();
  void printResult();
  void removeDuplicateResult(std::vector<Evidence> *vec);
  bool checkBetween(int32_t pos, int32_t targetPos, int32_t overlapped);
  int writeFile(Evidence vr);
  void printDeletion();
  void printDuplication();
  void printInversion();
  std::vector<Evidence> convertMapToEvidenceList(std::map<std::pair<int32_t, int32_t>, RefiningSV::MatchRead> *mapSV, std::string svtype);
  void mergeEvidence(std::vector<Evidence> *elist);
  void setAllCIEvidence(std::vector<Evidence> *elist);
  void filterEvidenceList(std::vector<Evidence> *elist);
  void filterLengthMinEvidenceList(std::vector<Evidence> *elist,int32_t min);
  void filterLengthMaxEvidenceList(std::vector<Evidence> *elist,int32_t min);

  
};

#endif