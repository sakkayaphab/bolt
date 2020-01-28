#ifndef SPLITREAD_H
#define SPLITREAD_H
#include <string>
#include "readparser.h"
#include "samplestat.h"
#include "refinesv.h"
#include "evidence.h"
#include "filemanager.h"

class SplitRead
{
private:
  ReadParser *readparser;
  SampleStat *samplestat;
  std::string chrname;
  std::vector<ReadParser::SATag> satag;
  FileManager *filepath;

  std::map<std::pair<int32_t, int32_t>, RefineSV::MatchRead> mapDUP;
  std::map<std::pair<int32_t, int32_t>, RefineSV::MatchRead> mapDEL;
   std::map<std::pair<int32_t, int32_t>, RefineSV::MatchRead> mapSmallDEL;
  std::map<std::pair<int32_t, int32_t>, RefineSV::MatchRead> mapINV;
  std::map<std::pair<int32_t, int32_t>, RefineSV::MatchRead> mapSmallINS;
  std::map<std::pair<int32_t, int32_t>, RefineSV::MatchRead> mapINS;
    std::map<std::pair<int32_t, int32_t>, RefineSV::MatchRead> mapTRA;

  std::vector<Evidence> vecINV;

  int vcfIdNumber = 0;

public:
  SplitRead(std::string chrname,ReadParser *readparser, SampleStat *samplestate,FileManager *filepath);
  void findInsertionInRead();
  void updateRead();
  void findTandemDuplication();
  void findDeletion();
  void findDeletionInRead();
  void findInversion();
  void printResult();
  void removeDuplicateResult(std::vector<Evidence> *vec);
  bool checkBetween(int32_t pos, int32_t targetPos, int32_t overlapped);
  int writeFile(Evidence vr);
  void printDeletion();
  void printSmallDeletion();
  void printSmallInsertion();
  void findInsertion();
  void printDuplication();
  void printInversion();
  void printInsertion();
  std::vector<Evidence> convertMapToEvidenceList(std::map<std::pair<int32_t, int32_t>, RefineSV::MatchRead> *mapSV, std::string svtype, std::string mark);
  void mergeEvidence(std::vector<Evidence> *elist);
  void setAllCIEvidence(std::vector<Evidence> *elist,int32_t rangePos);
  void filterEvidenceList(std::vector<Evidence> *elist);
  void filterLengthMinEvidenceList(std::vector<Evidence> *elist,int32_t min);
  void filterLengthMaxEvidenceList(std::vector<Evidence> *elist,int32_t min);
    void filterMapQLowerThan(uint8_t mapq,std::vector<Evidence> *elist);
     void filterFrequencyLowerThan(int number,std::vector<Evidence> *elist);
     bool haveSmallDeletion();
    int getDivider(int value, int top, int down, int minimum);
    void findTranslocation();
    void printTranslocation();
    std::vector<Evidence> convertMapToEvidenceListTRA(std::map<std::pair<int32_t, int32_t>, RefineSV::MatchRead> *mapSV, std::string svtype, std::string mark);

  
};

#endif