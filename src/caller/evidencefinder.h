#ifndef EVIDENCEFINDER_H
#define EVIDENCEFINDER_H
#include "samplestat.h"
#include "filemanager.h"
#include <stdint.h>
#include <string>
#include <map>
#include <vector>
#include "evidence.h"
#include "readparser.h"

class EvidenceFinder
{
  public:
  struct ReadDepthDetail
  {
    int RD = 0;
    int DEL1 = 0;
    int INS1 = 0;
    int INV1 = 0;
    int DUP1 = 0;
    int TRA1 = 0;
    int DEL2 = 0;
    int INS2 = 0;
    int INV2 = 0;
    int DUP2 = 0;
    int TRA2 = 0;
    int SCF = 0;
    int SCL = 0;
    int R1_MUN = 0;
    int R2_MUN = 0;
  };

private:
  hts_idx_t *bam_index = NULL;
  SampleStat *samplestat;
  FileManager *filepath;
  std::string *target_chromosome;
  std::map<int32_t, ReadDepthDetail> ReadDepthLineSegment;
  bam_hdr_t *bam_header;
  bam1_t *read;
  ReadParser readparser;
  void updateReadDepthSV(ReadDepthDetail *rdd);
  void checkNormalRead(ReadDepthDetail *rdd);
  int32_t insertSizeFirstRead=0;
  int32_t insertSizeSecondRead=0;

  void setupReadParser(bam1_t *alnT);

  void findEvidence();

  uint32_t currentPos = 0;
  uint32_t currentMPos = 0;

public:

  EvidenceFinder(SampleStat *samplepath_T, FileManager *filepath, std::string *target_chromosome);
  void setSampleStat(SampleStat *samplepath_T);
  void setBamHeader(bam_hdr_t *bam_header);
  void setFilePath(FileManager *filepath);
  void setTargetChromosome(std::string *target_chromosome);
  void execute();
  void setHtsIndex(hts_idx_t *index);
  bool isDEL();
  bool isDUP();
  bool isINV();
  bool isINS();
  bool isTRA();


  bool isDeletion();
  bool isTandemDuplication();
  bool isInversion();
  bool isInsertion();
  bool isTranslocation();
  bool isDefaultOrientation();

};

#endif