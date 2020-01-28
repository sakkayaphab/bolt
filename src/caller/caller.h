#ifndef CALLER_H
#define CALLER_H
#include <stdint.h>
#include "samplestat.h"
#include "filemanager.h"
#include "evidence.h"
#include <thread>
#include "samplestat.h"
#include <iostream>
#include <stdio.h>
#include <iomanip>
#include <string.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <cstdlib>
#include <iostream>
#include <unistd.h>
#include "task.h"
#include "evidenceprovider.h"
#include <fstream>
#include <unistd.h>
#include <mutex>
#include <iostream>
#include <dirent.h>
#include "variantresultfilter.h"
#include "readdepthhelper.h"
#include "readdepthanalysis.h"
#include "refinedepthblock.h"
#include "depthblockfile.h"
#include <fasta/fastareader.h>
#include "refinetandemduplication.h"
#include "refinetranslocation.h"
#include "refineinversion.h"
#include "refineinsertion.h"
#include "refinedeletion.h"
#include <tbb/tbb.h>

class Caller
{

private:
  SampleStat samplestat;
  FileManager filepath;
  bam_hdr_t bam_header;
  int numberofpair_stat = 1;
  int numberofparallel = 30;
  int vcfIdNumber = 0;
  samFile *inFile = NULL;
  hts_idx_t *bam_index = NULL;
  std::mutex mxWriteFile;

public:
  Caller(std::string samplepath_T, std::string referencepath_T, std::string outputpath_T);
  ~Caller();

  void execute();
  void execSampleStat();
  void showinfo();
  void deleteAll();
  void applyBamHeader();
  void setnumberofpair_stat(int n);
  bam_hdr_t *getBamHeader();
  void setParallel(int number);
  int findBreakPoint();
  void catfile();
  void catEvidenceFile();
  void mergeReadDepthFile();
  void removeResult();
  int writeFile(Evidence vr);
  void debugEvidenceProvider();
  void prepareHts();
  void refineDelpthBlock();
  void mergeSplitRead();



};

#endif