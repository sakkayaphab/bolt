#ifndef SAMPLESTAT_H
#define SAMPLESTAT_H
#include <htslib/sam.h>
#include <string>
#include <vector>

class SampleStat
{
private:
  std::string sample_path;

  int insertsize_median=0;
  int insertsize_sd=0;
  int avg_rd=0;
  int read_length=0;
  
  //config
  int countMax=100000;


  void findReadLength();
  void findMedianSampleStat(std::vector<int32_t> *insertlist);
  void findSDSampleStat(std::vector<int32_t> *insertlist);
  std::vector<int32_t> getInsertSizeList(int64_t numberofread);

public:
  SampleStat(std::string samplepath);
  SampleStat();
  void execute();
  void setSamplePath(std::string t_sample_path);
  void setNumberOfRead(int number);
  int32_t getReadLength();
  int32_t getAverageSampleStat();
  int32_t getSDSampleStat();
};

#endif