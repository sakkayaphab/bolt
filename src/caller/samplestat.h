#ifndef SAMPLESTAT_H
#define SAMPLESTAT_H
#include <htslib/sam.h>
#include <string>

class SampleStat
{
private:
  std::string sample_path;

  int insertsize_median=0;
  int insertsize_sd=0;
  int avg_rd=0;
  int read_length=0;
  
  //config
  int countMax;


  void findReadLength();
  void findMedianSampleStat();
  void findSDSampleStat();

public:
  SampleStat(std::string samplepath);
  SampleStat();
  void execute();
  void setSamplePath(std::string t_sample_path);
  void setNumberOfRead(int number);
  int32_t getReadLength();
  int getMedianSampleStat();
  int getSDSampleStat();
};

#endif