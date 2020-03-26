#ifndef READDEPTHANALYSIS_H
#define READDEPTHANALYSIS_H

#include <string>
#include <vector>
#include "filemanager.h"
#include "readdepthhelper.h"
#include <sstream>
#include "evidence.h"
#include "samplestat.h"
#include <map>
#include "readdepthstat.h"

class ReadDepthAnalysis {
private:
    FileManager *filemanager;
    std::string cachechr;
    int avgReadDepthFocus=0;
    int avgReadDepth = 0;

    std::map<int32_t, ReadDepthHelper::ReadDepthVector> mapReadDepthLineSegment;
    // std::vector<ReadDepthHelper::ReadDepthVector> cacheReadDepthFile;
    std::vector<ReadDepthHelper::ReadDepthVector> startFocusReadDepth;
    std::vector<ReadDepthHelper::ReadDepthVector> endFocusReadDepth;
    std::vector<std::string> split(const std::string &s, char delimiter);
    SampleStat *samplestat;
    uint32_t configRound = 250;

    ReadDepthStat readDepthStat;
    int minimumdivide = 2;
    int getDivider(int value,int top,int down,int minimum);

private:
    int sumStartDEL = 0;
    int sumStartDUP = 0;
    int sumStartINV = 0;
    int sumStartTRA = 0;
    int sumStartINS = 0;
    int sumStartSCF = 0;
    int sumStartSCL = 0;
    int sumStartR1_MUN = 0;
    int sumStartR2_MUN = 0;

    int sumEndDEL = 0;
    int sumEndDUP = 0;
    int sumEndINV = 0;
    int sumEndTRA = 0;
    int sumEndINS = 0;
    int sumEndSCF = 0;
    int sumEndSCL = 0;

    int sumEndR1_MUN = 0;
    int sumEndR2_MUN = 0;


public:
    void collectNewData();
    ReadDepthAnalysis(FileManager *filemanager);
    int32_t getRound(int32_t x, int32_t max);
    std::vector<std::int32_t> getVectorRange(int32_t pos, int32_t end);
    bool analyzeByEvidence(Evidence e);
    void setMedianReadDepth(int readdepth);
    int getReadDepthAverageFocusArea(std::vector<ReadDepthHelper::ReadDepthVector> *focusReadDepth);
    void setSampleStat(SampleStat *samplestat);
    void loadDataToCache(std::string filepath);
    void loadEvidenceFromFile();
    int getAvgReadDepth();
    void setFocusReadDepth(int32_t pos, int32_t end, std::vector<ReadDepthHelper::ReadDepthVector> *focusReadDepth);

    int getNumberReadDepthVector(std::vector<ReadDepthHelper::ReadDepthVector> focus,int blocknumber);
    bool filterDeletion(Evidence e);
    bool filterInversion(Evidence e);
    bool filterInsertion(Evidence e);
    bool filterDuplication(Evidence e);
    void loadavgReadDepthFocusStat();
    bool filterTranslocation(Evidence e);


    int getSCFFocusArea(std::vector<ReadDepthHelper::ReadDepthVector> *focusReadDepth);
    int getSCLFocusArea(std::vector<ReadDepthHelper::ReadDepthVector> *focusReadDepth);
    bool isDepthMoreThan(int depth);
};

#endif