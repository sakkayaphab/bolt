#ifndef REFINEDEPTHBLOCK_H
#define REFINEDEPTHBLOCK_H
#include "evidence.h"
#include "filemanager.h"
#include "depthblockfile.h"
#include "readdepthstat.h"

class RefineDepthBlock
{
private:
    FileManager *filemanager;
    int32_t vcfIdNumber = 0;
    DepthBlockFile rdf;
    int32_t roundConfig = 250;
    ReadDepthStat readDepthStat;
    SampleStat *samplestat;
    int minimumdivide = 4;

public:
    RefineDepthBlock();

    void execute();
    std::vector<std::string> getPathVCFFiles();
    void setFileManager(FileManager *filemanager);
    void setSampleStat(SampleStat *samplestat);
    std::vector<Evidence> getEvidenceByFilepath(std::string filepaht);
    std::vector<Evidence> getResultWithOutOverlapped(std::vector<Evidence> *master, std::vector<Evidence> *slave);
    std::vector<Evidence> getResultRemoveOverlapped(std::vector<Evidence> *master, std::vector<Evidence> *slave);
    std::vector<Evidence> getRefineResultDeletion(std::vector<Evidence> *master);
    std::vector<Evidence> getRefineResultDuplication(std::vector<Evidence> *master);
    std::vector<Evidence> getRefineResultInversion(std::vector<Evidence> *master);
    std::vector<Evidence> getRefineResultTranslocation(std::vector<Evidence> *master);
    std::vector<Evidence> getRefineResultInsertion(std::vector<Evidence> *master);

    void writeFile(std::vector<Evidence> *master);
    int32_t roundNumber(int32_t number,int32_t round);
    int32_t nextNumber(int32_t number,int32_t round);
    int32_t previousNumber(int32_t number,int32_t round);
    int getDivider(int value,int top,int down,int minimum);

    bool checkBetween(int32_t pos, int32_t targetPos, int32_t overlapped);
};

#endif