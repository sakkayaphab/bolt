#ifndef REFINEDEPTHBLOCK_H
#define REFINEDEPTHBLOCK_H
#include "evidence.h"
#include "filemanager.h"
#include "depthblockfile.h"

class RefineDepthBlock
{
private:
  FileManager *filemanager;
  int32_t vcfIdNumber = 0;
  DepthBlockFile rdf;
  int32_t roundConfig = 250;
public:
  RefineDepthBlock();
  void execute();
  std::vector<std::string> getPathVCFFiles();
  void setFileManager(FileManager *filemanager);
  std::vector<Evidence> getEvidenceByFilepath(std::string filepaht);
  std::vector<Evidence> getResultWithOutOverlapped(std::vector<Evidence> *master, std::vector<Evidence> *slave);
  std::vector<Evidence> getRefineResultDeletion(std::vector<Evidence> *master);
  std::vector<Evidence> getRefineResultDuplication(std::vector<Evidence> *master);
  std::vector<Evidence> getRefineResultInversion(std::vector<Evidence> *master);
  void writeFile(std::vector<Evidence> *master);
  int32_t roundNumber(int32_t number,int32_t round);
  int32_t nextNumber(int32_t number,int32_t round);
  int32_t previousNumber(int32_t number,int32_t round);
};

#endif