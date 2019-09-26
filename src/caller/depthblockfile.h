#ifndef DEPTHBLOCKFILE_H
#define DEPTHBLOCKFILE_H
#include <string>
#include "filemanager.h"
#include "readdepthhelper.h"

class DepthBlockFile
{
private:
  std::string filepath;
  FileManager *filemanager;
  // int32_t avgReadDepthFocus=0;
  std::map<int32_t, ReadDepthHelper::ReadDepthVector> mapReadDepthLineSegment;

public:
  std::string getFilePath();
  DepthBlockFile();
  void setFileManager(FileManager *filemanager);
  void execute();
  void loadDataToCache(std::string filepath);
  ReadDepthHelper::ReadDepthVector getBlock(int32_t number);
  std::vector<std::string> split(const std::string &s, char delimiter);
  ReadDepthHelper::ReadDepthVector findBlockWithFile(int32_t number,std::string filepath,int32_t scope);
  void addToMapReadDepthLineSegment(std::vector<std::string> *lineBuffer);
};

#endif