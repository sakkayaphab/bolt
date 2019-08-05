#ifndef WRITER_H
#define WRITER_H
#include <string>
#include <vector>
#include "vcf.h"

class Writer: public Vcf
{
  private:
    std::string outputPath;
    
  public:
    Writer();
    void setOutputPath(std::string path);
    void writeBufferAndClear();
    void writeHeader();
};

#endif