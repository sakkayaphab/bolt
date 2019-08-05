#ifndef FAIDX_H
#define FAIDX_H
#include <string>
#include <vector>
#include <map>
#include <stdint.h>

class Faidx
{
private:
  std::string indexfilepath;

  struct FileFormat
  {
    std::string NAME;
    uint64_t LENGTH;
    uint64_t OFFSET;
    uint64_t LINEBASES;
    uint64_t LINEWIDTH;
    uint64_t QUALOFFSET;
  };

  struct OffsetFile {
    uint64_t start;
    uint64_t end;
  };

  std::vector<FileFormat> faidxlist;

  uint64_t getApproximateLineStart(std::string chromosome,uint64_t start);
  uint64_t getApproximateLineEnd(std::string chromosome,uint64_t end);

  uint64_t getApproximateOffsetStart(std::string chromosome,uint64_t start);
  uint64_t getApproximateOffsetEnd(std::string chromosome,uint64_t end);
public:
  Faidx();
  void initialize();
  void setIndexFilePath(std::string indexfilepath);
  std::string getIndexFilePath();

  bool hasChromosomeName(std::string chr);
  uint64_t getLengthbyChromosome(std::string chr);
  uint64_t getOffsetbyChromosome(std::string chr);
  uint64_t getLinebasesbyChromosome(std::string chr);
  uint64_t getLinewidthbyChromosome(std::string chr);
  uint64_t getQualoffsetbyChromosome(std::string chr);

  uint64_t getOffsetEndByPosition(std::string chromosome, uint64_t end);
  uint64_t getOffsetStartByPosition(std::string chromosome, uint64_t start);
};

#endif