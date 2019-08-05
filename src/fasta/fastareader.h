#ifndef FASTAREADER_H
#define FASTAREADER_H
#include <string>
#include <vector>
#include <map>
#include <stdint.h>
#include "faidx.h"
#include "fasta.h"
#include <algorithm>

class FastaReader
{
private:
  Faidx faidx;
  Fasta fasta;

public:
  // Fasta(std::string t_filepath, std::string t_indexfilepath);
  FastaReader();
  void setFilePath(std::string filepath);
  std::string getFilePath();
  void setIndexFilePath(std::string indexfilepath);
  std::string getIndexFilePath();

  void initialize();
  std::string getSeqbyPosition(std::string chromosome, uint64_t start, uint64_t end);

  void replaceSeqToUppercase(std::string *str)
  {
    std::transform((*str).begin(), (*str).end(), (*str).begin(), ::toupper);
  }
};

#endif