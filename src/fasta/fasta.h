#ifndef FASTA_H
#define FASTA_H
#include <string>
#include <vector>
#include <map>
#include <stdint.h>

class Fasta
{
private:
  std::string filepath;

public:
  // Fasta(std::string t_filepath, std::string t_indexfilepath);
  Fasta();
  void setFilePath(std::string filepath);
  std::string getFilePath();
  void initialize();

  std::string  getSeqByOffset(uint64_t start,uint64_t end);
};

#endif 