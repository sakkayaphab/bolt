#ifndef CORE_H
#define CORE_H
#include <iostream>
#include <stdio.h>
#include <stdint.h>
#include <vector>
#include <map>
#include "bwtscore.h"
#include <algorithm>

class Core
{
private:
  struct BwtStore
  {
    char first;
    uint32_t firstindex;
    char last;
    uint32_t lastindex;
    uint32_t lasttonext;
    uint32_t nexttolast;
    uint32_t positionoffirst;
  };

  struct RotationString
  {
    std::string text;
    uint32_t position;

    bool operator<(const RotationString &rhs) const
    {
      return (text < rhs.text);
    }
  };

  std::vector<BwtStore> bwtstore;

  uint32_t startPosition;

  std::map<char, std::vector<int>> mapRankFirst;
  std::map<char, std::vector<int>> mapRankLast;

  void setMapFirst();
  void setMapLast();
  std::vector<BwtScore> runNext(std::string *str_m, int startindex, int32_t allowMissMatch, int32_t allowSkip, int32_t minimumMatch);
  std::vector<BwtScore> runBackAllowOnlyMissmatch(std::string *str_m, int startindex, int32_t allowMissMatch,int32_t allowSkip, int32_t minimumMatch);
  bool compareString(const RotationString &a, const RotationString &b);

public:
  
  void printBwtStore();
  Core();
  void build(const char *seq);
  void setStartPosition(uint32_t position);

  std::vector<int> getIndexFromMapFirst(std::string *str_m);
  std::vector<int> getIndexFromMapLast(std::string *str_m);
  void bwtviaBwm(const char *seq);
  void rotataions(const char *seq);
  // F2L
  // main = hello
  // find = hel
  // xxx--
  std::vector<BwtScore> findStartToEnd(std::string *str_m, int allowMissMatch, int allowSkip, int minimumMatch);

  //L2F
  // main = hello
  std::vector<BwtScore> findEndToStart(std::string *str_m, int allowMissMatch, int allowSkip,int minimumMatch);

};

#endif