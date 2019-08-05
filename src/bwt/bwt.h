#ifndef BWT_H
#define BWT_H

#include "bwtscore.h"
#include "core.h"

class Bwt {

private:
    Core core;
public:
    Bwt();
    void build(const char *seq);
    void setStartPosition(int32_t position);
    std::vector<BwtScore> autoFindStartToEnd(std::string *str_m,int minimumMatch);
    std::vector<BwtScore> autoFindEndToStart(std::string *str_m,int minimumMatch);

    BwtScore getResultMaxHit(std::vector<BwtScore> *bwtscorelist);
};


#endif
