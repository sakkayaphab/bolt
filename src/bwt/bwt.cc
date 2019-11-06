#include "bwt.h"

Bwt::Bwt()
{

}

void Bwt::build(const char *seq)
{
    core.build(seq);
}

void  Bwt::setStartPosition(int32_t position)
{
    core.setStartPosition(position);
}

std::vector<BwtScore> Bwt::autoFindStartToEnd(std::string *str_m,int minimumMatch)
{
    core.printBwtStore();
    return core.findStartToEnd(str_m,2,2,minimumMatch);
}

std::vector<BwtScore> Bwt::autoFindEndToStart(std::string *str_m,int minimumMatch)
{
    return core.findEndToStart(str_m,2,2,minimumMatch);
}


BwtScore Bwt::getResultMaxHit(std::vector<BwtScore> *bwtscorelist)
{

    int maxHit = 0;
    int position = 0;
    for (int i = 0; i < (*bwtscorelist).size(); i++)
    {
        if (bwtscorelist->at(i).getHit() > maxHit)
        {
            maxHit = bwtscorelist->at(i).getHit();
            position = i;
        }
    }

    return bwtscorelist->at(position);
}