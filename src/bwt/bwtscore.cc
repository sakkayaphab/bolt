#include "bwtscore.h"

BwtScore::BwtScore()
{
}

void BwtScore::setPos(int32_t pos)
{
    position = pos;
}

void BwtScore::setHit(int hit_m)
{
    hit = hit_m;
}

void BwtScore::setMiss(int miss_m)
{
    miss = miss_m;
}

void BwtScore::setLastMatch(bool value)
{
    lastmatch = value;
}

bool BwtScore::isLastMatch()
{
    lastmatch = false;
}

void BwtScore::setCurrentIndex(int32_t currentIndex_m)
{
    currentIndex = currentIndex_m;
}

int32_t BwtScore::getPos()
{
    return position;
}

void BwtScore::setEnd(int32_t m_end)
{
    end = m_end;
}

int32_t BwtScore::getEnd()
{
    return end;
}

int BwtScore::getHit()
{
    return hit;
}

int BwtScore::getMiss()
{
    return miss;
}

int32_t BwtScore::getCurrentIndex()
{
    return currentIndex;
}

void BwtScore::increaseOneOnHit()
{
    hit++;
}

void BwtScore::decreaseOneOnHit()
{
    hit--;
}

void BwtScore::decreaseOnHit(int number)
{
    hit -= number;
}

void BwtScore::setMatchSeq(std::string seq)
{
    seqMatch = seq;
}

std::string BwtScore::getMatchSeq()
{
    return seqMatch;
}