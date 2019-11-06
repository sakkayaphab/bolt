#ifndef BWTSCORE_H
#define BWTSCORE_H
#include <iostream>
#include <stdio.h>
#include <stdint.h>
#include <vector>

class BwtScore
{
  private:
    
    int32_t position=0;
    int32_t end=0;
    int hit=0;
    int miss=0;
    int32_t currentIndex=0;
    
    bool lastmatch=false;
    std::string seqMatch;

  public:

    void setMatchSeq(std::string seq);
    std::string getMatchSeq();

    void setEnd(int32_t m_end);
    int32_t getEnd();
    void setPos(int32_t pos);
    void setHit(int hit);
    void setMiss(int miss);
    void setCurrentIndex(int32_t currentIndex);
    void increaseOneOnHit();
    void setLastMatch(bool value);
    bool isLastMatch();
    void decreaseOneOnHit();
    void decreaseOnHit(int number);
    int32_t getPos();
    int getHit();
    int getMiss();
    int32_t getCurrentIndex();

    BwtScore();
};

#endif