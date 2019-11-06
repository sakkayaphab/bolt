#ifndef STRINGSEARCH_H
#define STRINGSEARCH_H
#include <string>
#include <map>
#include <vector>
#include <deque> 
#include "stringsearchconfig.h"
class StringSearch
{
private:
    std::string reference;
    std::map<std::string, std::vector<int32_t>> map;
    int mapsplitsize = 3;

public:
    struct Score
    {
        int32_t pos = 0;
        int32_t end = 0;
        int32_t posseq = 0;
        int32_t endseq = 0;
        std::deque <char> matchSeqPattern;
        int32_t matchCount = 0;
        int32_t missmatchCount = 0;
        int32_t insertCount = 0;
        int32_t deleteCount = 0;
    };

    struct KeepPos
    {
        int32_t pos = 0;
        int32_t end = 0;
        int32_t count = 0;
    };

    StringSearch();
    void setReference(const char *reference);
    void buildHashTable();
    std::vector<StringSearch::Score> searchStartToEnd(std::string *seq,StringSearchConfig *ssc);
    std::vector<StringSearch::Score> searchEndToStart(std::string *seq,StringSearchConfig *ssc);
    std::vector<int> findMap(std::string seq);
    std::string getStringRange(std::string *seq, int32_t pos, int32_t range);
    std::vector<Score> runNext(std::string *seq, int32_t pos, int32_t posseq,StringSearchConfig *ssc);
    std::vector<Score> runBack(std::string *seq, int32_t pos, int32_t posseq,StringSearchConfig *ssc);
    void calculateScoreRunNext(Score *score);
    void calculateScoreRunBack(Score *score);
    void removeMissMatchEndPattern(Score *score);
    void removeMissMatchStartPattern(Score *score);
    bool conditionSaveResult(Score *score);
    StringSearch::Score getBestMatchFragment(std::string *seq);
};

#endif