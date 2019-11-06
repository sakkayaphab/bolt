#include "stringsearch.h"
#include <iostream>
#include <string.h>
#include <algorithm>
StringSearch::StringSearch()
{
}

void StringSearch::setReference(const char *reference)
{
    StringSearch::reference = reference;
}

void StringSearch::buildHashTable()
{
    map.clear();
    // std::cout << reference.length() << std::endl;
    std::string temp;
    for (int32_t i = 0; i < reference.length() - mapsplitsize + 1; i++)
    {
        temp = "";
        for (int j = 0; j < mapsplitsize; j++)
        {
            temp += reference[i + j];
        }
        map[temp].push_back(i);
        // std::cout << temp << std::endl;
    }

    // std::cout << "-----------MAP-----------" << std::endl;
    // for (auto n : map)
    // {
    //     std::cout << n.first << " = ";
    //     for (auto m : n.second)
    //     {
    //         std::cout << m;
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << "------------------------" << std::endl;
}

std::vector<int32_t> StringSearch::findMap(std::string seq)
{
    return map[seq];
}

std::string StringSearch::getStringRange(std::string *seq, int32_t pos, int32_t range)
{
    return seq->substr(pos, range);
}

std::vector<StringSearch::Score> StringSearch::searchStartToEnd(std::string *seq, StringSearchConfig *ssc)
{
    std::vector<StringSearch::Score> scoreresult;
    if (seq->length() < mapsplitsize)
    {
        return scoreresult;
    }

    std::string data = getStringRange(seq, 0, mapsplitsize);

    // std::cout << "-----------INDEX-----------" << std::endl;
    std::vector<int32_t> index = findMap(data);

    for (auto n : index)
    {
        std::vector<StringSearch::Score> scorelist = runNext(seq, n, 0, ssc);
        for (auto m : scorelist)
        {
            scoreresult.push_back(m);
        }
    }

    return scoreresult;
}

std::vector<StringSearch::Score> StringSearch::searchEndToStart(std::string *seq, StringSearchConfig *ssc)
{
    std::vector<StringSearch::Score> scoreresult;
    if (seq->length() < mapsplitsize)
    {
        return scoreresult;
    }

    std::string data = getStringRange(seq, seq->length() - mapsplitsize, mapsplitsize);

    // std::cout << "-----------INDEX-----------" << data << std::endl;
    std::vector<int32_t> index = findMap(data);

    for (auto n : index)
    {
        std::vector<StringSearch::Score> scorelist = runBack(seq, n + mapsplitsize - 1, seq->size() - mapsplitsize - 1, ssc);
        for (auto m : scorelist)
        {
            scoreresult.push_back(m);
        }
    }

    return scoreresult;
}

bool StringSearch::conditionSaveResult(Score *score)
{
    if (score->matchCount <= 3)
    {
        return false;
    }

    if (score->matchCount <= 6)
    {
        if (score->missmatchCount >= 2)
        {
            return false;
        }
        return true;
    }

    if (score->matchCount > 6)
    {
        return true;
    }

    return false;
}

void StringSearch::calculateScoreRunBack(StringSearch::Score *score)
{
    removeMissMatchStartPattern(score);
    // for (auto a : score->matchSeqPattern)
    // {
    //     std::cout << a;
    // }

    int32_t i = score->matchSeqPattern.size() - 1;
    score->pos = score->end;

    int32_t tempposseq = 0;

    for (; i >= 0; i--)
    {
        // std::cout << score->matchSeqPattern[i] << " = " <<
        if (score->matchSeqPattern[i] == 'M')
        {
            score->pos--;
            score->posseq--;
            score->matchCount++;
        }
        else if (score->matchSeqPattern[i] == 'X')
        {
            score->pos--;
            score->posseq--;
            score->missmatchCount++;
        }
    }
}

void StringSearch::removeMissMatchStartPattern(StringSearch::Score *score)
{
    int32_t i = 0;

    for (; i < score->matchSeqPattern.size(); i++)
    {
        if (score->matchSeqPattern[i] != 'M')
        {
            score->matchSeqPattern.pop_front();
            i--;
        }
        else
        {
            break;
        }
    }
}

std::vector<StringSearch::Score> StringSearch::runBack(std::string *seq, int32_t pos, int32_t posseq, StringSearchConfig *ssc)
{
    // std::cout << "-----------runBack-----------" << std::endl;
    std::vector<StringSearch::Score> scorelist;
    Score score;
    score.end = pos + 1;
    for (int j = 0; j < mapsplitsize; j++)
    {
        score.matchSeqPattern.push_front('M');
    }

    int32_t i = pos - mapsplitsize;
    int32_t iseq = posseq;
    int32_t lastContinueMissMatch = 0;
    int32_t matchCount = 0;
    int32_t accMissMatch = 0;
    for (; i >= 0; i--)
    {
        // std::cout << reference[i] << " = " << (*seq)[iseq] << std::endl;
        if (reference[i] == (*seq)[iseq])
        {
            score.matchSeqPattern.push_front('M');
            lastContinueMissMatch = 0;
            matchCount++;
        }
        else
        {

            if (lastContinueMissMatch + 1 >= ssc->getAllowMissMatch() || accMissMatch + 1 >= ssc->getAllowMissMatch())
            {
            exitlabel:
                for (int j = 0; j < lastContinueMissMatch; j++)
                {
                    score.matchSeqPattern.pop_front();
                }
                break;
            }
            score.matchSeqPattern.push_front('X');
            lastContinueMissMatch++;
            accMissMatch++;
            matchCount++;
        }

        //condition loop
        iseq--;
        if (iseq < 0)
        {
            break;
        }

        if (ssc->getMaxAllowAlign() < matchCount && ssc->getMaxAllowAlign() != 0)
        {
            goto exitlabel;
        }
    }

    score.endseq = seq->size();
    score.posseq = score.endseq;
    calculateScoreRunBack(&score);

    if (conditionSaveResult(&score))
    {
        scorelist.push_back(score);
    }

    return scorelist;
}

std::vector<StringSearch::Score> StringSearch::runNext(std::string *seq, int32_t pos, int32_t posseq, StringSearchConfig *ssc)
{
    // std::cout << "-----------runNext-----------" << std::endl;

    Score score;
    score.pos = pos;

    int32_t i = pos + mapsplitsize;
    for (int j = 0; j < mapsplitsize; j++)
    {
        score.matchSeqPattern.push_back('M');
    }

    int32_t iseq = posseq + mapsplitsize;
    int32_t lastContinueMissMatch = 0;
    int32_t matchCount = 0;
    int32_t accMissMatch = 0;
    for (; i < reference.length(); i++)
    {
        if (reference[i] == (*seq)[iseq])
        {
            score.matchSeqPattern.push_back('M');
            matchCount++;
            lastContinueMissMatch = 0;
        }
        else
        {

            if (lastContinueMissMatch + 1 > ssc->getAllowMissMatch() || accMissMatch + 1 > ssc->getAllowMissMatch())
            {
            exitlabel:
                for (int j = 0; j < lastContinueMissMatch; j++)
                {
                    score.matchSeqPattern.pop_back();
                }

                break;
            }
            score.matchSeqPattern.push_back('X');
            lastContinueMissMatch++;
            accMissMatch++;
            matchCount++;
        }

        // std::cout << (*seq)[iseq] << iseq << std::endl;
        //condition loop
        iseq++;
        if (iseq > (seq->size() - 1))
        {
            break;
        }

        if (ssc->getMaxAllowAlign() < matchCount && ssc->getMaxAllowAlign() != 0)
        {
            goto exitlabel;
        }
    }
    // for (auto a : score.matchSeqPattern)
    // {
    //     std::cout << a;
    // }
    // std::cout << std::endl;
    calculateScoreRunNext(&score);

    // std::cout << std::endl;
    std::vector<StringSearch::Score> scorelist;
    scorelist.push_back(score);

    return scorelist;
}

void StringSearch::calculateScoreRunNext(StringSearch::Score *score)
{
    removeMissMatchEndPattern(score);
    int32_t tempend = 0;
    int32_t i = 0;

    for (; i < score->matchSeqPattern.size(); i++)
    {
        // std::cout << score->matchSeqPattern[i];
        if (score->matchSeqPattern[i] == 'M')
        {
            score->endseq++;
            tempend++;
            score->matchCount++;
        }
        else if (score->matchSeqPattern[i] == 'X')
        {
            score->endseq++;
            tempend++;
            score->missmatchCount++;
        }
    }
    // std::cout << std::endl;

    score->end = score->pos + tempend;
}

void StringSearch::removeMissMatchEndPattern(StringSearch::Score *score)
{
    int32_t i = score->matchSeqPattern.size() - 1;

    for (; i >= 0; i--)
    {
        if (score->matchSeqPattern[i] != 'M')
        {
            score->matchSeqPattern.pop_back();
        }
        else
        {
            break;
        }
    }
}

StringSearch::Score StringSearch::getBestMatchFragment(std::string *seq)
{
    std::string temp = "";
    std::vector<KeepPos> keeplist;
    int count = 0;
    for (int32_t i = 0; i < seq->length() - mapsplitsize + 1; i++)
    {
        temp = "";
        for (int j = 0; j < mapsplitsize; j++)
        {
            temp += (*seq)[i + j];
        }

        std::vector<int> tempMatch = findMap(temp);
        for (auto n : tempMatch)
        {
            KeepPos tempkeep;
            tempkeep.pos = n;
            tempkeep.count = count;
            tempkeep.end = tempkeep.pos+mapsplitsize;
            keeplist.push_back(tempkeep);
            std::cout << n << std::endl;
        }
        std::cout << " ----- " << std::endl;
        // MatchList.insert(MatchList.end(), tempMatch.begin(), tempMatch.end());
        
        if (count>=5) {
             break;
        }
        count++;
       
    }

    // std::sort(MatchList.begin(), MatchList.end());
    // MatchList.erase(std::unique(MatchList.begin(), MatchList.end()), MatchList.end());

    // for (auto n : MatchList)
    // {
    //     std::cout << n << std::endl;
    // }
}