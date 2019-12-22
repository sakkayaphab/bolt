#include <iostream>
#include "core.h"
#include <algorithm>
#include <sstream>

void Core::setMapLast()
{

    BwtStore *p;
    std::map<std::string, int> rankOrder;
    std::string temp;

    for (int i = 0; i < bwtstore.size(); i++)
    {
        temp = "";
        p = &bwtstore.at(i);
        mapRankLast[(*p).last].push_back(i);

        temp += (*p).last;
        std::stringstream ss;
        ss << (*p).lastindex;
        temp += ss.str();

        rankOrder[temp] = i;

    }

    // return;

    for (int i = 0; i < bwtstore.size(); i++)
    {
        temp = "";
        p = &bwtstore.at(i);
        temp += (*p).first;
        std::stringstream ss;
        ss << (*p).firstindex;
        temp += ss.str();
        (*p).nexttolast = rankOrder[temp];
        // std::cout << temp << "=" << (*p).lastindex << "[" << rankOrder[temp] << "]" << std::endl;

        //a0=0[1]
        // b0=0[5]
        // b1=1[6]
        // a1=1[2]
        // $0=0[0]
        // a2=2[3]
        // a3=3[4]
    }
}

void Core::setMapFirst()
{
    BwtStore *p;
    std::map<std::string, int> rankOrder;
    std::string temp;

    for (int i = 0; i < bwtstore.size(); i++)
    {
        temp = "";
        p = &bwtstore.at(i);
        mapRankFirst[(*p).first].push_back(i);

        //for make map
        temp += (*p).first;
        std::stringstream ss;
        ss << (*p).firstindex;
        temp += ss.str();

        rankOrder[temp] = i;

        // std::cout << temp << "=" << i << std::endl;
        // std::cout << mapRankFirst[p.first].size() << std::endl;
    }

    // std::cout << "[" << rankOrder["a0"] << "]" << std::endl;
    // return;

    //for make map
    for (int i = 0; i < bwtstore.size(); i++)
    {
        temp = "";
        p = &bwtstore.at(i);
        temp += (*p).last;
        std::stringstream ss;
        ss << (*p).lastindex;
        temp += ss.str();
        (*p).lasttonext = rankOrder[temp];
        // std::cout << temp << "=" << (*p).lastindex << "[" << rankOrder[temp] << "]" << std::endl;

        //a0=0[1]
        // b0=0[5]
        // b1=1[6]
        // a1=1[2]
        // $0=0[0]
        // a2=2[3]
        // a3=3[4]
    }
}

std::vector<int> Core::getIndexFromMapFirst(std::string *str_m)
{
    return mapRankFirst[(*str_m).at(str_m->size() - 1)];
}

std::vector<int> Core::getIndexFromMapLast(std::string *str_m)
{
    return mapRankLast[(*str_m).at(0)];
}

// bool Core::compareString(const RotationString &a, const RotationString &b)
// {
//     return a.text < b.text;
// }

void Core::rotataions(const char *seq)
{
    std::string buf(seq);
    int32_t bufsize = buf.size();
    buf = buf + "$" + buf;

    std::vector<RotationString> myvector;

    RotationString rs;
    for (int i = 0; i < bufsize + 1; i++)
    {
        std::string out = buf.substr(i, buf.size() - bufsize);
        rs.text = out;
        rs.position = i;
        myvector.push_back(rs);
    }

    // if ("ba" < "ab") {
    //     std::cout << "ba > ab" << std::endl;
    // }
    //sort string
    std::sort(myvector.begin(), myvector.end());

    // qsort()
    //     auto sortFunc = [](RotationString a, RotationString b) -> bool {
    //     return a.text < b.text;
    // };

    // std::sort(myvector.begin(), myvector.end(), sortFunc);
    // sort(myvector.begin(), myvector.end(), [](const RotationString &lhs, const RotationString &rhs) {
    //     return lhs.text < rhs.text;
    // });

    // char Bwt[bufsize+1];
    std::map<char, int> mRankFirst;
    std::map<char, int> mRankLast;
    RotationString rsF;
    std::string p;
    for (int i = 0; i < myvector.size(); i++)
    {
        rsF = myvector.at(i);
        p = rsF.text;

        //     // p.find('$');

        //     // std::cout << bufsize-p.find('$') << std::endl;

        //     //Last char
        char letterLast = p[bufsize];
        int32_t vLast = mRankLast[letterLast];
        mRankLast[letterLast]++;

        //First char
        char letterFirst = p[0];
        int32_t vFirst = mRankFirst[letterFirst];
        mRankFirst[letterFirst]++;

        BwtStore bwtstemp;
        bwtstemp.first = letterFirst;
        bwtstemp.firstindex = vFirst;
        bwtstemp.last = letterLast;
        bwtstemp.lastindex = vLast;

        bwtstemp.positionoffirst = rsF.position;

        bwtstore.push_back(bwtstemp);

        // std::cout << "F : " << letterFirst << "=" << vFirst << std::endl;
        // std::cout << "L : " << letterLast << "=" << vLast << std::endl;
        // std::cout << "> " << p << " index = " << rsF.position << std::endl;
        // std::cout << "----" << std::endl;

        //output
        // F : $=0
        // L : a=0
        // > $abaaba index = 6
        // ----
        // F : a=0
        // L : b=0
        // > a$abaab index = 5
        // ----
        // F : a=1
        // L : b=1
        // > aaba$ab index = 2
        // ----
        // F : a=2
        // L : a=1
        // > aba$aba index = 3
        // ----
        // F : a=3
        // L : $=0
        // > abaaba$ index = 0
        // ----
        // F : b=0
        // L : a=2
        // > ba$abaa index = 4
        // ----
        // F : b=1
        // L : a=3
        // > baaba$a index = 1
        // ----
    }
}

void Core::bwtviaBwm(const char *seq)
{
    rotataions(seq);
}

void Core::build(const char *seq)
{
    bwtviaBwm(seq);
    setMapFirst();
    setMapLast();
}

void Core::setStartPosition(uint32_t position_m)
{
    startPosition = position_m;
}

Core::Core()
{
}

std::vector<BwtScore> Core::runNext(std::string *str_m, int startindex, int allowMissMatch,
                                    int allowMissMatchExtended, int minimumMatch)
{
    std::vector<BwtScore> bwtstateList;
    BwtScore bwtscore;

    int sizeStr = (*str_m).size();
    BwtStore *bwtStoreTemp;

    bwtStoreTemp = &bwtstore.at(startindex);
    char charfirst;
    int32_t indexlasttonext;
    bwtscore.setPos(startPosition + bwtstore.at(startindex).positionoffirst + 1);
    // bwtscore.setHit(0);
    // return;
    int i = sizeStr - 1;
    int missmatchcount = 0;
    int sumMissmatchcount = 0;
    for (;;)
    {
        // std::cout << ":::::" << i << "/" << sizeStr << std::endl;
        // std::cout << "[" << bwtstore.at(startindex).first << "/" << (*str_m).at(i) << " " << startindex << "] ";
        charfirst = bwtstore.at(startindex).first;
        indexlasttonext = bwtstore.at(startindex).lasttonext;
        // std::cout << "charlast : " << charfirst << " , (*str_m).at(i)) : " << (*str_m).at(i) << std::endl;

        if (charfirst == (*str_m).at(i))
        {
            bwtscore.increaseOneOnHit();
            // std::cout << ",[" << charfirst << "/" << (*str_m).at(i) << " " << indexlasttonext << "] ";
            startindex = indexlasttonext;
            missmatchcount = 0;
        }
        else if (missmatchcount < allowMissMatchExtended)
        {
            if (sumMissmatchcount >= 2)
            {
                break;
            }
            sumMissmatchcount++;

            bwtscore.increaseOneOnHit();
            // std::cout << ",[" << charfirst << "/" << (*str_m).at(i) << " " << indexlasttonext << "] ";
            startindex = indexlasttonext;
            missmatchcount++;
        }
        else
        {
            if (missmatchcount > 0)
            {
                bwtscore.decreaseOneOnHit();
            }
            // std::cout << "EXIT >>>> " << charfirst << "/" << (*str_m).at(i) << " <<<< ";
            break;
        }

        i--;
        if (i < 0)
        {
            if (missmatchcount > 0)
            {
                bwtscore.decreaseOneOnHit();
            }
            // std::cout << "EXIT if (i >= sizeStr) " << std::endl;
            break;
        }
        else if (charfirst == '$')
        {
            if (missmatchcount > 0)
            {
                bwtscore.decreaseOneOnHit();
            }
            // std::cout << "EXIT else if (charfirst == '$') " << std::endl;
            break;
        }
    }

    bwtscore.setEnd(bwtscore.getPos());
    bwtscore.setPos(bwtscore.getPos() - bwtscore.getHit());

    // start : 1 start at position : 5
    // [a/a 1] ,[b/b 5] ,[a/a 3]

    // --------- = 3
    // size of str : 3 = 3

    // start : 2 start at position : 2
    // [a/a 2] ,[b/b 6] ,[a/a 4]

    // --------- = 3
    // size of str : 3 = 3

    // start : 3 start at position : 3
    // [a/a 3]
    // >>>> a/b <<<<
    // --------- = 1
    // size of str : 3 = 3

    // start : 4 start at position : 0
    // [a/a 4]
    // >>>> $/b <<<<
    // --------- = 1

    if (bwtscore.getHit() >= minimumMatch)
    {
        bwtscore.setMatchSeq(str_m->substr(str_m->size() - bwtscore.getHit(), str_m->size()));
        bwtstateList.push_back(bwtscore);
        // std::cout << "--------- = " << (bwtscore).getHit() << std::endl;
    }

    return bwtstateList;
}

std::vector<BwtScore> Core::runBackAllowOnlyMissmatch(std::string *str_m, int startindex, int allowMissMatch, int allowMissMatchExtended, int minimumMatch)
{
    
    BwtScore bwtscore;

    BwtStore *bwtStoreTemp;

    bwtStoreTemp = &bwtstore.at(startindex);
    char charlast;
    int32_t indexlasttonext;
    bwtscore.increaseOneOnHit();

    bwtscore.setPos(startPosition + bwtstore.at(startindex).positionoffirst - 1);
    // return;
    int i = 1;

    int missmatchcount = 0;
    for (;;)
    {
        charlast = bwtstore.at(startindex).first;
        indexlasttonext = bwtstore.at(startindex).nexttolast;
        // std::cout << "charlast : " << charlast << " , (*str_m).at(i)) : " << (*str_m).at(i) << std::endl;

        if (charlast == (*str_m).at(i))
        {
            bwtscore.increaseOneOnHit();
            startindex = indexlasttonext;
            missmatchcount = 0;
        }
        else if (missmatchcount < allowMissMatchExtended)
        {
            bwtscore.increaseOneOnHit();
            startindex = indexlasttonext;
            missmatchcount++;
        }
        else
        {
            if (missmatchcount > 0)
            {
                // bwtscore.de
                bwtscore.decreaseOneOnHit();
            }
            // std::cout << "break else" << std::endl;
            break;
        }

        i++;
        if (i >= minimumMatch)
        {
            // if (missmatchcount > 0)
            // {
            //     bwtscore.decreaseOneOnHit();
            // }
            // // std::cout << "if (i >= sizeStr)" << std::endl;
            // break;
        }
        else if (charlast == '$')
        {
            if (missmatchcount > 0)
            {
                bwtscore.decreaseOneOnHit();
            }
            // std::cout << "else if (charlast == '$')" << std::endl;
            break;
        }
    }

    

    bwtscore.setEnd(bwtscore.getPos() + bwtscore.getHit());

    // start : 1 start at position : 5
    // [a/a 1] ,[b/b 5] ,[a/a 3]

    // --------- = 3
    // size of str : 3 = 3

    // start : 2 start at position : 2
    // [a/a 2] ,[b/b 6] ,[a/a 4]

    // --------- = 3
    // size of str : 3 = 3

    // start : 3 start at position : 3
    // [a/a 3]
    // >>>> a/b <<<<
    // --------- = 1
    // size of str : 3 = 3

    // start : 4 start at position : 0
    // [a/a 4]
    // >>>> $/b <<<<
    // --------- = 1
    std::vector<BwtScore> bwtstateList;
    if (bwtscore.getHit() >= minimumMatch)
    {
        bwtscore.setMatchSeq(str_m->substr(0, bwtscore.getHit()));
        bwtstateList.push_back(bwtscore);
        //  std::cout << "--------- = " << (bwtscore).getHit() << std::endl;
    }
    // std::cout << "------------------"<< std::endl;
    return bwtstateList;
}

/*
    if Bwt(abaaba)
    $ a b a a b a3          (0)
    a0 $ a b a a b1         (1)
    a1 a b a $ a b0         (2)
    a2 b a $ a b a1         (3)
    a3 b a a b a $          (4)
    b0 a $ a b a a2         (5)
    b1 a a b a $ a0         (6)
*/

/*
    Bwt = "helloworld";
    findEndToStart = "world";
    position : 9, hit : 5
*/

std::vector<BwtScore> Core::findEndToStart(std::string *str_m, int32_t allowMissMatch, int32_t allowMissMatchExtended, int32_t minimumMatch)
{
    // std::cout << "found size : " << getIndexFromMapFirst(str_m).size() << " find by : " << (*str_m).at(str_m->size() - 1) << std::endl;
    std::vector<int> p = getIndexFromMapFirst(str_m);

    std::vector<BwtScore> results;
    std::vector<BwtScore> bwtstatelist;

    for (int i = 0; i < p.size(); i++)
    {
        // std::cout << "start : " << p[i] << " start at position : " << std::endl;
        bwtstatelist = runNext(str_m, p[i], allowMissMatch, allowMissMatchExtended, minimumMatch);
        if (bwtstatelist.size() == 0)
        {
        }
        else
        {
            for (BwtScore n : bwtstatelist)
            {
                results.push_back(n);
            }
        }
    }

    return results;
}

/*
# findStartToEnd
    Bwt = "helloworld";
    findEndToStart = "world";
    position : 5, hit : 5
*/
std::vector<BwtScore> Core::findStartToEnd(std::string *str_m, int32_t allowMissMatch, int32_t allowMissMatchExtended, int32_t minimumMatch)
{
    std::vector<int> p = getIndexFromMapLast(str_m);

    std::vector<BwtScore> results;

    for (int i = 0; i < p.size(); i++)
    {
        std::vector<BwtScore> bwtstatelist;
        // std::cout << "++++++++" << std::endl;
        // std::cout << "start index : " << p[i] << " find by : " << (*str_m).at(0) << std::endl;

        bwtstatelist = runBackAllowOnlyMissmatch(str_m, p[i], allowMissMatch, allowMissMatchExtended, minimumMatch);

        for (BwtScore n : bwtstatelist)
        {
            results.push_back(n);
        }

        bwtstatelist.clear();
    }

    return results;
}

// F                L
// $ a0 b0 a1 a2 b1 a3
// a3 $ a0 b0 a1 a2 b1
// a1 a2 b1 a3 $ a0 b0
// a2 b1 a3 $ a0 b0 a1
// a0 b0 a1 a2 b1 a3 $
// b1 a3 $ a0 b0 a1 a2
// b0 a1 a2 b1 a3 $ a0

void Core::printBwtStore()
{
    //     char first;
    // uint32_t firstindex;
    // char last;
    // uint32_t lastindex;
    // uint32_t lasttonext;
    // uint32_t nexttolast;
    // uint32_t positionoffirst;

    for (BwtStore n : bwtstore)
    {
        std::cout
            << "[[ "
            << "first : " << n.first << ","
            << "last : " << n.last << ","
            << "firstindex : " << n.firstindex << ","
            << "lastindex : " << n.lastindex << ","
            << "lasttonext : " << n.lasttonext << ","
            << "nexttolast : " << n.nexttolast << ","
            << "positionoffirst : " << n.positionoffirst << ","
            << " ]]"
            << std::endl;
    }
    // bwtstore
}
