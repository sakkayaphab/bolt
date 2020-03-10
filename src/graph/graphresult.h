//
// Created by Sakkayaphab Piwluang on 11/3/20.
//

#ifndef BOLT_GRAPHRESULT_H
#define BOLT_GRAPHRESULT_H

#include <string>

class GraphResult {
private:
    std::string maxLeft;
    std::string maxRight;
    std::string maxConcordant;
public:
    GraphResult();
    std::string getMaxLeft();
    std::string getMaxRight();
    std::string getMaxConcordant();
    void setMaxLeft(std::string text);
    void setMaxRight(std::string text);
    void setMaxConcordant(std::string text);
};


#endif //BOLT_GRAPHRESULT_H
