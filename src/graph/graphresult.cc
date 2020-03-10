//
// Created by Sakkayaphab Piwluang on 11/3/20.
//

#include "graphresult.h"

GraphResult::GraphResult() {

}

std::string GraphResult::getMaxLeft() {
    return GraphResult::maxLeft;
}

std::string GraphResult::getMaxRight() {
    return GraphResult::maxRight;
}

std::string GraphResult::getMaxConcordant() {
    return GraphResult::maxConcordant;
}

void GraphResult::setMaxLeft(std::string text) {
    if (text.size()>GraphResult::maxLeft.size()) {
        GraphResult::maxLeft = text;
    }
}

void GraphResult::setMaxRight(std::string text) {
    if (text.size()>GraphResult::maxRight.size()) {
        GraphResult::maxRight = text;
    }
}

void GraphResult::setMaxConcordant(std::string text) {
    if (text.size()>GraphResult::maxConcordant.size()) {
        GraphResult::maxConcordant = text;
    }
}
