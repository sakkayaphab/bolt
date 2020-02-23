//
// Created by Sakkayaphab Piwluang on 23/2/20.
//

#include "edge.h"

Edge::Edge() {

}

void Edge::setMasked(bool masked) {
    Edge::masked = masked;
}

bool Edge::getMasked() {
    return Edge::masked;
}

void Edge::setFromNodeIndex(long long fromNodeIndex) {
    Edge::fromNodeIndex = fromNodeIndex;
}

void Edge::setToNodeIndex(long long toNodeIndex) {
    Edge::toNodeIndex = toNodeIndex;
}

long long Edge::getFromNodeIndex() {
    return Edge::fromNodeIndex;
}

long long Edge::getToNodeIndex() {
    return Edge::toNodeIndex;
}

