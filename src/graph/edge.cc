//
// Created by Sakkayaphab Piwluang on 22/2/20.
//

#include "edge.h"

Edge::Edge() {

}

void Edge::setEdgeTo(int pos) {
    edgeto = pos;
}

int Edge::getEdgeTo() {
    return edgeto;
}

void Edge::setMarked(bool marked) {
    Edge::marked = marked;
}

bool Edge::getMarked() {
    return Edge::marked;
}
