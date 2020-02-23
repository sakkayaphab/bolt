//
// Created by Sakkayaphab Piwluang on 23/2/20.
//

#ifndef BOLT_EDGE_H
#define BOLT_EDGE_H

#include <string>

class Edge {
    bool masked = false;
    long long fromNodeIndex;
    long long toNodeIndex;

public:
    Edge();
    void setMasked(bool masked);
    bool getMasked();
    void setFromNodeIndex(long long fromNodeIndex);
    void setToNodeIndex(long long toNodeIndex);
    long long getFromNodeIndex();
    long long getToNodeIndex();

};


#endif //BOLT_EDGE_H
