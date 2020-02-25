//
// Created by Sakkayaphab Piwluang on 23/2/20.
//

#ifndef BOLT_NODE_H
#define BOLT_NODE_H

#include <string>
#include <vector>
#include "edge.h"
#include <stack>

class Node {
    std::string name;
    std::vector<Edge> vEdge;

public:
    Node();
    void setName(std::string name);
    std::string getName();
    void addEdgeOut(Edge edge);
    std::vector<Edge> *getEdges();
    long long getNumberOfEdges();
    long long getNumberOfUnmaskedEdges();
    long long getNextNodeIndexWithUnmaskedEdgeAndSetMaskedEdge(std::stack<Edge*> *stackEdgePath);
    void clearMaskedEdges();
    bool setMaskedEdgeToNodeIndexWithUnmasked(long long nodeindex);
};


#endif //BOLT_NODE_H
