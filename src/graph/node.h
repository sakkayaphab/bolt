//
// Created by Sakkayaphab Piwluang on 23/2/20.
//

#ifndef BOLT_NODE_H
#define BOLT_NODE_H

#include <string>
#include <vector>
#include "edge.h"

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
    long long getNextNodeIndexWithUnmaskedEdgeAndSetMaskedEdge();
    void clearMaskedEdges();

};


#endif //BOLT_NODE_H
