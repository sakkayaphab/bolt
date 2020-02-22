//
// Created by Sakkayaphab Piwluang on 21/2/20.
//

#ifndef BOLT_NODE_H
#define BOLT_NODE_H

#include <string>
#include <map>
#include "edge.h"

class Node {
    std::string text;
    std::vector<Edge> edges;

public:
    Node();
    void setText(std::string text);
    std::string getText();
    void addEdge(Edge edge);
    std::vector<Edge> *getEdges();
    void clearEdges();
    int getEdgeNumberUnmasked();
    int getEdgeIndexWithUnmasked();
};

#endif //BOLT_NODE_H
