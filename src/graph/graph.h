//
// Created by Sakkayaphab Piwluang on 21/2/20.
//

#ifndef BOLT_GRAPH_H
#define BOLT_GRAPH_H

#include <string>
#include "node.h"
#include <vector>
#include <stack>
#include <iostream>

class Graph {
    std::vector<Node> nodes;

public:
    Graph();
    void addNodeToNode(std::string nodetext1,std::string nodetext2);
    int addNodeText(std::string nodetext);
    int getNodePosByText(std::string nodetext);
    void showGraph();
    void findDFS();
    int getNumberEdgeInNode(std::string text);
    int getNumberEdgeOutNode(std::string text);
    void showstack(std::stack <std::string> s);
    void buildGraph(std::string text,int kmer);
};


#endif //BOLT_GRAPH_H
