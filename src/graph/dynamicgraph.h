//
// Created by Sakkayaphab Piwluang on 11/3/20.
//

#ifndef BOLT_DYNAMICGRAPH_H
#define BOLT_DYNAMICGRAPH_H

#include <string>
#include <vector>
#include "graphresult.h"
#include "graph.h"

class DynamicGraph {
private:
    std::vector<std::string> vRawTextFirst;
    std::vector<std::string> vRawTextSecond;
    int kmer = 0;
public:
    DynamicGraph();
    void buildGraph(std::string text);
    void addNodeToNode(std::string fromNodeName, std::string toNodeName);
    std::string reverseString(std::string text);
    void setKmer(int kmer);
    int getKmer();
    GraphResult getGraphResult(std::string begin,std::string end);
    GraphResult ReverseStringGraphResult(GraphResult *gr);
};

#endif //BOLT_DYNAMICGRAPH_H
