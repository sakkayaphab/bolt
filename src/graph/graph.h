//
// Created by Sakkayaphab Piwluang on 21/2/20.
//

#ifndef BOLT_GRAPH_H
#define BOLT_GRAPH_H

#include <string>
#include <vector>
#include <stack>
#include <iostream>
#include <map>
#include "edge.h"
#include "node.h"
#include <algorithm>

class Graph {
    std::vector<Node> vNode;
    std::map<std::string,long long> mNode;

public:
    Graph();
    void buildGraph(std::string text,int kmer);
    void addNode(std::string addText);
    void addEdge(std::string fromText,std::string toText);
    void addNodeToNode(std::string fromNodeName,std::string toNodeName);
    std::string getNameNodeByIndex(long long index);
    Node *getNodeByIndex(long long index);
    void showAllNode();
    std::vector<std::stack<std::string>> findDFS(std::string begin,std::string end);
    void reverseVectorString(std::vector<std::string> *s);
    void showVectorString(std::vector<std::string> *s);
    std::string getTextFromVectorString(std::vector<std::string> *s);
    bool RunNextNode(std::stack<std::string> *stackPath,std::stack<Edge*> *stackEdgePath, std::string *currentNodeName);
    void clearMaskedEdges();
    std::vector<std::string> covertStackToVector(std::stack<std::string> mStack);
    void setEdgeByThisPath(std::vector<std::string> vPath);
    void findEulerPath(std::stack<std::string> stackpath);
    void showStack(std::stack<std::string> mStack);
    bool popAllStack(std::stack<std::string> *stackPath,std::stack<Edge*> *stackEdgePath, std::string *currentNodeName);
};

#endif //BOLT_GRAPH_H
