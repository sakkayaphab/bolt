//
// Created by Sakkayaphab Piwluang on 21/2/20.
//

#include "node.h"
#include <vector>

Node::Node() {

}

void Node::setText(std::string text) {
    Node::text = text;
}

std::string Node::getText() {
    return text;
}

void Node::addEdge(Edge edge) {
    edges.push_back(edge);
}

std::vector<Edge> *Node::getEdges() {
    return &edges;
}

void Node::clearEdges() {
    for (int i=0;i<edges.size();i++) {
        edges.at(i).setMarked(false);
    }
}

int Node::getEdgeNumberUnmasked() {
    int count = 0;
    for (int i=0;i<edges.size();i++) {
        if (edges.at(i).getMarked()== false) {
            count += 1;
        }
    }
    return count;
}

int Node::getEdgeIndexWithUnmasked() {
    for (int i=0;i<edges.size();i++) {
        if (edges.at(i).getMarked()==false) {
            return i;
        }
    }
    return -1;
}