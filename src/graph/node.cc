//
// Created by Sakkayaphab Piwluang on 23/2/20.
//

#include "node.h"

Node::Node() {

}

void Node::setName(std::string name) {
    Node::name = name;
}

std::string Node::getName() {
    return name;
}

void Node::addEdgeOut(Edge edge) {
    vEdge.push_back(edge);
}

std::vector<Edge> *Node::getEdges() {
    return &vEdge;
}

long long Node::getNumberOfEdges() {
    return vEdge.size();
}

long long Node::getNumberOfUnmaskedEdges() {
    long long count = 0;
    for (Edge e:vEdge) {
        if (e.getMasked()==false) {
            count++;
        }
    }

    return count;
}

long long Node::getNextNodeIndexWithUnmaskedEdgeAndSetMaskedEdge(std::stack<Edge*> *stackEdgePath) {
    for (int i=0;i<vEdge.size();i++) {
        if (vEdge.at(i).getMasked()== false) {
            vEdge.at(i).setMasked(true);
            stackEdgePath->push(&vEdge.at(i));
            return vEdge.at(i).getToNodeIndex();
        }
    }
    return -1;
}

void Node::clearMaskedEdges() {
    for (long long i=0;i<vEdge.size();i++) {
        vEdge.at(i).setMasked(false);
    }
}

bool Node::setMaskedEdgeToNodeIndexWithUnmasked(long long nodeindex) {
    for (long long i=0;i<vEdge.size();i++) {
        if (vEdge.at(i).getMasked()== true) {
            continue;
        }

        if (vEdge.at(i).getToNodeIndex()==nodeindex) {
            vEdge.at(i).setMasked(true);
            return true;
        }
    }

    return false;
}
