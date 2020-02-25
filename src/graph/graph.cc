//
// Created by Sakkayaphab Piwluang on 21/2/20.
//

#include "graph.h"

Graph::Graph() {

}

void Graph::buildGraph(std::string text, int kmer) {
    std::vector<std::string> mString;
    for (int i=0;i<text.size()-kmer+1;i++) {
        std::string r = text.substr(i, kmer);
        mString.push_back(r);
    }

    for (int i=0;i<mString.size()-1;i++) {
        addNodeToNode(mString.at(i),mString.at(i+1));
    }
}

void Graph::addNode(std::string addText) {
    long long nodeindex = mNode[addText];
    long long vNodeIndexSize = vNode.size();
    if (vNodeIndexSize==0) {

    } else {
        if (getNodeByIndex(nodeindex)->getName()==addText) {
            return;
        }
    }


    Node n;
    n.setName(addText);
    vNode.push_back(n);
    mNode[addText] = vNode.size()-1;

    return;
}

void Graph::addEdge(std::string fromText, std::string toText) {
    long long nodeFromIndex = mNode[fromText];
    std::string nodeFromName = getNameNodeByIndex(nodeFromIndex);
    if (nodeFromName=="") {
        std::cerr << "nodeFrom is NULL" << std::endl;
        exit(EXIT_FAILURE);
    }

    long long nodeToIndex = mNode[toText];
    std::string nodeToName = getNameNodeByIndex(nodeToIndex);
    if (nodeToName=="") {
        std::cerr << "nodeTo is NULL" << std::endl;
        exit(EXIT_FAILURE);
    }

    Edge e;
    e.setFromNodeIndex(nodeFromIndex);
    e.setToNodeIndex(nodeToIndex);
    Node *n = getNodeByIndex(nodeFromIndex);
    n->addEdgeOut(e);
}

void Graph::showAllNode() {
    long long count = 0;
    for (Node n:vNode) {
        std::cout << "#" << n.getName()
                  << "->" << count
                  << " N.Edge: " << n.getEdges()->size() << std::endl;
        count++;
    }
}

void Graph::addNodeToNode(std::string fromNodeName, std::string toNodeName) {
    Graph::addNode(fromNodeName);
    Graph::addNode(toNodeName);
    Graph::addEdge(fromNodeName,toNodeName);
}

std::string Graph::getNameNodeByIndex(long long index) {
    if (vNode.size()>=index) {
        return vNode[index].getName();
    }
    return "";
}

Node *Graph::getNodeByIndex(long long index) {
    if (vNode.size()>=index) {
        return &vNode[index];
    }
    return NULL;
}

std::vector<std::stack<std::string>> Graph::findDFS(std::string begin,std::string end) {
    long long indexNodeCurrent = mNode[begin];
    std::string currentNodeName = getNodeByIndex(indexNodeCurrent)->getName();

    std::vector<std::stack<std::string>> output;
    if (currentNodeName!=begin) {
        return output;
    }
    std::stack<Edge*> stackEdgePath;
    std::stack<std::string> stackPath;
    stackPath.push(begin);
    bool pop = false;
    for (;;) {
        bool stop = RunNextNode(&stackPath,&stackEdgePath,&currentNodeName);
        if (currentNodeName==end) {
            showStack(stackPath);
            output.push_back(stackPath);
            bool poppass = popAllStack(&stackPath,&stackEdgePath,&currentNodeName);
            if (poppass== false) {
                break;
            }
        }

        if (stop== false) {
            bool poppass = popAllStack(&stackPath, &stackEdgePath, &currentNodeName);
            if (poppass == false) {
                break;
            }
        }

    }

    return output;
}


bool Graph::RunNextNode(std::stack<std::string> *stackPath,std::stack<Edge*> *stackEdgePath, std::string *currentNodeName) {
    long long indexNodeCurrent = mNode[*currentNodeName];
    Node *currentNode = getNodeByIndex(indexNodeCurrent);
    if (currentNode->getNumberOfUnmaskedEdges()==0) {
        return false;
    }

    long long nextNodeIndex = currentNode->getNextNodeIndexWithUnmaskedEdgeAndSetMaskedEdge(stackEdgePath);
    stackPath->push(getNodeByIndex(nextNodeIndex)->getName());
    *currentNodeName = getNodeByIndex(nextNodeIndex)->getName();

//    std::cout << "Next => " << *currentNodeName << std::endl;

    return true;
}

void Graph::reverseVectorString(std::vector<std::string> *s) {
    std::reverse((*s).begin(), (*s).end());
}

void Graph::showVectorString(std::vector<std::string> *s) {
    std::cout << getTextFromVectorString(s) << std::endl;
}

std::string Graph::getTextFromVectorString(std::vector<std::string> *s) {
    std::string output = "";
    for (long long i=0;i<s->size();i++) {
        if ((*s).at(i).size()==1) {
            output += (*s).at(i);
        } else {
            if (output=="") {
                output += (*s).at(i);
            } else {
                output += (*s).at(i).substr((*s).at(i).size()-1);
            }

           if (i+1>(*s).at(i).size()) {
               break;
           }
        }
    }

    return output;
}

void Graph::clearMaskedEdges() {
    for (long long i=0;i<vNode.size();i++) {
        vNode.at(i).clearMaskedEdges();
    }
}

std::vector<std::string> Graph::covertStackToVector(std::stack<std::string> stackPath) {

    std::vector<std::string> vResult;
    while (!stackPath.empty()) {
        vResult.push_back(stackPath.top());
        stackPath.pop();
    }

    return vResult;
}

void Graph::findEulerPath(std::stack<std::string> stackpath) {

}

void Graph::setEdgeByThisPath(std::vector<std::string> vPath) {
    for (long long i=0;i<vPath.size()-1;i++) {
        std::string beginNameNode = vPath.at(i);
        std::string endNameNode = vPath.at(i+1);
        Node *node = getNodeByIndex(mNode[beginNameNode]);
//        std::cout << "-------------------" << std::endl;
//        std::cout << beginNameNode << std::endl;
//        std::cout << endNameNode << std::endl;
        bool masked = node->setMaskedEdgeToNodeIndexWithUnmasked(mNode[endNameNode]);
        if (masked == false) {
            std::cout << "setEdgeByThisPath error" << std::endl;
        }
    }
}

void Graph::showStack(std::stack<std::string> stackPath) {
    while (!stackPath.empty()) {
        std::cout << stackPath.top() << " , ";
        stackPath.pop();
    }
    std::cout << std::endl;
}

bool Graph::popAllStack(std::stack<std::string> *stackPath, std::stack<Edge *> *stackEdgePath,std::string *currentNodeName) {

    stackPath->pop();
    if (stackPath->size()==0) {
        return false;
    }

    if (stackPath->top()==getNodeByIndex(stackEdgePath->top()->getFromNodeIndex())->getName()) {

    } else {
        while (!stackPath->empty()) {
            if (stackPath->top()==getNodeByIndex(stackEdgePath->top()->getFromNodeIndex())->getName()) {
                break;
            }
            stackEdgePath->top()->setMasked(false);
            stackEdgePath->pop();
        }
    }

    *currentNodeName = stackPath->top();

    return true;
}
