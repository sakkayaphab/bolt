//
// Created by Sakkayaphab Piwluang on 21/2/20.
//

#include "graph.h"

Graph::Graph() {

}

void Graph::showGraph() {
    for (Node n:nodes) {
        std::cout << "# " << n.getText() << std::endl;
        std::cout << "out : ";
        std::vector<Edge> *edges = n.getEdges();
        for (Edge en:*edges) {
            std::cout << nodes.at(en.getEdgeTo()).getText() << " , ";
        }
        std::cout << std::endl;
        std::cout << "------------------" << std::endl;
    }
}

void Graph::addNodeToNode(std::string nodeText1, std::string nodeText2) {
    int pos1 = addNodeText(nodeText1);
    int pos2 = addNodeText(nodeText2);
    Edge e;
    e.setEdgeTo(pos2);
    nodes.at(pos1).addEdge(e);
}

int Graph::addNodeText(std::string nodetext) {
    int pos = Graph::getNodePosByText(nodetext);
    if (pos>=0) {
        return pos;
    }

    Node node;
    node.setText(nodetext);
    nodes.push_back(node);
    return nodes.size()-1;
}

int Graph::getNodePosByText(std::string nodetext) {
    for (int i=0;i<nodes.size();i++) {
        if (nodetext==nodes[i].getText()) {
            return i;
        }
    }

    return -1;
}

void Graph::findDFS() {
    std::string starttext = "1";
    int postart = getNodePosByText(starttext);
    std::cout << ">> " << nodes.at(postart).getText() << std::endl;

    int count = 0;
    std::stack<std::string> stack;
//    stack.push(nodes.at(postart).getText());
    for (;;) {
        count++;
        int getnodetoedgeindex = nodes.at(postart).getEdgeIndexWithUnmasked();
        if (getnodetoedgeindex<0) {
            break;
        }

        stack.push(nodes.at(postart).getText());
        if (nodes.at(postart).getText()=="6") {
            break;
        }

        int getedgetonodeindex = nodes.at(postart).getEdges()->at(getnodetoedgeindex).getEdgeTo();
        std::cout << "getedgeto : " << getnodetoedgeindex << " at node : " << nodes.at(getedgetonodeindex).getText() << std::endl;
        nodes.at(postart).getEdges()->at(getnodetoedgeindex).setMarked(true);
        postart = getedgetonodeindex;
    }

    showstack(stack);

    std::cout << "=== edge remain ===" << std::endl;
    for (Node n:nodes) {
        if (n.getEdgeIndexWithUnmasked()>=0) {
            std::cout << "# " << n.getText() << std::endl;
            std::cout << "out : " << n.getEdgeNumberUnmasked() << std::endl;
            std::cout << "------------------" << std::endl;
        }
    }
}

void Graph::showstack(std::stack <std::string> s)
{
    while (!s.empty())
    {
        std::cout << ',' << s.top();
        s.pop();
    }
    std::cout << '\n';
}

int Graph::getNumberEdgeInNode(std::string text) {
    int posNode = getNodePosByText(text);
    if (posNode<0) {
        return -1;
    }

    int count = 0;
    for (Node n:nodes) {
        for (Edge e:*n.getEdges()) {
            if (e.getEdgeTo()==posNode) {
                count++;
            }
        }
    }

    return count;
}

int Graph::getNumberEdgeOutNode(std::string text) {
    int posNode = getNodePosByText(text);
    if (posNode<0) {
        return -1;
    }

    Node node = nodes.at(posNode);
    return node.getEdges()->size();
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



