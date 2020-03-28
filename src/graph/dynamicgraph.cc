//
// Created by Sakkayaphab Piwluang on 11/3/20.
//

#include "dynamicgraph.h"

DynamicGraph::DynamicGraph() {

}

void DynamicGraph::buildGraph(std::string text) {
    if (kmer==0) {
        std::cerr << "please add kmer first build graph" << std::endl;
        exit(EXIT_FAILURE);
    }

    std::vector<std::string> mString;
    for (int i = 0; i < text.size() - kmer + 1; i++) {
        std::string r = text.substr(i, kmer);
        mString.push_back(r);
    }

    for (int i = 0; i < mString.size() - 1; i++) {
        addNodeToNode(mString.at(i), mString.at(i + 1));
    }

}

void DynamicGraph::addNodeToNode(std::string fromNodeName, std::string toNodeName) {
    DynamicGraph::vRawTextFirst.push_back(fromNodeName);
    DynamicGraph::vRawTextSecond.push_back(toNodeName);
}

std::string DynamicGraph::reverseString(std::string text) {
    std::reverse((text).begin(), (text).end());
    return text;
}

void DynamicGraph::setKmer(int kmer) {
    DynamicGraph::kmer = kmer;
}

int DynamicGraph::getKmer() {
    return DynamicGraph::kmer;
}

GraphResult DynamicGraph::getGraphResult(std::string begin,std::string end) {
//    std::cout << "start" << std::endl;
    GraphResult ResultGR;

    Graph graph1;
    int countread = 0;
    for (int i=0;i<vRawTextFirst.size();i++) {
        graph1.addNodeToNode(vRawTextFirst.at(i),vRawTextSecond.at(i));
        countread++;
    }

//    std::cout << "countread : " << countread << std::endl;
//    graph1.showAllNode();

//    std::cout << "add node to node" << std::endl;
//    std::cout << "begin : " << begin << std::endl;
//    std::cout << "end : " << end << std::endl;
//    graph1.showAllNode();

    GraphResult gr1 = graph1.findCustom(begin,end);

    ResultGR.setMaxLeft(gr1.getMaxLeft());
    ResultGR.setMaxConcordant(gr1.getMaxConcordant());

    Graph graph2;
    for (int i=0;i<vRawTextFirst.size();i++) {
        graph2.addNodeToNode(reverseString(vRawTextSecond.at(i)),reverseString(vRawTextFirst.at(i)));
    }

    GraphResult gr2 = graph2.findCustom(reverseString(end),reverseString(begin));
    ResultGR.setMaxRight(reverseString(gr2.getMaxLeft()));
    if (ResultGR.getMaxConcordant().size()<gr2.getMaxConcordant().size()) {
        ResultGR.setMaxConcordant(reverseString(gr2.getMaxConcordant()));
    }

    return ResultGR;
}

GraphResult DynamicGraph::ReverseStringGraphResult(GraphResult *gr) {
    GraphResult outGR;
//    std::cout << " > " << gr->getMaxLeft() << std::endl;
    outGR.setMaxLeft(reverseString(gr->getMaxLeft()));
//    std::cout << " > " << reverseString(gr->getMaxLeft()) << std::endl;
    outGR.setMaxRight(reverseString(gr->getMaxRight()));
    outGR.setMaxConcordant(reverseString(gr->getMaxConcordant()));
    return outGR;
}
