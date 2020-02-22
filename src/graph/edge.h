//
// Created by Sakkayaphab Piwluang on 22/2/20.
//

#ifndef BOLT_EDGE_H
#define BOLT_EDGE_H

#include <string>
#include <vector>

class Edge {
    int edgeto = -1;
    bool marked = false;
public:
    Edge();
    void setEdgeTo(int pos);
    int getEdgeTo();
    void setMarked(bool marked);
    bool getMarked();
};

#endif //BOLT_EDGE_H
