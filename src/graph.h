/*
 * Declare simple graph structure
 * Assume the graph is
 *   1. simple
 *   2. undirected
 */
#ifndef GRAPH_H
#define GRAPH_H


#include <vector>
#include <iostream>


#include "config.h"


struct v_dat
{
    int nEdges;
    std::vector<VtxType> edges; // each edges[0, ..., k] represent the vertex
                                // connects to this one
    std::vector<WgtType> pWgts; // prefered edges weights corresponds to edges[0, ..., k]
};

class Graph
{
    private:
        std::vector<v_dat> vtxs;
    public:
        Graph(VtxType nNodes);
        ~Graph();
        void add_edge(VtxType u, VtxType v, WgtType wgt);
        std::vector<VtxType> adj(VtxType s);
        std::vector<WgtType> adj_wgts(VtxType s);
        int get_num_vtxs();
        void print_graph();
};

#endif