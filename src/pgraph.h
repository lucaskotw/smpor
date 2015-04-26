#ifndef PGRAPH_H
#define PGRAPH_H


#include <vector>
#include <iostream>

#include "config.h"

#define PORT_BOUNDARY_OFFSET 1


struct p_node
{
                     int nEdges;  // number of edges
                 VtxType pCenter; // corresponding center vertex in given graph
                 WgtType pRadius; // radius of the partition
    std::vector<VtxType> edges;   // each edges[0, ..., k] represent the vertex
                                  // connects to this one
    std::vector<WgtType> pWgts;   // prefered edges weights corresponds to
                                  // edges[0, ..., k]
};

/*
 * ID of vtxs will be the partition number
 */
class PGraph
{
    private:
        std::vector<p_node> vtxs;
    public:
        PGraph();
        ~PGraph();
        void add_node(VtxType pCenter, WgtType pRadius);
        void add_edge(VtxType u, VtxType v, WgtType pWgt);
        VtxType get_center_id(VtxType pID);
        WgtType get_radius(VtxType pID);
        std::vector<VtxType> adj(VtxType s);
        std::vector<WgtType> adj_wgts(VtxType s);
        int get_num_vtxs();
        void print_graph();
};



#endif