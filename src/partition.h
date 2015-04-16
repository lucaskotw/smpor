#ifndef PATITION_H
#define PATITION_H

extern "C" {
    #include <stdlib.h> // for srand()
}

#include <vector>
#include <algorithm> // for random_shuffle
#include "graph.h"
#include "cgraph.h"


#define THRESHOLD_VTX    15
#define THRESHOLD_RATIO  0.8
#define MAX_WGT_INIT_VAL -1


#define SUCCESS_MATCHING   0
#define SUCCESS_COARSENING 0
#define SUCCESS_PARTITION  0


inline void select_vertex_with_heaviest_edge(std::vector<VtxType>& adj,\
                                            std::vector<WgtType>& adj_wgts,\
                                            VtxType& max_wgt_vtx)
{
    double max_wgt = MAX_WGT_INIT_VAL;
    max_wgt_vtx = 0;
    for (int i=0; i<adj.size(); ++i)
    {
        if (adj_wgts.at(i) > max_wgt)
        {
            max_wgt_vtx = adj.at(i);
            max_wgt = adj_wgts.at(i);
        }
    }

}

int heavy_edge_matching(CGraph& pcg, CGraph& cg);
bool is_coarsening_terminate(CGraph::CGraph& pcg, CGraph::CGraph& cg, int partNum);
int coarsening_phase (Graph::Graph& g, std::vector<CGraph::CGraph>& cgs, int partNum);
int partition_graph (Graph::Graph& g, std::vector<int>& partition, int partNum);


#endif