#ifndef PATITION_H
#define PATITION_H

extern "C" {
    #include <stdlib.h> // for srand()
}

#include <vector>
#include <algorithm> // for random_shuffle
#include "graph.h"
#include "cgraph.h"


#define THRESHOLD_VTX     15
#define THRESHOLD_RATIO   0.8
#define MAX_WGT_INIT_VAL  -10000000
#define MAX_GAIN_INIT_VAL -10000000





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


/********************
 * Coarsening Phase *
 ********************/
int heavy_edge_matching(CGraph& pcg, CGraph& cg);
bool is_coarsening_terminate(CGraph::CGraph& pcg, CGraph::CGraph& cg, int partNum);
int coarsening_phase (Graph::Graph& g, std::vector<CGraph::CGraph>& cgs, int partNum);

/******************************
 * Initial Partitioning Phase *
 ******************************/
bool is_partition_num_exceed(std::vector<int>& partition, int partID, int partNum);
bool is_adj_lock(std::vector<VtxType>& part_order, \
                std::vector<bool>& lock, \
                VtxType adjVtx);
void compute_gain(CGraph::CGraph& coarsetCG, \
    std::vector<VtxType>& part_order, \
    std::vector<WgtType>& gain, \
    std::vector<bool>& lock);
void find_max_gain(CGraph::CGraph& coarsetCG, \
                    std::vector<WgtType>& gain, \
                    std::vector<bool>& lock, \
                    WgtType& d, \
                    VtxType& a, \
                    VtxType& b);
void find_k_max_diff(std::vector<WgtType>& diff, \
                    WgtType& diff_max, \
                    int& k);
void exchange_av_with_bv(std::vector<VtxType>& partOrder, \
                        std::vector<VtxType>& av, \
                        std::vector<VtxType>& bv);
void kl_bisection(CGraph::CGraph& coarsetCG, \
    std::vector<int>& coarsetPartition, \
    int partID, int nextPartID);
int partitioning_phase(CGraph::CGraph& coarsetCG, \
    std::vector<int>& coarsetPartition, \
    int partNum);


/**********************
 * Uncoarsening Phase *
 **********************/
void project_partition(CGraph::CGraph& prevCG, \
                       CGraph::CGraph& currCG, \
                       std::vector<int>& prev_partition, \
                       std::vector<int>& curr_partition);
int uncoarsening_phase(std::vector<CGraph::CGraph>& cgs, \
    std::vector<int>& partition, \
    std::vector<int>& coarset_partition);


/********************************
 * Overall Partitioning Process *
 ********************************/
int partition_graph (Graph::Graph& g, std::vector<int>& partition, int partNum);


#endif