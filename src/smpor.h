#ifndef SMPOR_H
#define SMPOR_H

#include "config.h"
#include "eigenmat.h"
#include "graph.h"
#include "pgraph.h"


#include "bfs.h"
#include "lap.h"
#include "sm.h"


#include <vector>
#include <algorithm>  // std::find, std::max_element


#define INIT_COORD       0
#define WEIGHT_PARAMETER 2



/******************
 * Create P-Graph *
 ******************/
void find_p_center_radius(Graph::Graph& g, std::vector<Graph>& sg_vec,\
                          std::vector<int>& clusters,\
                          int partID, VtxType& pCenter, WgtType& pRadius);
void add_p_graph_edges(PGraph::PGraph& pg, Graph::Graph& g);
int create_small_graph_list(PGraph::PGraph& pg, std::vector<Graph>& sg_vec,\
    Graph::Graph& g, std::vector<int>& clusters, int nCluster);

void match_partition_coord(std::vector<int>& clusters_vtxs,\
    std::vector< std::vector<CoordType> >& partition_coords,\
    std::vector< std::vector<CoordType> >& pg_coord,\
    std::vector< std::vector<CoordType> >& coord,\
    int partID,\
    int pCenter);


/****************
 * Main Process *
 ****************/
int smpor(int graphSize, DenseMat& distMat,\
    std::vector< std::vector<CoordType> >& coord,\
    std::vector< std::vector<CoordType> >& center_coord,\
    std::vector< WgtType >& radius,\
    std::vector<int>& clusters, int nCluster);

#endif