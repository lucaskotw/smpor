#ifndef BFS_H
#define BFS_H


#include "config.h"
#include "graph.h"
#include <vector>



#define VTX_NOT_CONNECTED -1
#define DISCONNECTED_OFFSET 10



int bfs(Graph::Graph& g, int gSize, VtxType s, std::vector<WgtType>& dist);
int bfs_create_clusters_graph(Graph::Graph& g, int gSize, VtxType s,\
    std::vector< WgtType >& radii,\
    std::vector<int>& clusters, int nCluster, Graph::Graph& cg);


#endif