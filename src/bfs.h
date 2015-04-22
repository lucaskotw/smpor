#ifndef BFS_H
#define BFS_H


#include "config.h"
#include "graph.h"
#include "pgraph.h"
#include <vector>
#include <queue>


#define VTX_NOT_CONNECTED -1



int bfs(Graph::Graph& g, int gSize, VtxType s, std::vector<WgtType>& dist);
int bfs_pg(PGraph::PGraph& pg, int pgSize, VtxType s, std::vector<WgtType>& dist);


#endif