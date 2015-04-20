#ifndef BFS_H
#define BFS_H


#include "config.h"
#include "graph.h"
#include <vector>
#include <queue>


#define VTX_NOT_CONNECTED -1



int bfs(Graph::Graph& g, int gSize, VtxType s, std::vector<WgtType>& dist);


#endif