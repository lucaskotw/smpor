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
void find_p_center_radius(Graph::Graph& g, std::vector<int>& partition,\
                   int partID, VtxType& pCenter, WgtType& pRadius);
void add_p_graph_edges(PGraph::PGraph& pg, Graph::Graph& g);
int create_pgraph(PGraph::PGraph& pg, Graph::Graph& g,\
    std::vector<int>& partition, int partNum);

/****************
 * Main Process *
 ****************/
int smpor(Graph::Graph& g,\
          std::vector< std::vector<CoordType> >& coord,\
          std::vector<int>& partition,\
          int partNum);

#endif