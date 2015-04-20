#ifndef SMPOR_H
#define SMPOR_H

#include "config.h"
#include "eigenmat.h"
#include "graph.h"
#include "pgraph.h"


#include "bfs.h"


#include <vector>
#include <algorithm>  // std::find, std::max_element



/******************
 * Create P-Graph *
 ******************/
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