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




/****************
 * Main Process *
 ****************/
int smpor(Graph::Graph& g, int graphSize, DenseMat& distMat,\
    std::vector< std::vector<CoordType> >& coord,\
    std::vector< std::vector<CoordType> >& center_coord,\
    std::vector< WgtType >& radii,\
    std::vector<int>& clusters, int nCluster);

#endif