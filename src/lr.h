#ifndef LR_H
#define LR_H


#include "config.h"
#include "graph.h"
#include "eigenmat.h"


#include <vector>


int layout_refinement(Graph::Graph& g,\
    DenseMat& distMat,\
    std::vector<int>& clusters, int nCluster,\
    std::vector< std::vector<CoordType> >& center_coord,\
    std::vector< WgtType >& radii,\
    std::vector< std::vector<CoordType> >& coord);

#endif