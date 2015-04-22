#ifndef LAP_H
#define LAP_H

#include "config.h"
#include "eigenmat.h"
#include "graph.h"
#include "pgraph.h"
#include "bfs.h"


double inv_norm(DenseMat& coord, int iR, int jR);
void w_lap_normal(Graph::Graph& g, DenseMat& dist, DenseMat& lap);
void iter_lap_normal(Graph::Graph& g, DenseMat& dist,\
        DenseMat& lap, DenseMat& coord);
void w_lap_with_vertex_radius(PGraph::PGraph& pg, DenseMat& dist, DenseMat& lap);
void iter_lap_with_vertex_radius(PGraph::PGraph& pg, DenseMat& dist,\
        DenseMat& lap, DenseMat& coord);

#endif