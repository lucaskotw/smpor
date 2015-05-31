#ifndef LAP_H
#define LAP_H

#include "config.h"
#include "eigenmat.h"
#include "graph.h"
#include "bfs.h"


double inv_norm(DenseMat& coord, int iR, int jR);
void w_lap_normal(int graphSize, DenseMat& dist, DenseMat& lap);
void iter_lap_normal(int graphSize, DenseMat& dist,\
        DenseMat& lap, DenseMat& coord);

#endif