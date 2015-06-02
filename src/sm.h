#ifndef SM_H
#define SM_H

#include "config.h"
#include "eigenmat.h"
#include "graph.h"
#include "lap.h"
#include <math.h>


#define SM_THRESHOLD 0.0001


// double stress(DenseMat & distMat, std::vector< std::vector<CoordType> >& coord);
int stress_majorization(int graphSize,\
                        DenseMat& distMat,\
                        std::vector< std::vector<CoordType> >& coord);

int stress_majorization_radial_refinement(int graphSize,\
                        DenseMat& distMat,\
                        std::vector< std::vector<CoordType> >& coord);

#endif