#ifndef SM_H
#define SM_H

#include "config.h"
#include "eigenmat.h"
#include "graph.h"
#include "pgraph.h"
#include "lap.h"


#define SM_THRESHOLD 0.0001


// double stress(DenseMat & distMat, std::vector< std::vector<CoordType> >& coord);
int stress_majorization(int graphSize,\
                        DenseMat& distMat,\
                        std::vector< std::vector<CoordType> >& coord);
int stress_majorization_with_pg(PGraph::PGraph & pg,\
                                DenseMat& distMat,\
                                DenseMat& wLap,\
                                std::vector< std::vector<CoordType> >& pg_coord);

#endif