#include "distance.h"
#include "bfs.h"
#include "config.h"

/******************************************************************************
 *                       Stress Majorization                                  *
 *                    of small graph components                               *
 ******************************************************************************/
void distance_matrix(Graph::Graph& g, DenseMat& distMat)
{
    
    int g_size = g.get_num_vtxs();
    std::vector<WgtType> dist(g_size);

    for (int v=0; v<g_size; ++v)
    {
        bfs(g, g_size, v, dist);
        distMat.row(v) = MapVec(&dist[0], dist.size());
    }


}