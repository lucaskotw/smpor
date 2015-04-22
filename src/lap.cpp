#include "lap.h"



double inv_norm(DenseMat& coord, int iR, int jR)
{
    DenseVec diff(coord.cols());
    diff = (coord.row(iR) - coord.row(jR));
    return 1/diff.norm();
}



void w_lap_normal(Graph::Graph& g, DenseMat& dist, DenseMat& lap)
{
    lap.fill(0);
    DenseVec dist_row;
    for (int i=0; i<g.get_num_vtxs(); ++i) {

        // if node j is node i's neighbor
        dist_row = dist.row(i);
        lap.row(i) = -1 * dist_row.array().pow(-2);

        lap(i, i) = 0;
        for (int j=0; j<dist_row.size(); ++j) {
            if (i != j) {
                lap(i, i) += pow(dist_row(j), -2);
            }
        }

    }

}


void iter_lap_normal(Graph::Graph& g, DenseMat& dist,\
        DenseMat& lap, DenseMat& coord)
{
    int n_vtxs = g.get_num_vtxs();
    lap.fill(0);
    DenseVec dist_row(n_vtxs);
    for (int i=0; i<n_vtxs; ++i) {

        // if node j is node i's neighbor
        dist_row = dist.row(i);

        for (int j=0; j<dist_row.size(); ++j) {
            if (i != j) {
                lap(i, j) = std::pow(dist_row(j), -2); // w_{ij} = d_{ij}^-2
                lap(i, j) = -1 * lap(i, j) * dist_row(j) * inv_norm(coord, i, j);

                lap(i, i) -= lap(i, j); // added to diagonal
            }
        }
    }

}


void w_lap_with_vertex_radius(PGraph::PGraph& pg, DenseMat& dist, DenseMat& lap)
{
    lap.fill(0);
    DenseVec distRow;
    for (int i=0; i<pg.get_num_vtxs(); ++i) {

        // if node j is node i's neighbor
        distRow = dist.row(i);
        lap.row(i) = -1 * distRow.array().pow(-2);

        lap(i, i) = 0;
        for (int j=0; j<distRow.size(); ++j) {
            if (i != j) {
                lap(i, i) += pow(distRow(j), -2);
            }
        }

    }

}


void iter_lap_with_vertex_radius(PGraph::PGraph& pg, DenseMat& dist,\
        DenseMat& lap, DenseMat& coord)
{
    int n_vtxs = pg.get_num_vtxs();
    lap.fill(0);
    DenseVec dist_row(n_vtxs);
    for (int i=0; i<n_vtxs; ++i) {

        // if node j is node i's neighbor
        dist_row = dist.row(i);

        for (int j=0; j<dist_row.size(); ++j) {
            if (i != j) {
                lap(i, j) = std::pow(dist_row(j), -2); // w_{ij} = d_{ij}^-2
                lap(i, j) = -1 * lap(i, j) * dist_row(j) * inv_norm(coord, i, j);

                lap(i, i) -= lap(i, j); // added to diagonal
            }
        }
    }

}