#include "sm.h"


double coord_norm(std::vector<CoordType>& x1, std::vector<CoordType>& x2)
{
    double x_val = std::pow( x1.at(0) - x2.at(0), 2 );
    double y_val = std::pow( x1.at(1) - x2.at(1), 2 );
    return std::sqrt(x_val+y_val);
}

double stress(DenseMat & distMat, DenseMat& coord)
{
    double stress = 0.0;
    int nNodes = coord.rows();
    for (int i=0; i<nNodes; ++i) {
        for (int j=i+1; j<nNodes; ++j) {
            if (i != j) {
                stress += std::pow(distMat(i, j), -2) * \
                std::pow( ((coord.row(i) - coord.row(j)).norm() - distMat(i, j)), 2 );
            }
        }
    }

    return stress;
}

int stress_majorization(Graph::Graph & g,\
                        DenseMat& distMat,\
                        DenseMat& wLap,\
                        std::vector< std::vector<CoordType> >& coord,
                        int pCenter)
{
    // Transform the pg_coord to matrix format
    int g_size = g.get_num_vtxs();
    DenseMat coord_mat(g_size, 2);
    for (int r=0; r<g_size; ++r)
    {
        for (int c=0; c<2; c++)
        {
            coord_mat(r, c) = coord.at(r).at(c);
        }
    }

    // switch distance and wlap with the pCenter ones
    distMat.row(0).swap(distMat.row(pCenter));
    distMat.col(0).swap(distMat.col(pCenter));
    wLap.row(0).swap(wLap.row(pCenter));
    wLap.col(0).swap(wLap.col(pCenter));

    // Create weight laplacian matrix
    double stress_ratio = 0.0;
    double pre_stress = stress(distMat, coord_mat);
    DenseMat iter_coord(g_size-1, 2);
    iter_coord = coord_mat.block(1, 0, coord_mat.rows()-1, coord_mat.cols());
    iter_coord = wLap.block(1, 1, wLap.rows()-1, wLap.cols()-1).ldlt().solve(iter_coord);
    coord_mat.fill(0);
    coord_mat.block(1, 0, coord_mat.rows()-1, coord_mat.cols()) = iter_coord;
    double aft_stress = stress(distMat, coord_mat);
    stress_ratio = (pre_stress-aft_stress)/pre_stress;
    std::cout << "round 0 previous stress = " << pre_stress << std::endl;
    std::cout << "round 0 after stress = " << aft_stress << std::endl;
    std::cout << "ratio = " << stress_ratio << std::endl;


    // First Iteration
    double epsl = SM_THRESHOLD;
    
    DenseMat i_lap(wLap.rows(), wLap.cols());
    DenseMat p_sol(wLap.rows(), wLap.cols());
    int i = 1;
    
    while ( true ) {
        iter_lap_normal(g, distMat, i_lap, coord_mat);

        p_sol = i_lap * coord_mat;
        iter_coord = wLap.block(1, 1, i_lap.rows()-1, i_lap.cols()-1)\
            .ldlt().solve(p_sol.block(1, 0, coord_mat.rows()-1, coord_mat.cols()));

        // assign preStress based on current iteration
        pre_stress = aft_stress;
        coord_mat.fill(0);
        coord_mat.block(1, 0, coord_mat.rows()-1, coord_mat.cols()) = iter_coord;
        // assign aftStress after coord assign
        aft_stress = stress(distMat, coord_mat);
        
        stress_ratio = (pre_stress-aft_stress)/pre_stress;
        std::cout << "round " << i << " stress = " << aft_stress << std::endl;
        std::cout << "ratio = " << (pre_stress-aft_stress)/pre_stress << std::endl;

        if (stress_ratio < epsl) {
            break;
        }
        ++i;
    }

    // switch distance and wlap back
    // also switch coord_mat (since the operation is after swaping)
    distMat.row(0).swap(distMat.row(pCenter));
    distMat.col(0).swap(distMat.col(pCenter));
    wLap.row(0).swap(wLap.row(pCenter));
    wLap.col(0).swap(wLap.col(pCenter));
    coord_mat.row(0).swap(coord_mat.row(pCenter));

    // transform coord_mat back to pg_coord
    for (int r=0; r<g_size; ++r)
    {
        for (int c=0; c<2; c++)
        {
            coord.at(r).at(c) = coord_mat(r, c);
        }
    }

    return SUCCESS_SM;
}


int stress_majorization_with_pg(PGraph::PGraph & pg,\
                                DenseMat& distMat,\
                                DenseMat& wLap,\
                                std::vector< std::vector<CoordType> >& pg_coord)
{
    // Transform the pg_coord to matrix format
    int pg_size = pg.get_num_vtxs();
    DenseMat coord_mat(pg_size, 2);
    for (int r=0; r<pg_size; ++r)
    {
        for (int c=0; c<2; c++)
        {
            coord_mat(r, c) = pg_coord.at(r).at(c);
        }
    }

    // Create weight laplacian matrix
    double stress_ratio = 0.0;
    double pre_stress = stress(distMat, coord_mat);
    DenseMat iter_coord(pg_size-1, 2);
    iter_coord = coord_mat.block(1, 0, coord_mat.rows()-1, coord_mat.cols());
    iter_coord = wLap.block(1, 1, wLap.rows()-1, wLap.cols()-1).ldlt().solve(iter_coord);
    coord_mat.fill(0);
    coord_mat.block(1, 0, coord_mat.rows()-1, coord_mat.cols()) = iter_coord;
    double aft_stress = stress(distMat, coord_mat);
    stress_ratio = (pre_stress-aft_stress)/pre_stress;
    std::cout << "round 0 previous stress = " << pre_stress << std::endl;
    std::cout << "round 0 after stress = " << aft_stress << std::endl;
    std::cout << "ratio = " << stress_ratio << std::endl;


    // First Iteration
    double epsl = SM_THRESHOLD;
    
    DenseMat i_lap(wLap.rows(), wLap.cols());
    DenseMat p_sol(wLap.rows(), wLap.cols());
    int i = 1;
    
    while ( true ) {
        iter_lap_with_vertex_radius(pg, distMat, i_lap, coord_mat);

        p_sol = i_lap * coord_mat;
        iter_coord = wLap.block(1, 1, i_lap.rows()-1, i_lap.cols()-1)\
            .ldlt().solve(p_sol.block(1, 0, coord_mat.rows()-1, coord_mat.cols()));

        // assign preStress based on current iteration
        pre_stress = aft_stress;
        coord_mat.fill(0);
        coord_mat.block(1, 0, coord_mat.rows()-1, coord_mat.cols()) = iter_coord;
        // assign aftStress after coord assign
        aft_stress = stress(distMat, coord_mat);
        
        stress_ratio = (pre_stress-aft_stress)/pre_stress;
        std::cout << "round " << i << " stress = " << aft_stress << std::endl;
        std::cout << "ratio = " << (pre_stress-aft_stress)/pre_stress << std::endl;

        if (stress_ratio < epsl) {
            break;
        }
        ++i;
    }

    // transform coord_mat back to pg_coord
    for (int r=0; r<pg_size; ++r)
    {
        for (int c=0; c<2; c++)
        {
            pg_coord.at(r).at(c) = coord_mat(r, c);
        }
    }

    return SUCCESS_SM_PG;

}