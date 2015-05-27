#include "sm.h"


static
double coord_norm(std::vector<CoordType>& x1, std::vector<CoordType>& x2)
{
    double x_val = std::pow( x1.at(0) - x2.at(0), 2 );
    double y_val = std::pow( x1.at(1) - x2.at(1), 2 );
    return std::sqrt(x_val+y_val);
}


static
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


int stress_majorization(int graphSize,\
                        DenseMat& distMat,\
                        std::vector< std::vector<CoordType> >& coord)
{
    // Steps
    // 1. Initially assign the coordinates (random)
    // 2. Create weight laplacian matrix
    // 3. Majorization
    //     * Pick the first node as reference point
    // 4. map Eigen matrix to coord

    int g_size = graphSize;
    // Step 1
    DenseMat coord_mat(g_size, 2);
    for (int c=0; c<coord_mat.cols(); c++)
    {
        for (int r=0; r<coord_mat.rows(); ++r)
        {
            coord_mat(r, c) = rand()%100/50.0;
        }
    }

    // Step 2
    DenseMat w_lap(g_size, g_size);
    w_lap_normal(g_size, distMat, w_lap);


    // Step 3
    // First Iteration
    double stress_ratio = 0.0;
    double pre_stress = stress(distMat, coord_mat);
    DenseMat iter_coord(graphSize-1, 2);
    iter_coord = coord_mat.block(1, 0, coord_mat.rows()-1, coord_mat.cols());
    iter_coord = w_lap.block(1, 1, w_lap.rows()-1, w_lap.cols()-1).ldlt().solve(iter_coord);
    coord_mat.fill(0);  // refill coord_mat for next iteration
    coord_mat.block(1, 0, coord_mat.rows()-1, coord_mat.cols()) = iter_coord;
    double aft_stress = stress(distMat, coord_mat);
    stress_ratio = (pre_stress-aft_stress)/pre_stress;
    // std::cout << "round 0 previous stress = " << pre_stress << std::endl;
    // std::cout << "round 0 after stress = " << aft_stress << std::endl;
    // std::cout << "ratio = " << stress_ratio << std::endl;


    double epsl = SM_THRESHOLD;
    
    DenseMat i_lap(w_lap.rows(), w_lap.cols());
    DenseMat p_sol;
    int i = 1;
    
    while ( true ) {
        iter_lap_normal(g_size, distMat, i_lap, coord_mat);

        p_sol = i_lap * coord_mat;
        iter_coord = w_lap.block(1, 1, i_lap.rows()-1, i_lap.cols()-1)\
            .ldlt().solve(p_sol.block(1, 0, coord_mat.rows()-1, coord_mat.cols()));

        // assign preStress based on current iteration
        pre_stress = aft_stress;
        coord_mat.fill(0);
        coord_mat.block(1, 0, coord_mat.rows()-1, coord_mat.cols()) = iter_coord;
        // assign aftStress after coord assign
        aft_stress = stress(distMat, coord_mat);
        
        stress_ratio = (pre_stress-aft_stress)/pre_stress;
        // std::cout << "round " << i << " stress = " << aft_stress << std::endl;
        // std::cout << "ratio = " << (pre_stress-aft_stress)/pre_stress << std::endl;

        if (stress_ratio < epsl) {
            break;
        }
        ++i;
    }


    // transform coord_mat back to pg_coord
    for (int c=0; c<coord_mat.cols(); c++)
    {
        for (int r=0; r<coord_mat.rows(); ++r)
        {
        
            coord.at(r).at(c) = coord_mat(r, c);
        }
    }

    return SUCCESS_SM;
}