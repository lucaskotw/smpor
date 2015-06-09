#include "sm.h"


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
    // 1. Create Corresponding Coordinates Matrix
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
            coord_mat(r, c) = coord.at(r).at(c);
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

        if (aft_stress == 0) { // means that the stress after init is done
            break;
        }

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
        if (isnan(stress_ratio)) {
            std::cout << "graph size = " << g_size << std::endl;
            std::cout << "stress ratio nan show up" << std::endl;
            std::cout << "distance matrix" << std::endl;
            std::cout << distMat << std::endl;
            std::cout << "coordinates" << std::endl;
            std::cout << coord_mat << std::endl;
            std::cout << "previous stress=" << pre_stress << ", after stress=" << aft_stress << std::endl;
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


int stress_majorization_radial_refinement(int graphSize,\
                        DenseMat& distMat,\
                        std::vector< std::vector<CoordType> >& coord)
{
    // Steps
    // 1. Create Corresponding Coordinates Matrix
    // 1. Create weight laplacian matrix
    // 2. Majorization
    //     * Pick the first node (center) as reference point
    // 3. map Eigen matrix to coord

    int g_size = graphSize;
    // Step 1
    DenseMat coord_mat(g_size, 2);
    for (int c=0; c<coord_mat.cols(); c++)
    {
        for (int r=0; r<coord_mat.rows(); ++r)
        {
            coord_mat(r, c) = coord.at(r).at(c);
        }
    }

    // std::cout << "coord_mat" << std::endl;
    // std::cout << coord_mat << std::endl;


    // Step 2
    DenseMat w_lap(g_size, g_size);
    w_lap_normal(g_size, distMat, w_lap);
    // std::cout << "first w lap" << std::endl;
    // std::cout << w_lap << std::endl;


    // Step 3
    // First Iteration
    double stress_ratio = 0.0;
    double pre_stress = stress(distMat, coord_mat);
    DenseMat iter_coord(graphSize, 2);
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

    // std::cout << "first step coord mat" << std::endl;
    // std::cout << coord_mat << std::endl;
    
    while ( true ) {
        // std::cout << "after stress=" << aft_stress << std::endl;
        if (aft_stress == 0) { // means that the stress after init is done
            break;
        }

        iter_lap_normal(g_size, distMat, i_lap, coord_mat);
        // if (i == 1)
        // {
        //     std::cout << "first step iter lap" << std::endl;
        //     std::cout << i_lap << std::endl;
        // }

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
        if (isnan(stress_ratio)) {
            std::cout << "graph size = " << g_size << std::endl;
            std::cout << "stress ratio nan show up" << std::endl;
            std::cout << "distance matrix" << std::endl;
            std::cout << distMat << std::endl;
            std::cout << "coordinates" << std::endl;
            std::cout << coord_mat << std::endl;
            std::cout << "previous stress=" << pre_stress << ", after stress=" << aft_stress << std::endl;
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