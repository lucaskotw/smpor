#include "sm.h"


static
double stress(DenseMat & distMat, DenseMat& coord)
{
    using namespace std;
    double stress = 0.0;
    int nNodes = coord.rows();

    // other node not central node
    for (int i=0; i<nNodes; ++i)
    {
        for (int j=i+1; j<nNodes; ++j)
        {
            if (i != j)
            {
                stress += pow(distMat(i, j), -2) * \
                pow( ((coord.row(i) - coord.row(j)).norm() - distMat(i, j)), 2 );
            }
        }
    }

    return stress;
}


static
double stress_with_radial(DenseMat & distMat, DenseMat& coord, double coeff)
{
    using namespace std;
    double stress = 0.0;
    int nNodes = coord.rows();
    // radial part
    for (int i=1; i<nNodes; ++i)
    {
        stress += coeff * pow(distMat(0, i), -2) * \
            pow( ((coord.row(0) - coord.row(i)).norm() - distMat(0, i)), 2 );
    }
    // other node not central node
    for (int i=1; i<nNodes; ++i)
    {
        for (int j=i+1; j<nNodes; ++j)
        {
            if (i != j)
            {
                stress += (1-coeff) * pow(distMat(i, j), -2) * \
                pow( ((coord.row(i) - coord.row(j)).norm() - distMat(i, j)), 2 );
            }
        }
    }

    return stress;
}


int stress_majorization(int graphSize,
                        DenseMat& distMat,
                        std::vector< std::vector<CoordType> >& coord)
{
    // Steps
    // 1. Create Corresponding Coordinates Matrix
    // 2. Create weight laplacian matrix
    // 3. Majorization
    //     * Pick the first node as reference point
    // 4. map Eigen matrix to coord
    using namespace std;
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
    double stress_ratio = 0.0;
    double epsl = SM_THRESHOLD;
    double pre_stress;
    double aft_stress;
    DenseMat i_lap(w_lap.rows(), w_lap.cols());
    DenseMat p_sol;
    DenseMat iter_coord;
    while (true)
    {
        pre_stress = stress(distMat, coord_mat);
        // cout << "previous stress=" << pre_stress << endl;
        // create i_lap
        i_lap.fill(0);
        for (int r=0; r<i_lap.rows(); ++r)
        {
            for (int c=0; c<i_lap.cols(); ++c)
            {
                if (r != c)
                {
                    i_lap(r, c) = w_lap(r, c) * distMat(r, c) / (coord_mat.row(r) - coord_mat.row(c)).norm();
                    i_lap(r, r) -= i_lap(r, c);
                }
            }
        }

        // solving part
        p_sol = i_lap * coord_mat;
        iter_coord = w_lap.block(1, 1, i_lap.rows()-1, i_lap.cols()-1)\
            .ldlt().solve(p_sol.block(1, 0, coord_mat.rows()-1, coord_mat.cols()));
        coord_mat.fill(0);  // refill coord_mat for next iteration
        coord_mat.block(1, 0, coord_mat.rows()-1, coord_mat.cols()) = iter_coord;
        aft_stress = stress(distMat, coord_mat);
        // cout << "after stress=" << aft_stress << endl;

        // decide whether to stop stress majorization
        stress_ratio = (pre_stress-aft_stress)/pre_stress;
        if ( (stress_ratio < epsl) || isnan(stress_ratio) ) break;

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


int stress_majorization_radial_refinement(
    int graphSize,
    DenseMat& distMat,
    double coeff,   // the coeff of linear combs of orig stress and constriants
    std::vector< std::vector<CoordType> >& coord)
{
    // Steps
    // 1. Create weight laplacian matrix
    // 2. Majorization
    //     * Pick the first node (center) as reference point
    // 3. map Eigen matrix to coord
    using namespace std;
    int g_size = graphSize;
    int coord_size = coord.size();

    
    DenseMat coord_mat(coord_size, 2);
    int coord_mat_cols = coord_mat.cols();
    int coord_mat_rows = coord_mat.rows();
    for (int c=0; c<coord_mat_cols; c++)
    {
        for (int r=0; r<coord_mat_rows; ++r)
        {
            coord_mat(r, c) = coord.at(r).at(c);
        }
    }
    cout << "init coord" << endl;
    cout << coord_mat << endl;



    // Step 1
    DenseMat w_lap(coord_size, coord_size);
    w_lap.fill(0);
    int w_lap_size = w_lap.rows();
    // deal the radius part
    for (int i=1; i<w_lap_size; ++i) {
        w_lap(0, i) = -1 * coeff * pow(distMat(0, i), -2);
        w_lap(i, 0) = -1 * coeff * pow(distMat(i, 0), -2);

        // i == j part
        w_lap(i, i) -= w_lap(0, i);
    }

    // deal the i != j other part
    for (int i=1; i<w_lap_size; ++i) {

        // i != j
        for (int j=1; j<w_lap_size; ++j)
        {
            if (i != j)
            {
                w_lap(i, j) = -1 * (1-coeff) * pow(distMat(i, j), -2);
                w_lap(i, i) -= w_lap(i, j); // i == j part
            } 

        }

    }


    // cout << "distance matrix" << endl;
    // cout << distMat << endl;
    // cout << "w lap" << endl;
    // cout << w_lap << endl;


    // Step 3
    double stress_ratio = 0.0;
    double dist_norm;
    double epsl = SM_THRESHOLD;
    double pre_stress;
    double aft_stress;
    DenseMat i_lap(w_lap.rows(), w_lap.cols());
    DenseMat p_sol;
    DenseMat iter_coord;
    int i = 0;
    while (true)
    {
        pre_stress = stress_with_radial(distMat, coord_mat, coeff);
        if (i == 0) cout << "previous stress=" << pre_stress << endl;
        // create i_lap
        i_lap.fill(0);
        for (int r=0; r<i_lap.rows(); ++r)
        {
            for (int c=0; c<i_lap.cols(); ++c)
            {
                if (r != c)
                {
                    dist_norm = (coord_mat.row(r) - coord_mat.row(c)).norm();
                    if (dist_norm==0) dist_norm = 1; // [Warning] incase norm = 0;
                    i_lap(r, c) = w_lap(r, c) * distMat(r, c) / dist_norm;
                    if (isinf(i_lap(r, c))) cout << "inf coord norm=" << (coord_mat.row(r) - coord_mat.row(c)).norm() << endl;
                    if (isinf(i_lap(r, c))) cout << "coord " << r << "=" << coord_mat.row(r) << endl;
                    if (isinf(i_lap(r, c))) cout << "coord " << c << "=" << coord_mat.row(c) << endl;
                    i_lap(r, r) -= i_lap(r, c);
                }
            }
        }

        // solving part
        p_sol = i_lap * coord_mat;
        iter_coord = w_lap.block(1, 1, i_lap.rows()-1, i_lap.cols()-1)\
            .ldlt().solve(p_sol.block(1, 0, coord_mat.rows()-1, coord_mat.cols()));
        if (i == 0) cout << "init coord_mat" << endl << coord_mat << endl;
        coord_mat.fill(0);  // refill coord_mat for next iteration
        coord_mat.block(1, 0, coord_mat.rows()-1, coord_mat.cols()) = iter_coord;
        aft_stress = stress_with_radial(distMat, coord_mat, coeff);
        if (i == 0) cout << "after stress=" << aft_stress << endl;
        if (isnan(aft_stress)) cout << "w_lap" << endl << w_lap << endl;
        if (isnan(aft_stress)) cout << "distMat" << endl << distMat << endl;
        if (isnan(aft_stress)) cout << "i_lap" << endl << i_lap << endl;
        if (isnan(aft_stress)) cout << "p_sol" << endl << p_sol << endl;
        if (isnan(aft_stress)) cout << "iter_coord" << endl << iter_coord << endl;
        if (isnan(aft_stress)) cout << "coord_mat" << endl << coord_mat << endl;

        // decide whether to stop stress majorization
        stress_ratio = (pre_stress-aft_stress)/pre_stress;
        if ( (stress_ratio < epsl) || isnan(stress_ratio) ) break;
        ++i;
    }

    // Step 3
    for (int c=0; c<coord_mat.cols(); c++)
    {
        for (int r=0; r<coord_mat.rows(); ++r)
        {
        
            coord.at(r).at(c) = coord_mat(r, c);
        }
    }


    return SUCCESS_SM;
}