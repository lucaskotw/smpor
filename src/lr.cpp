#include "lr.h"
#include "sm.h"
#include <math.h>
#include <boost/bind.hpp>


#define PI 3.14159


/******************************************************************************
 *                        Radial Refinement                                   *
 ******************************************************************************/
/* --------------------------------- Step 1 --------------------------------- */
/* calculate_vtxs_radii: Step 1 */
static
void calculate_connection(Graph::Graph& g,\
    std::vector<int>& clusters, int nCluster,\
    std::vector<VtxType>& inter_conn,\
    std::vector<VtxType>& intra_conn)
{
    int n_vtxs = g.get_num_vtxs();
    std::cout << "clusters" << std::endl;
    for (int i=0; i<clusters.size(); ++i)
    {
        std::cout << clusters.at(i) << " ";
    }
    std::cout << std::endl;
    std::vector<VtxType> nbors;
    for (int i=0; i<n_vtxs; ++i)
    {
        nbors = g.adj(i);
        for (int n=0; n<nbors.size(); ++n)
        {
            if (clusters.at(i) == clusters.at( nbors.at(n) ))
            {
                intra_conn.at(i) = intra_conn.at(i) + 1;
            }
            else
            {
                inter_conn.at(i) = inter_conn.at(i) + 1;
            }
        }
    }


}


/* calculate_vtxs_radii: Step 2 */
static
void calculate_primary_radii(std::vector< std::vector<CoordType> >& coord,\
    std::vector<WgtType>& radii)
{
    using namespace std;
    CoordType x;
    CoordType y;
    for (int i=0; i<coord.size(); ++i)
    {
        x = coord.at(i).at(0);
        y = coord.at(i).at(1);
        radii.at(i) = sqrt(pow(x, 2) + pow(y, 2));
    }
}

static
void calculate_rank(Graph::Graph& g,\
    std::vector<int>& clusters, int nCluster,\
    std::vector<VtxType>& inter_conn,\
    std::vector<VtxType>& intra_conn,\
    std::vector<VtxType>& rank)
{
    // Steps
    // 1. Give Initial Rank (less num tends to be inside)
    // 2. Calculate the score for ranking
    // 3. Interchange the rank based on score

    // Step 1
    int n_vtxs = g.get_num_vtxs();
    std::vector<int> n_vtxs_c(nCluster, 0);
    int vtx_c;
    for (int i=0; i<n_vtxs; ++i)
    {
        vtx_c = clusters.at(i);
        rank.at(i) = n_vtxs_c.at(vtx_c);
        n_vtxs_c.at(vtx_c) += 1;
    }

    // Step 2
    std::vector<double> score_vtxs(n_vtxs, 0);
    double score;
    double out_conn;
    double in_conn;
    for (int i=0; i<n_vtxs; ++i)
    {
        out_conn = (double) inter_conn.at(i);
        in_conn = (double) intra_conn.at(i);
        score = (out_conn+1) / in_conn / (out_conn+in_conn);
        score_vtxs.at(i) = score;
    }

    std::cout << "calculated score" << std::endl;
    for (int i=0; i<n_vtxs; ++i)
    {
        std::cout << score_vtxs.at(i) << " ";
    }
    std::cout << std::endl;


    // Step 3
    std::vector< std::vector< std::pair<VtxType, WgtType> > > rank_pair(nCluster,
                                    std::vector< std::pair<VtxType, WgtType> >(0));

    for (int i=0; i<n_vtxs; ++i)
    {
        rank_pair.at(clusters.at(i)).push_back( std::make_pair(i, score_vtxs.at(i)) );
    }

    // sort by pair method

    for (int c=0; c<nCluster; ++c)
    {
        std::sort(rank_pair.at(c).begin(), rank_pair.at(c).end(), 
            boost::bind(&std::pair<VtxType, WgtType>::second, _1) <
            boost::bind(&std::pair<VtxType, WgtType>::second, _2));
    }

    std::cout << "rank pair" << std::endl;
    for (int c=0; c<nCluster; ++c)
    {
        std::cout << "cluster #" << c << std::endl;
        for (int i=0; i<rank_pair.at(c).size(); ++i)
        {
            std::cout << rank_pair.at(c).at(i).first << ", " << rank_pair.at(c).at(i).second << std::endl;
        }
    }

    // combine the rank
    for (int c=0; c<nCluster; ++c)
    {
        for (int i=0; i<rank_pair.at(c).size(); ++i)
        {
            rank.at( rank_pair.at(c).at(i).first ) = i;
        }
    }

    std::cout << "rank" << std::endl;
    for (int i=0; i<rank.size(); ++i)
    {
        std::cout << rank.at(i) << " ";
    }
    std::cout << std::endl;
    

}


static
void sort_radii(std::vector<int>& clusters, int nCluster,\
    std::vector<WgtType>& radii,\
    std::vector< std::vector<WgtType> >& cluster_radii)
{
    // Steps
    // 1. Create vectors of clusters' radii
    // 2. Sort each vector of radii of clusters
    using namespace std;

    // Step 1
    for (int i=0; i<clusters.size(); ++i)
    {
        cluster_radii.at( clusters.at(i) ).push_back( radii.at(i) );
    }

    // Step 2
    for (int c=0; c<nCluster; ++c)
    {
        sort(cluster_radii.at(c).begin(), cluster_radii.at(c).end());
    }

    cout << "clusters' radii" << endl;
    for (int c=0; c<nCluster; ++c)
    {
        for (int i=0; i<cluster_radii.at(c).size(); ++i)
        {
            cout << cluster_radii.at(c).at(i) << " ";
        }
        cout << endl;
    }



}


static
void assign_radii_by_rank(std::vector<VtxType>& rank,\
    std::vector<int>& clusters,\
    std::vector< std::vector<WgtType> >& cluster_radii,\
    std::vector<WgtType>& radii)
{
    using namespace std;
    int vtx_r;
    int vtx_c;
    for (int i=0; i<rank.size(); ++i)
    {
        vtx_r = rank.at(i);
        vtx_c = clusters.at(i);
        radii.at(i) = cluster_radii.at(vtx_c).at(vtx_r);
    }
    cout << "radii" << endl;
    for (int i=0; i<radii.size(); ++i)
    {
        cout << radii.at(i) << " ";
    }
    cout << endl;
    
}


static
void assign_radii(Graph::Graph& g,\
    std::vector<int>& clusters, int nCluster,\
    std::vector<VtxType>& inter_conn,\
    std::vector<VtxType>& intra_conn,\
    std::vector< std::vector<CoordType> >& coord,\
    std::vector<WgtType>& radii)
{
    // Steps
    // 1. Get the non-processing radii: O(N)
    // 2. Calculate radii rank
    // 3. Sort the radii by value and by cluster
    // 4. Assign the radii by rank
    using namespace std;

    // Step 1
    calculate_primary_radii(coord, radii);

    // Step 2
    std::vector<VtxType> rank(g.get_num_vtxs(), 0);
    calculate_rank(g, clusters, nCluster, inter_conn, intra_conn, rank);

    // Step 3
    vector< vector<WgtType> > cluster_radii(nCluster, vector<WgtType>(0));
    sort_radii(clusters, nCluster, radii, cluster_radii);

    // Step 4
    assign_radii_by_rank(rank, clusters, cluster_radii, radii);
}


static
void calculate_vtxs_radii(Graph::Graph& g,\
    std::vector<int>& clusters, int nCluster,\
    std::vector< std::vector<CoordType> >& coord,\
    std::vector<WgtType>& radii)
{
    // Steps
    // 1. collect the number inter- and intra-connection of each vtx: O(E)
    //     <run double edges>
    // 2. assign the radii of vtxs

    // Step 1
    std::vector<VtxType> inter_conn(g.get_num_vtxs(), 0);
    std::vector<VtxType> intra_conn(g.get_num_vtxs(), 0);
    calculate_connection(g, clusters, nCluster, inter_conn, intra_conn);

    // Step 2
    assign_radii(g, clusters, nCluster, inter_conn, intra_conn, coord, radii);


}



/* --------------------------------- Step 2 --------------------------------- */
static
void create_cluster_dist_mat(DenseMat& distMat,\
                            std::vector<VtxType>& cluster_vtxs,\
                            std::vector<WgtType>& radii,\
                            DenseMat& clusterDistMat)
{
    // Steps
    // Shift clusterDistMat 1 id
    // 1. resize cluster distance matrix
    // 2. match the distance from distMat to clusterDistMat
    using namespace std;

    // Step 1
    clusterDistMat.resize(cluster_vtxs.size()+1, cluster_vtxs.size()+1);
    clusterDistMat.fill(0);
    cout << "cluster rows=" << clusterDistMat.rows() << endl;
    cout << "cluster cols=" << clusterDistMat.cols() << endl;

    // Step 2
    int rowIdx;
    int colIdx;
    for (int c=0; c<cluster_vtxs.size(); ++c)
    {
        for (int r=0; r<cluster_vtxs.size(); ++r)
        {
            rowIdx = cluster_vtxs.at(r);
            colIdx = cluster_vtxs.at(c);
            clusterDistMat(r+1, c+1) = 0.5*distMat(rowIdx, colIdx);
        }
        clusterDistMat(0, c+1) = 0.5*radii.at(cluster_vtxs.at(c));
        clusterDistMat(c+1, 0) = 0.5*radii.at(cluster_vtxs.at(c));
    }


    cout << "cluster dist mat" << endl;
    cout << clusterDistMat << endl;

}



static
void match_cluster_coord(std::vector<VtxType>& cluster_vtxs,\
                        std::vector< std::vector<CoordType> >& intra_coord,
                        std::vector< std::vector<CoordType> >& coord)
{
    // Assign intro_coord to coord based on cluster_vtxs
    int vtx_id;
    for (int i=0; i<cluster_vtxs.size(); ++i)
    {
        vtx_id = cluster_vtxs.at(i);
        coord.at(vtx_id).at(0) = intra_coord.at(i).at(0);
        coord.at(vtx_id).at(1) = intra_coord.at(i).at(1);
    }
}



static
void radial_stress_majorization(DenseMat& distMat,\
    std::vector<WgtType>& radii,\
    std::vector<int>& clusters, int nCluster,\
    std::vector< std::vector<CoordType> >& coord)
{
    // Steps
    // Forloop from 0 to nCluster-1
    //   1. Get each clusters' distance matrix
    //   2. run the stress majorization of each cluster with radial refinement

    // [check] whether the coord will deviate while not considering the centers'
    //         coord: yes -> so put the center to the coord
    //                and the stress majorization should not be random
    //                (assigned already)
    using namespace std;
    vector< vector<CoordType> > intra_coord;
    vector<VtxType> cluster_vtxs; // id: vtx id of clustered graph
                                       // val: vtx id in whole graph
    DenseMat cluster_dist_mat;
    vector<CoordType> center(2);
    vector< vector<CoordType> > cluster_center(nCluster, vector<CoordType>(2));
    
    for (int p=0; p<nCluster; ++p)
    {
        // loop step 1
        center.at(0) = 0;
        center.at(1) = 0;
        cluster_vtxs.clear();
        cluster_vtxs.resize(0);
        
        
        for (int i=0; i<clusters.size(); ++i)
        {
            if (clusters.at(i) == p)
            {
                cluster_vtxs.push_back(i);
            }

        }
        intra_coord.clear();
        intra_coord.resize(0, std::vector<CoordType>(2));
        intra_coord.push_back(center);
        for (int i=0; i<clusters.size(); ++i)
        {
            if (clusters.at(i) == p)
            {
                cout << i << endl;
                intra_coord.push_back(coord.at(i));
                center.at(0) += coord.at(i).at(0);
                center.at(1) += coord.at(i).at(1);
            }

        }


        center.at(0) /= cluster_vtxs.size();
        center.at(1) /= cluster_vtxs.size();
        cout << "center: (" << center.at(0) << ", " << center.at(0) << ")" << endl;
        intra_coord.at(0) = center;
        cluster_center.at(p) = center;

        // shift intra coordinates with center coordinates;
        // for (int i=0; i<intra_coord.size(); ++i)
        // {
        //     intra_coord.at(i).at(0) -= center.at(0);
        //     intra_coord.at(i).at(1) -= center.at(1);
        // }


        

        // loop step 2
        create_cluster_dist_mat(distMat, cluster_vtxs, radii, cluster_dist_mat);
        std::cout << "cluster dist mat rows = " << cluster_dist_mat.rows() << std::endl;
        std::cout << "intra coord size = " << intra_coord.size() << std::endl;
        // stress_majorization_radial_refinement(cluster_dist_mat.rows(),\
        //                     cluster_dist_mat,\
        //                     intra_coord);

        

        // loop step 3
        for (int i=0; i<intra_coord.size(); ++i)
        {
            intra_coord.at(i).at(0) += center.at(0);
            intra_coord.at(i).at(1) += center.at(1);
        }
        std::cout << "coord" << std::endl;
        for (std::vector< std::vector<CoordType> >::iterator it1=intra_coord.begin();\
            it1!=intra_coord.end();\
            ++it1)
        {
            for (std::vector<CoordType>::iterator it2=(*it1).begin();\
            it2!=(*it1).end();
            ++it2)
            {
                std::cout << *it2 << " ";
            }
            std::cout << std::endl;
        }
        intra_coord.erase(intra_coord.begin());

        std::cout << "cluster #" << p << std::endl;
        std::cout << "coord" << std::endl;
        for (std::vector< std::vector<CoordType> >::iterator it1=intra_coord.begin();\
            it1!=intra_coord.end();\
            ++it1)
        {
            for (std::vector<CoordType>::iterator it2=(*it1).begin();\
            it2!=(*it1).end();
            ++it2)
            {
                std::cout << *it2 << " ";
            }
            std::cout << std::endl;
        }
        std::cout << intra_coord.size() << std::endl;
        std::cout << "_____________ " << std::endl;
        match_cluster_coord(cluster_vtxs, intra_coord, coord);
       
    }
}



static
void radial_refinement(Graph::Graph& g,\
    DenseMat& distMat,\
    std::vector<int>& clusters, int nCluster,\
    std::vector< std::vector<CoordType> >& coord)
{
    // Steps
    // 1. Calculate the radii of vertices within clutser based on the origin
    // 2. run majorization of new stress model (linear combination of original
    //     stress and radial constraints)

    // Step 1
    std::vector<WgtType> radii(g.get_num_vtxs());
    calculate_vtxs_radii(g, clusters, nCluster, coord, radii);

    // Step 2
    radial_stress_majorization(distMat, radii, clusters, nCluster, coord);
}



/******************************************************************************
 *                        Rotate Refinement                                   *
 ******************************************************************************/
static
void polar_transformation(std::vector<int>& clusters, int nCluster,\
    std::vector< std::vector<CoordType> >& center_coord,\
    std::vector< std::vector<CoordType> >& coord,\
    std::vector<WgtType>& radius_coord,\
    std::vector<WgtType>& theta_coord)
{
    // Steps
    // For each vertex
    //   1. Calculate the radius
    //   2. Calculate the theta
    using namespace std;
    int vtx_c;
    CoordType vtx_x;
    CoordType vtx_y;
    CoordType center_x;
    CoordType center_y;
    for (int i=0; i<coord.size(); ++i)
    {
        vtx_c = clusters.at(i);
        vtx_x = coord.at(i).at(0);
        vtx_y = coord.at(i).at(1);

        center_x = center_coord.at(vtx_c).at(0);
        center_y = center_coord.at(vtx_c).at(1);

        // radius
        radius_coord.at(i) =\
            sqrt(pow( (vtx_x - center_x) , 2) + pow( (vtx_y - center_y) , 2));

        // theta
        theta_coord.at(i)=\
            atan2( (vtx_y - center_y), (vtx_x - center_x) ) * 180.0 / PI;

    }

}


static
void calculate_rotate_degree(Graph::Graph& g,\
    std::vector<int>& clusters, int nCluster,\
    std::vector< std::vector<CoordType> >& center_coord,\
    std::vector< WgtType >& radii,\
    std::vector< std::vector<CoordType> >& coord,\
    std::vector<WgtType>& theta_coord,\
    std::vector<WgtType>& rotate_degree)
{
    // Steps
    // 1. For each vertex
    //     1. Calculate each vertex's inter-connections
    //     2. measure the corresponding ideal points of the vertices
    // 2. add up the difference that need to rotate

    // Step 1
    using namespace std;
    vector< pair<int, WgtType> > theta_diff; // first: cluster id, second: diff
    vector<VtxType> nbors;
    int vtx_c;
    CoordType center_x;
    CoordType center_y;
    CoordType nbor_x;
    CoordType nbor_y;
    WgtType   nbor_theta;
    for (int i=0; i<coord.size(); ++i)
    {
        vtx_c = clusters.at(i);
        center_x = center_coord.at(vtx_c).at(0);
        center_y = center_coord.at(vtx_c).at(1);
        nbors.clear();
        nbors = g.adj(i);
        for (int n=0; n<nbors.size(); ++n)
        {
            if (vtx_c != clusters.at(nbors.at(n)))
            {
                nbor_x = coord.at(nbors.at(n)).at(0);
                nbor_y = coord.at(nbors.at(n)).at(1);
                // calculate ideal theta
                nbor_theta = atan2((nbor_y - center_y), (nbor_x - center_x));
                theta_diff.push_back( make_pair(vtx_c,\
                                                nbor_theta - theta_coord.at(i))
                                    );


            }
        }
    }

    // Step 2
    int c_id;
    vector<int> rotate_cnt(nCluster, 0);
    for (int i=0; i<theta_diff.size(); ++i)
    {
        c_id = theta_diff.at(i).first;
        rotate_cnt.at(c_id) += 1;
        rotate_degree.at(c_id) += theta_diff.at(i).second;
    }

    for (int r=0; r<rotate_cnt.size(); ++r)
    {
        rotate_degree.at(r) /= rotate_cnt.at(r);
    }

}


static
void rotate_polar_coord(std::vector<int>& clusters, int nCluster,\
    std::vector<WgtType>& rotate_degree,\
    std::vector<WgtType>& theta_coord)
{
    int vtx_c;
    for (int i=0; i<theta_coord.size(); ++i)
    {
        vtx_c = clusters.at(i);
        theta_coord.at(i) += rotate_degree.at(vtx_c);
    }
}


static
void euclidean_transformation(std::vector<int>& clusters, int nCluster,\
    std::vector< std::vector<CoordType> >& center_coord,\
    std::vector<WgtType>& radius_coord,\
    std::vector<WgtType>& theta_coord,\
    std::vector< std::vector<CoordType> >& coord)
{
    using namespace std;
    CoordType center_x;
    CoordType center_y;
    CoordType vtx_x;
    CoordType vtx_y;
    int vtx_c;
    for (int i=0; i<coord.size(); ++i)
    {
        vtx_c = clusters.at(i);
        center_x = center_coord.at(vtx_c).at(0);
        center_y = center_coord.at(vtx_c).at(1);
        vtx_x = center_x + radius_coord.at(i) * cos(theta_coord.at(i) * PI / 180);
        vtx_y = center_y + radius_coord.at(i) * sin(theta_coord.at(i) * PI / 180);
        coord.at(i).at(0) = vtx_x;
        coord.at(i).at(1) = vtx_y;

    }
}


static
void rotate_refinement(Graph::Graph& g,\
    std::vector<int>& clusters, int nCluster,\
    std::vector< std::vector<CoordType> >& center_coord,\
    std::vector< WgtType >& radii,\
    std::vector< std::vector<CoordType> >& coord)
{
    // Steps
    // 1. get the radius and theta coordinates
    // 2. calculate the sum of difference
    // 3. rotate the coordinates based on centers
    // 4. transform polar coordinates to euclidean coordinates

    // Step 1
    using namespace std;
    vector<WgtType> radius_coord(coord.size());
    vector<WgtType> theta_coord(coord.size());
    polar_transformation(clusters, nCluster, center_coord, coord,\
                        radius_coord, theta_coord);

    // Step 2
    vector<WgtType> rotate_degree(nCluster, 0); // in radian
    calculate_rotate_degree(g, clusters, nCluster,\
                            center_coord, radii, coord, theta_coord,\
                            rotate_degree);
    
    cout << "rotate degree" << endl;
    for (int c=0; c<rotate_degree.size(); ++c)
    {
        cout << "cluster #" << c << endl;
        cout << rotate_degree.at(c) << endl;
    }
    cout << endl;
    cout << "theta before rotate" << endl;
    for (int c=0; c<theta_coord.size(); ++c)
    {
        cout << theta_coord.at(c) << " ";
    }
    cout << endl;

    // Step 3
    rotate_polar_coord(clusters, nCluster, rotate_degree, theta_coord);
    cout << "theta after rotate" << endl;
    for (int c=0; c<theta_coord.size(); ++c)
    {
        cout << theta_coord.at(c) << " ";
    }
    cout << endl;

    // Step 3
    euclidean_transformation(clusters, nCluster,\
        center_coord, radius_coord, theta_coord, coord);
}


/******************************************************************************
 *                           Main Process                                     *
 ******************************************************************************/
int layout_refinement(Graph::Graph& g,\
    DenseMat& distMat,\
    std::vector<int>& clusters, int nCluster,\
    std::vector< std::vector<CoordType> >& center_coord,\
    std::vector< WgtType >& radii,\
    std::vector< std::vector<CoordType> >& coord)
{
    // Steps
    // 1. Push the vertices which have more inter-connections outside; ones
    //     which have less inter-connection inside
    // 2. Shift (rotate) the vertices to the ideal inter clusters connection
    //     position

    // Step 1
    // refine with input center_coord
    // radial_refinement(g, distMat, clusters, nCluster, coord);

    rotate_refinement(g, clusters, nCluster, center_coord, radii, coord);

    return SUCCESS_LAYOUT_REFINEMENT;

}