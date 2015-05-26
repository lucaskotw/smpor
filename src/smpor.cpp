#include <math.h>
#include "smpor.h"
#include "distance.h"


/******************************************************************************
 *                       Stress Majorization                                  *
 *                         inside clusters                                    *
 ******************************************************************************/
static
void create_cluster_dist_mat(DenseMat& distMat,\
                            std::vector<VtxType>& cluster_vtxs,\
                            DenseMat& clusterDistMat)
{
    // Steps
    // 1. resize cluster distance matrix
    // 2. match the distance from distMat to clusterDistMat

    // Step 1
    clusterDistMat.resize(cluster_vtxs.size(), cluster_vtxs.size());

    // Step 2
    int rowIdx;
    int colIdx;
    for (int c=0; c<cluster_vtxs.size(); ++c)
    {
        for (int r=0; r<cluster_vtxs.size(); ++r)
        {
            rowIdx = cluster_vtxs.at(r);
            colIdx = cluster_vtxs.at(c);
            clusterDistMat(r, c) = distMat(rowIdx, colIdx);
        }
    }
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
int intra_stress_majorization(std::vector<int>& clusters, int nCluster,\
                            DenseMat& distMat,\
                            std::vector< std::vector<CoordType> >& coord)
{
    // Steps
    // Forloop from 0 to nCluster-1
    //     1. get the nodes inside cluster O(N)
    //     2. layout inside cluster: O(stress majorization layout)
    //     3. put the coordinates to the coord: O(N)

    std::vector< std::vector<CoordType> > intra_coord;
    std::vector<VtxType> cluster_vtxs; // id: vtx id of clustered graph
                                       // val: vtx id in whole graph
    DenseMat cluster_dist_mat;
    
    for (int p=0; p<nCluster; ++p)
    {
        // loop step 1
        cluster_vtxs.resize(0);
        for (int i=0; i<clusters.size(); ++i)
            if (clusters.at(i) == p) cluster_vtxs.push_back(i);
        intra_coord.resize(cluster_vtxs.size(), std::vector<CoordType>(2));

        // loop step 2
        create_cluster_dist_mat(distMat, cluster_vtxs, cluster_dist_mat);
        stress_majorization(cluster_dist_mat.rows(), cluster_dist_mat,\
                            intra_coord);

        // loop step 3
        match_cluster_coord(cluster_vtxs, intra_coord, coord);
       
    }

    return SUCCESS_INTRA_STRESS_MAJORIZATION;
}

/******************************************************************************
 *                   Get the radii and centers coordinates                                  *
 ******************************************************************************/
static
WgtType normTwo(std::vector< CoordType >& pt1, std::vector< CoordType >& pt2)
{
    // process 2-norm of given input
    WgtType res;
    res = std::pow(pt1.at(0)- pt2.at(0), 2) + std::pow(pt1.at(1)- pt2.at(1), 2);
    res = std::sqrt(res);

    return res;
}


static
int get_radii_centers(std::vector<int>& clusters, int nCluster,\
    std::vector< std::vector<CoordType> >& coord,\
    std::vector< WgtType >& radii,\
    std::vector< std::vector<CoordType> >& center_coord)
{
    // Steps
    // Forloop from 0 to nCluster-1
    //     1. extract the clustered vtx coord
    //     2. get center
    //     3. get the radii
    std::vector< std::vector<CoordType> > intra_coord;
    std::vector<CoordType> vtx_coord(2);
    std::vector<CoordType> center(2);
    WgtType radius = -1; // radius default value
    WgtType candi_radius; // candidate radius value

    for (int p=0; p<nCluster; ++p)
    {
        // loop step 1
        intra_coord.resize(0, std::vector<CoordType>(2));
        for (int i=0; i<clusters.size(); ++i)
        {
            if (clusters.at(i) == p)
            {
                vtx_coord.at(0) = coord.at(i).at(0);
                vtx_coord.at(1) = coord.at(i).at(1);
                intra_coord.push_back(vtx_coord);

            }
        }

        // loop step 2
        center.clear();
        center.resize(2, 0.0);
        std::cout << "current center" << std::endl;
        for (std::vector<CoordType>::iterator it=center.begin();\
    it!=center.end();
    ++it)
    {
        std::cout << *it << " ";
    }
    std::cout << std::endl;
        for (int i=0; i<intra_coord.size(); ++i)
        {
            center.at(0) += intra_coord.at(i).at(0);
            center.at(1) += intra_coord.at(i).at(1);
        }
        
        center.at(0) /= intra_coord.size();
        center.at(1) /= intra_coord.size();
        center_coord.at(p) = center;
        
        // loop step 3
        for (int i=0; i<intra_coord.size(); ++i)
        {
            candi_radius = normTwo(intra_coord.at(i), center);
            if (candi_radius > radius) radius = candi_radius;
        }
        radii.at(p) = radius;


    }

    return SUCCESS_GET_RADII_CENTERS;
}


/******************************************************************************
 *                           Main Process                                     *
 ******************************************************************************/
int smpor(int graphSize, DenseMat& distMat,\
    std::vector< std::vector<CoordType> >& coord,\
    std::vector< std::vector<CoordType> >& center_coord,\
    std::vector< WgtType >& radii,\
    std::vector<int>&  clusters, int nCluster)
{

    std::vector<Graph> sg_vec(nCluster, Graph(0));

    // Steps
    // 1. stress majorization inside clusters.
    // 2. get the radii and centers coordinates.
    // 3. stress majorization between clusters.
    // 4. shift nodes with cluster respect to the centers.

    // Step 1
    intra_stress_majorization(clusters, nCluster, distMat, coord);

    // Step 2
    get_radii_centers(clusters, nCluster, coord, radii, center_coord);



    // // draw the layout of p-graph: give the coordinates
    // std::vector< std::vector<CoordType> > pg_coord(pg.get_num_vtxs(),\
    //                                     std::vector<CoordType> (2, INIT_COORD));
    // stress_majorization_with_node_overlap_removal(pg, pg_coord);


    // matching the clusters coordinates to the overall coordinates
    // clusters graph for easily p-center
    // intra_stress_majorization(sg_vec, pg, clusters, pg_coord,\
    //                                    coord, nCluster);

    return SUCCESS_SMPOR;


}