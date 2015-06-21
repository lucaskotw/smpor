#include <math.h>
#include <algorithm>  // std::find, std::max_element
#include "smpor.h"
#include "distance.h"


#define RADIAL_COEFF 1


/******************************************************************************
 *                 Stress Majorization inside clusters                        *
 ******************************************************************************/
static
void create_cluster_dist_mat(
    DenseMat& distMat,
    std::vector<VtxType>& cluster_vtxs,
    std::vector<WgtType>& vtx_radii,
    DenseMat& clusterDistMat)
{
    // Steps
    // 1. resize cluster distance matrix
    // 2. match the distance from distMat to clusterDistMat
    // 3. add central node influence

    // Step 1
    clusterDistMat.resize(cluster_vtxs.size()+1, cluster_vtxs.size()+1);
    clusterDistMat.fill(0);

    // Step 2
    int rowIdx;
    int colIdx;
    for (int c=0; c<cluster_vtxs.size(); ++c)
    {
        for (int r=0; r<cluster_vtxs.size(); ++r)
        {
            rowIdx = cluster_vtxs.at(r);
            colIdx = cluster_vtxs.at(c);
            clusterDistMat(r+1, c+1) = distMat(rowIdx, colIdx);
        }
    }

    // Step 3
    for (int i=0; i<vtx_radii.size(); ++i)
    {
        clusterDistMat(0, i+1) = vtx_radii.at(i);
        clusterDistMat(i+1, 0) = vtx_radii.at(i);

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
        coord.at(vtx_id).at(0) = intra_coord.at(i+1).at(0);
        coord.at(vtx_id).at(1) = intra_coord.at(i+1).at(1);
    }
}


static
void calculate_vtx_radii(Graph::Graph& g,
    DenseMat& distMat,
    std::vector<int>& clusters,
    std::vector<VtxType>& cluster_vtxs,
    std::vector< std::vector<CoordType> >& intra_coord,
    std::vector<WgtType>& vtx_radii)
{
    // Steps
    // 1. Calculate cluster radius
    // 2. get the inter/intra cluster degree of each vtxs
    // 3. calculate radial constriants
    using namespace std;

    // Step 1
    // get the corresponding cluster distance
    // minus central node
    DenseMat clsDistMat(intra_coord.size()-1, intra_coord.size()-1);
    VtxType rr;
    VtxType cc;
    for (int c=0; c<clsDistMat.cols(); ++c)
    {
        for (int r=0; r<clsDistMat.rows(); ++r)
        {
            rr = cluster_vtxs.at(r);
            cc = cluster_vtxs.at(c);
            clsDistMat(c, r) = distMat(cc, rr);
        }
    }
    cout << "cls distance matrix" << endl;
    cout << clsDistMat << endl;
    // find the maximum pair distance
    double cls_radius = clsDistMat.maxCoeff()/2;

    // Step 2
    vector<VtxType> intra_deg(cluster_vtxs.size(), 0);
    vector<VtxType> inter_deg(cluster_vtxs.size(), 0);
    vector<VtxType> nbors;
    for (int i=0; i<cluster_vtxs.size(); ++i)
    {
        nbors = g.adj( cluster_vtxs.at(i) );
        for (int n=0; n<nbors.size(); ++n)
        {
            if (clusters.at( cluster_vtxs.at(i) ) == clusters.at( nbors.at(n) ))
            {
                intra_deg.at(i) += 1;
            }
            else
            {
                inter_deg.at(i) += 1;
            }
        }
    }
    cout << "intra_deg=" << intra_deg.size() << endl;
    cout << "inter_deg=" << inter_deg.size() << endl;


    // Step 3
    // calculate centrality
    vector<WgtType> centrality(vtx_radii.size());
    const int offset = 1;
    for (int i=0; i<vtx_radii.size(); ++i)
    {
        centrality.at(i) = (intra_deg.at(i)+offset)/(inter_deg.at(i)+offset);
    }
    cout << "centrality size=" << centrality.size() << endl;

    // calculate radius
    double min_centrality = *min_element(centrality.begin(), centrality.end());
    double max_centrality = *max_element(centrality.begin(), centrality.end());
    for (int i=0; i<vtx_radii.size(); ++i)
    {
        vtx_radii.at(i) = cls_radius*(1 - (centrality.at(i)-min_centrality)/(max_centrality-min_centrality+offset));
    }


}

    
static
int intra_stress_majorization(
    Graph::Graph& g,
    std::vector<int>& clusters,
    int nCluster,
    DenseMat& distMat,
    std::vector< std::vector<CoordType> >& coord,
    double interpolation)
{
    // Steps (with central node added)
    // Forloop from 0 to nCluster-1
    //     1. get the nodes inside cluster: O(V)
    //     2. layout inside cluster:
    //       - calculate vertices distance
    //       - create distance matrix
    //       - run stress majorization
    //     3. put the coordinates to the coord: O(N)
    using namespace std;
    vector< vector<CoordType> > intra_coord;
    vector<VtxType> cluster_vtxs; // id: vtx id of clustered graph
                                       // val: vtx id in whole graph
    vector<WgtType> vtx_radii;
    DenseMat cluster_dist_mat;
    
    for (int p=0; p<nCluster; ++p)
    {
        // loop step 1
        cluster_vtxs.resize(0);
        for (int i=0; i<clusters.size(); ++i)
            if (clusters.at(i) == p) cluster_vtxs.push_back(i);
        intra_coord.resize(cluster_vtxs.size()+1, std::vector<CoordType>(2));

        // random initialization
        intra_coord.at(0).at(0) = 0;
        intra_coord.at(0).at(1) = 0;
        for (int c=0; c<2; c++)
        {
            for (int r=1; r<intra_coord.size(); ++r)
            {
                intra_coord.at(r).at(c) = rand()%100/50.0;
            }
        }
        

        // loop step 2
        vtx_radii.clear();
        vtx_radii.resize(cluster_vtxs.size(), 0);
        std::cout << "intra coord size=" << intra_coord.size() << std::endl;
        calculate_vtx_radii(g, distMat, clusters, cluster_vtxs, intra_coord, vtx_radii);
        cout << "vtx radii size=" << vtx_radii.size() << endl;
        for (int i=0; i<vtx_radii.size(); ++i)
        {
            cout << vtx_radii.at(i) << " ";
        }
        cout << endl;
        std::cout << "finish calculate vtx radii" << std::endl;
        create_cluster_dist_mat(distMat, cluster_vtxs, vtx_radii, cluster_dist_mat);
        stress_majorization_radial_refinement(g.get_num_vtxs(),
                                              cluster_dist_mat,
                                              interpolation,
                                              intra_coord);
        cout << "intra coord" << endl;
        for (vector< vector<CoordType> >::iterator it1=intra_coord.begin();\
            it1!=intra_coord.end();
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
        std::cout << "cluster #" << p << " finish sm" << std::endl;

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
    WgtType radius; // radius default value
    WgtType candi_radius; // candidate radius value

    for (int c=0; c<nCluster; ++c)
    {
        // loop step 1
        intra_coord.clear();
        intra_coord.resize(0, std::vector<CoordType>(2));
        for (int i=0; i<clusters.size(); ++i)
        {
            if (clusters.at(i) == c)
            {
                vtx_coord.at(0) = coord.at(i).at(0);
                vtx_coord.at(1) = coord.at(i).at(1);
                intra_coord.push_back(vtx_coord);

            }
        }

        // loop step 2
        center.clear();
        center.resize(2, 0.0);
        radius = -1;

        for (int i=0; i<intra_coord.size(); ++i)
        {
            center.at(0) += intra_coord.at(i).at(0);
            center.at(1) += intra_coord.at(i).at(1);
        }
        
        center.at(0) /= intra_coord.size();
        center.at(1) /= intra_coord.size();
        center_coord.at(c) = center;
        
        // loop step 3
        for (int i=0; i<intra_coord.size(); ++i)
        {
            candi_radius = normTwo(intra_coord.at(i), center);
            if (candi_radius > radius) radius = candi_radius;
        }
        radii.at(c) = radius;
        std::cout << "radius = " << radii.at(c);


    }

    return SUCCESS_GET_RADII_CENTERS;
}


/******************************************************************************
 *                 Stress Majorization Between Clusters                       *
 ******************************************************************************/
static
int inter_stress_majorization(Graph::Graph& g,\
    std::vector< std::vector<CoordType> >& center_coord,\
    std::vector< WgtType >& radii,\
    std::vector<int>& clusters, int nCluster,\
    std::vector< std::vector<CoordType> >& new_center_coord)
{
    // Steps (whether the first step has improvement? sort node first?)
    // 1. Create the Graph correspond to the centers
    //     * node id will be the clusters id: O(N^2)
    // 2. Calculate distance matrix
    // 3. run the stress majorization of this graph
    // 4. shift the new center coordinates to new position, where the origin is
    //    the `center' point of the centers

    int cg_size = center_coord.size();
    // Step 1
    Graph cg(cg_size);
    bfs_create_clusters_graph(g, g.get_num_vtxs(), 0, radii,\
        clusters, nCluster, cg);

    // Step 2
    DenseMat cDistMat(cg_size, cg_size);
    distance_matrix(cg, cDistMat);
    // std::cout << "cluster graph distance matrix" << std::endl;
    // std::cout << cDistMat << std::endl;

    // Step 3
    // stress_majorization(cDistMat.rows(), cDistMat,\
    //                     new_center_coord);
    new_center_coord = center_coord;
    stress_majorization(cDistMat.rows(), cDistMat,\
                        new_center_coord);

    // Step 4
    std::vector<CoordType> new_center(2);
    for (int i=0; i<new_center_coord.size(); ++i)
    {
        new_center.at(0) += new_center_coord.at(i).at(0);
        new_center.at(1) += new_center_coord.at(i).at(1);
    }    
    new_center.at(0) /= new_center_coord.size();
    new_center.at(1) /= new_center_coord.size();

    for (int i=0; i<center_coord.size(); ++i)
    {
        new_center_coord.at(i).at(0) -= new_center.at(0);
        new_center_coord.at(i).at(1) -= new_center.at(1);
    }

    //
    // reduce_crossing_torque()


    return SUCCESS_INTER_STRESS_MAJORIZATION;
}


/******************************************************************************
 *                 Stress Majorization Between Clusters                       *
 ******************************************************************************/
static
int shift_intra_cluster_vertices(std::vector<int>& clusters,\
    std::vector< std::vector<CoordType> >& new_center_coord,\
    std::vector< std::vector<CoordType> >& center_coord,\
    std::vector< std::vector<CoordType> >& coord)
{
    // Steps
    // 1. get the distance difference between old center of each vtxs
    // 2. shift based on new center_coord

    // Step 1
    int vtx_c;
    std::vector< std::vector<CoordType> > coord_diff(coord.size(),\
                                                std::vector<CoordType>(2, 0));

    for (int i=0; i<coord.size(); ++i)
    {
        vtx_c = clusters.at(i);
        coord_diff.at(i).at(0) = coord.at(i).at(0) - center_coord.at(vtx_c).at(0);
        coord_diff.at(i).at(1) = coord.at(i).at(1) - center_coord.at(vtx_c).at(1);
    }

    // Step 2
    for (int i=0; i<coord.size(); ++i)
    {
        vtx_c = clusters.at(i);
        coord.at(i).at(0) = new_center_coord.at(vtx_c).at(0) + coord_diff.at(i).at(0);
        coord.at(i).at(1) = new_center_coord.at(vtx_c).at(1) + coord_diff.at(i).at(1);
    }

    return SUCCESS_SHIFT_INTRA_CLUSTER;
}


/******************************************************************************
 *                           Main Process                                     *
 ******************************************************************************/
int smpor(Graph::Graph& g, int graphSize, DenseMat& distMat,
    std::vector< std::vector<CoordType> >& coord,
    std::vector< std::vector<CoordType> >& center_coord,
    std::vector< WgtType >& radii,
    std::vector<int>& clusters, int nCluster,
    double interpolation)
{

    // Steps
    // 1. stress majorization inside clusters.
    // 2. get the radii and centers coordinates.
    // 3. stress majorization between clusters.
    // 4. shift vertices within cluster respect to the centers.
    // 5. replace the old center coord by new center coord

    // Step 1
    intra_stress_majorization(g, clusters, nCluster, distMat, coord, interpolation);
    std::cout << "intra sm finish" << std::endl;

    // Step 2
    get_radii_centers(clusters, nCluster, coord, radii, center_coord);
    std::cout << "radius and center calculation finish" << std::endl;
    std::cout << "radius" << std::endl;
    for (std::vector<CoordType>::iterator it2=radii.begin();\
        it2!=radii.end();
        ++it2)
    {
        std::cout << *it2 << " ";
    }
    std::cout << std::endl;

    // Step 3
    std::vector< std::vector<CoordType> > new_center_coord(nCluster,\
                                                std::vector<CoordType>(2));
    inter_stress_majorization(g, center_coord, radii, clusters, nCluster,\
        new_center_coord);
    std::cout << "inter sm finish" << std::endl;

    // Step 4
    shift_intra_cluster_vertices(clusters, new_center_coord, center_coord, coord);
    std::cout << "shift vtxs finish" << std::endl;

    // Step 5
    for (int i=0; i<center_coord.size(); ++i)
    {
        center_coord.at(i).at(0) = new_center_coord.at(i).at(0);
        center_coord.at(i).at(1) = new_center_coord.at(i).at(1);
    }
    std::cout << "center coord" << std::endl;
    for (std::vector< std::vector<CoordType> >::iterator it1=center_coord.begin();\
        it1!=center_coord.end();
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


    return SUCCESS_SMPOR;


}