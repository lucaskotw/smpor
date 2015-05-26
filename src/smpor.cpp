#include "smpor.h"
#include "distance.h"


/******************************************************************************
 *                          Create P-Graph                                    *
 ******************************************************************************/
/*
 * Use naive way to find the center of the clusters
 * minimal of maximal distance
 */
void find_p_center_radius(Graph::Graph& g, std::vector<Graph>& sg_vec,\
                          std::vector<int>& clusters,\
                          int partID, VtxType& pCenter, WgtType& pRadius)
{
    
    std::vector<VtxType> cluster_vtxs; // id: vtx id of clusters graph 
                                         // val: vtx id in whole graph
    for (int i=0; i<clusters.size(); ++i)
    {
        if (clusters.at(i) == partID) cluster_vtxs.push_back(i);
    }

    // create small graph
    Graph sg(cluster_vtxs.size());
    VtxType n_i;
    std::vector<VtxType> adj;
    std::vector<WgtType> adj_wgts;
    for (int i=0; i<sg.get_num_vtxs(); ++i)
    {
        adj      = g.adj( cluster_vtxs.at(i) );
        adj_wgts = g.adj_wgts( cluster_vtxs.at(i) );
        for (int ai=0; ai<adj.size(); ++ai)
        {
            if ( clusters.at( cluster_vtxs.at(i) ) == clusters.at( adj.at(ai) ) )
            {
                n_i = std::find( cluster_vtxs.begin(), cluster_vtxs.end(), adj.at(ai) ) - cluster_vtxs.begin();
                if (n_i > i) sg.add_edge( i, n_i, adj_wgts.at(ai) );
            }
        }
        
    }
    sg.print_graph();
    sg_vec.at(partID)=sg;   // push back the current small graph

    // run bfs (other distance algorithm) on vertex
    int sg_size = sg.get_num_vtxs();

    std::vector<WgtType> dist(sg_size);
    std::vector<WgtType> max_dist;


    for (int i=0; i<sg.get_num_vtxs(); ++i)
    {
        bfs(sg, sg_size, i, dist);
        WgtType max_elem = *std::max_element(dist.begin(), dist.end());
        max_dist.push_back(max_elem);
        // // print out dist val
        // std::cout << "bfs dist from source = " << i << std::endl;
        // for (std::vector<WgtType>::iterator it=dist.begin(); it!=dist.end(); ++it)
        // {
        //     std::cout << *it << " ";
        // }
        // std::cout << std::endl;
    }

    pCenter = std::min_element(max_dist.begin(), max_dist.end()) - max_dist.begin();
    pCenter = cluster_vtxs.at(pCenter);
    pRadius = *std::min_element(max_dist.begin(), max_dist.end());
}


void add_p_graph_edges(PGraph::PGraph& pg, Graph::Graph& g)
{
    int orig_g_size = g.get_num_vtxs();
    std::vector<WgtType> dist(orig_g_size);

    for (VtxType p_v=0; p_v<pg.get_num_vtxs(); ++p_v)
    {
        bfs(g, orig_g_size, pg.get_center_id(p_v), dist);
        for (std::vector<WgtType>::iterator it=dist.begin(); it!=dist.end(); ++it)
        {
            std::cout << *it << " ";
        }
        std::cout << std::endl;
        for (VtxType p_nbor=p_v+1; p_nbor<pg.get_num_vtxs(); ++p_nbor)
        {
            pg.add_edge( p_v, p_nbor, dist.at(pg.get_center_id(p_nbor)) );
            pg.add_edge( p_nbor, p_v, dist.at(pg.get_center_id(p_nbor)) );
        }
    }
}


int create_small_graph_list(Graph::Graph& g, std::vector<Graph>& sg_vec,\
    std::vector<int>& clusters, int nCluster)
{
    // Steps
    // 1) with p-center and p-radius assign
    // 2) graph-theoriatical distance will based on the p-center's distance
    //    on given graph

    // add p-graph's p-node
    // VtxType p_center;
    // WgtType p_radius;
    // for (int p=0; p<nCluster; ++p)
    // {
    //     find_p_center_radius(g, sg_vec, clusters, p, p_center, p_radius);
    //     pg.add_node(p_center, p_radius);
    // }

    // // add p-graph's edges
    // add_p_graph_edges(pg, g);

    // std::cout << "pg size = " << nCluster << std::endl;
    // pg.print_graph();
    return SUCCESS_CREATE_SMALL_GRAPH_LIST;
}


/******************************************************************************
 *                       Stress Majorization                                  *
 *                    with node overlap removal                               *
 ******************************************************************************/
void distance_matrix_with_modifed_radius(PGraph::PGraph& pg, DenseMat& distMat)
{
    
    int pg_size = pg.get_num_vtxs();
    std::vector<WgtType> dist(pg_size);

    for (int v=0; v<pg_size; ++v)
    {
        bfs_pg(pg, pg_size, v, dist);

        std::cout << "dist vec of " << v << std::endl;
        for (std::vector<WgtType>::iterator it=dist.begin();\
            it!=dist.end();
            ++it)
        {
            std::cout << *it << " ";
        }
        std::cout << std::endl;

        for (int n_v=0; n_v<pg_size; ++n_v)
        {
            if (v != n_v)
            {
                if ( (pg.get_radius(v)+pg.get_radius(n_v)) > dist.at(n_v) )
                    dist.at(n_v) = pg.get_radius(v)+pg.get_radius(n_v)\
                                 + 2*PORT_BOUNDARY_OFFSET;
                else if ((pg.get_radius(v)+pg.get_radius(n_v)+2*PORT_BOUNDARY_OFFSET) > dist.at(n_v))
                    dist.at(n_v) += 2*PORT_BOUNDARY_OFFSET;
            }

        }
        distMat.row(v) = MapVec(&dist[0], dist.size());
    }


}


int stress_majorization_with_node_overlap_removal(PGraph::PGraph& pg,\
    std::vector< std::vector<CoordType> >& pg_coord)
{
    DenseMat dist_mat(pg.get_num_vtxs(), pg.get_num_vtxs());
    DenseMat w_lap(pg.get_num_vtxs(), pg.get_num_vtxs());

    distance_matrix_with_modifed_radius(pg, dist_mat);
    std::cout << "dist_mat" << std::endl;
    std::cout << dist_mat << std::endl;
    w_lap_with_vertex_radius(pg, dist_mat, w_lap); // put to lap.h
    std::cout << w_lap << std::endl;

    // randomly start coordinates assignment
    srand(1); // set seed to 1
    for (int i=0; i<pg.get_num_vtxs(); ++i) {
        pg_coord.at(i).at(0) = rand()%100/50.0;
        pg_coord.at(i).at(1) = rand()%100/50.0;
    }

    std::cout << "initial coord" << std::endl;
    for (std::vector< std::vector<CoordType> >::iterator it1=pg_coord.begin();\
        it1!=pg_coord.end();
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
    stress_majorization_with_pg(pg, dist_mat, w_lap, pg_coord);
    std::cout << "after sm coord" << std::endl;
    for (std::vector< std::vector<CoordType> >::iterator it1=pg_coord.begin();\
        it1!=pg_coord.end();
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


    return SUCCESS_CLUSTER_PLACEMENT;
}


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
 *                           Main Process                                     *
 ******************************************************************************/
int smpor(int graphSize, DenseMat& distMat,\
    std::vector< std::vector<CoordType> >& coord,\
    std::vector< std::vector<CoordType> >& center_coord,\
    std::vector< WgtType >& radius,\
    std::vector<int>& clusters, int nCluster)
{

    std::vector<Graph> sg_vec(nCluster, Graph(0));

    // Steps
    // 1. stress majorization inside clusters.
    // 2. get the radii and centers coordinates.
    // 3. stress majorization between clusters.
    // 4. shift nodes with cluster respect to the centers.

    // Step 1
    intra_stress_majorization(clusters, nCluster, distMat, coord);



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