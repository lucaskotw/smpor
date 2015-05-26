#include "smpor.h"


/******************************************************************************
 *                          Create P-Graph                                    *
 ******************************************************************************/
/*
 * Use naive way to find the center of the partition
 * minimal of maximal distance
 */
void find_p_center_radius(Graph::Graph& g, std::vector<Graph>& sg_vec,\
                          std::vector<int>& partition,\
                          int partID, VtxType& pCenter, WgtType& pRadius)
{
    
    std::vector<VtxType> partition_vtxs; // id: vtx id of partition graph 
                                         // val: vtx id in whole graph
    for (int i=0; i<partition.size(); ++i)
    {
        if (partition.at(i) == partID) partition_vtxs.push_back(i);
    }

    // create small graph
    Graph sg(partition_vtxs.size());
    VtxType n_i;
    std::vector<VtxType> adj;
    std::vector<WgtType> adj_wgts;
    for (int i=0; i<sg.get_num_vtxs(); ++i)
    {
        adj      = g.adj( partition_vtxs.at(i) );
        adj_wgts = g.adj_wgts( partition_vtxs.at(i) );
        for (int ai=0; ai<adj.size(); ++ai)
        {
            if ( partition.at( partition_vtxs.at(i) ) == partition.at( adj.at(ai) ) )
            {
                n_i = std::find( partition_vtxs.begin(), partition_vtxs.end(), adj.at(ai) ) - partition_vtxs.begin();
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
    pCenter = partition_vtxs.at(pCenter);
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


int create_pgraph(PGraph::PGraph& pg, std::vector<Graph>& sg_vec,\
    Graph::Graph& g, std::vector<int>& partition, int partNum)
{
    // [Create pgraph]
    // 1) with p-center and p-radius assign
    // 2) graph-theoriatical distance will based on the p-center's distance
    //    on given graph

    // add p-graph's p-node
    VtxType p_center;
    WgtType p_radius;
    for (int p=0; p<partNum; ++p)
    {
        find_p_center_radius(g, sg_vec, partition, p, p_center, p_radius);
        pg.add_node(p_center, p_radius);
    }

    // add p-graph's edges
    add_p_graph_edges(pg, g);

    std::cout << "pg size = " << partNum << std::endl;
    pg.print_graph();
    return SUCCESS_CREATE_PGRAPH;
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


    return SUCCESS_PARTITION_PLACEMENT;
}


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


void match_partition_coord(std::vector<int>& partition_vtxs,\
    std::vector< std::vector<CoordType> >& partition_coords,\
    std::vector< std::vector<CoordType> >& pg_coord,\
    std::vector< std::vector<CoordType> >& coord,\
    int partID,\
    int pCenter)
{
    // assign the partition center coord
    int p_center = partition_vtxs.at(pCenter);
    coord.at(p_center).at(0) = pg_coord.at(partID).at(0);
    coord.at(p_center).at(1) = pg_coord.at(partID).at(1);

    // assign the rest
    for (int i=0; i<partition_vtxs.size(); ++i)
    {
        if (partition_vtxs.at(i) != p_center)
        {
            coord.at(partition_vtxs.at(i)).at(0) =\
                partition_coords.at(i).at(0) + coord.at(p_center).at(0);
            coord.at(partition_vtxs.at(i)).at(1) =\
                partition_coords.at(i).at(1) + coord.at(p_center).at(1);
        }
    }
}


int stress_majorization_of_small_graph(std::vector<Graph>& sg_vec,\
                                PGraph::PGraph& pg,\
                                std::vector<int>& partition,\
                                std::vector< std::vector<CoordType> >& pg_coord,\
                                std::vector< std::vector<CoordType> >& coord,\
                                int partNum)
{
    std::vector<VtxType> partition_vtxs; // id: vtx id of partition graph 
                                         // val: vtx id in whole graph
    std::vector< std::vector<CoordType> > partition_coords;
    DenseMat dist_mat;
    DenseMat w_lap;
    VtxType p_center;
    

    for (int p=0; p<partNum; ++p)
    {
        partition_vtxs.resize(0);

        for (int i=0; i<partition.size(); ++i)
            if (partition.at(i) == p) partition_vtxs.push_back(i);
        partition_coords.resize(partition_vtxs.size(), std::vector<CoordType>(2));

        srand(1); // set seed to 1
        for (int i=0; i<partition_coords.size(); ++i) {
            partition_coords.at(i).at(0) = rand()%100/50.0;
            partition_coords.at(i).at(1) = rand()%100/50.0;
        }
        dist_mat.resize(partition_vtxs.size(), partition_vtxs.size());
        w_lap.resize(partition_vtxs.size(), partition_vtxs.size());
        distance_matrix(sg_vec.at(p), dist_mat);
        w_lap_normal(sg_vec.at(p), dist_mat, w_lap);
        if (p==0)
        {
            std::cout << "dist mat" << std::endl;
            std::cout << dist_mat << std::endl;    
        }

        // find pcenter in partition_vtxs
        p_center = std::find( partition_vtxs.begin(), partition_vtxs.end(), pg.get_center_id(p) ) - partition_vtxs.begin();
        std::cout << "p-center in original graph = " << pg.get_center_id(p) << std::endl;
        std::cout << "p-center in partition vtx = " << p_center << std::endl;

        // here, pCenter is the vtx id in partition_vtxs
        
        stress_majorization(sg_vec.at(p), dist_mat, w_lap, partition_coords, p_center);
        if (p==1)
        {
            std::cout << "after sm coord" << std::endl;
            for (std::vector< std::vector<CoordType> >::iterator it1=partition_coords.begin();\
                it1!=partition_coords.end();
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
        }

        match_partition_coord(partition_vtxs, partition_coords, pg_coord, coord,\
                              p, p_center);
    }

    std::cout << "overall sm coord" << std::endl;
    for (std::vector< std::vector<CoordType> >::iterator it1=coord.begin();\
        it1!=coord.end();
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

    return SUCCESS_SMALL_GRAPH_PLACEMENT;
}


/******************************************************************************
 *                           Main Process                                     *
 ******************************************************************************/
int smpor(Graph::Graph& g, PGraph::PGraph& pg,\
    std::vector< std::vector<CoordType> >& coord,\
    std::vector<int>& partition, int partNum)
{

    std::vector<Graph> sg_vec(partNum, Graph(0));

    // stage 
    create_pgraph(pg, sg_vec, g, partition, partNum);



    // // draw the layout of p-graph: give the coordinates
    // std::vector< std::vector<CoordType> > pg_coord(pg.get_num_vtxs(),\
    //                                     std::vector<CoordType> (2, INIT_COORD));
    // stress_majorization_with_node_overlap_removal(pg, pg_coord);


    // matching the partition coordinates to the overall coordinates
    // partition graph for easily p-center
    stress_majorization_of_small_graph(sg_vec, pg, partition, pg_coord,\
                                       coord, partNum);

    return SUCCESS_SMPOR;


}