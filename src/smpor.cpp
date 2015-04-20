#include "smpor.h"


/******************************************************************************
 *                          Create P-Graph                                    *
 ******************************************************************************/
/*
 * Use naive way to find the center of the partition
 * minimal of maximal distance
 */
void find_p_center_radius(Graph::Graph& g, std::vector<int>& partition,\
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
        // print out dist val
        std::cout << "bfs dist from source = " << p_v << std::endl;
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


int create_pgraph(PGraph::PGraph& pg, Graph::Graph& g,\
    std::vector<int>& partition, int partNum)
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
        find_p_center_radius(g, partition, p, p_center, p_radius);
        std::cout << "partition #" << p << "'s center = " << p_center << std::endl;
        std::cout << "partition #" << p << "'s radius = " << p_radius << std::endl;
        pg.add_node(p_center, p_radius);
    }

    // add p-graph's edges
    add_p_graph_edges(pg, g);

    std::cout << "pg size = " << partNum << std::endl;
    pg.print_graph();
    return SUCCESS_CREATE_PGRAPH;
}


/******************************************************************************
 *                           Main Process                                     *
 ******************************************************************************/
int smpor(Graph::Graph& g, std::vector< std::vector<CoordType> >& coord,\
    std::vector<int>& partition, int partNum)
{
    PGraph pg;
    
    create_pgraph(pg, g, partition, partNum);

    return SUCCESS_SMPOR;


}