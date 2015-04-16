/*
 * Implemenet multilevel k-way partition based on the idea, which proposed by 
 * KaryPis and Kumar in 1998
 *
 * Implement
 *   1) heavy-edge matching
 *   2) KL-heuristic
 * 
 */
#include "partition.h"



int heavy_edge_matching(CGraph& pcg, CGraph& cg)
{
    // 1. given random order to traverse pcg
    std::vector<VtxType> visited_order;
    for (int i=0; i<pcg.get_num_vtxs(); ++i) visited_order.push_back(i);
    srand(1); // set random seed to 1
    std::random_shuffle(visited_order.begin(), visited_order.end());


    // 2. traverse pcg with given order and finding the matching
    //    to create m_node in cg
    std::vector<bool> matched(pcg.get_num_vtxs());     // initialize matched flag
    std::fill(matched.begin(), matched.end(), false);  // vector

    std::vector<VtxType> matching(pcg.get_num_vtxs()); // matching vector to
                                                       // record current matching
                                                       // node
    int                  m_node_idx = 0;
    std::vector<VtxType> adj;
    int                  adj_matched_cnt;
    std::vector<WgtType> adj_wgts;
    std::vector<VtxType> selected_adj;
    int                  selected_adj_idx;
    std::vector<WgtType> selected_adj_wgts;
    m_node vtx;
    VtxType max_wgt_vtx;

    for (std::vector<VtxType>::iterator it1=visited_order.begin();\
        it1!=visited_order.end();\
        ++it1)
    {
        if ( !matched.at( (*it1) ) )
        {

            adj = pcg.adj( (*it1) );
            adj_wgts = pcg.adj_wgts( (*it1) );
            selected_adj.resize(0);
            selected_adj_wgts.resize(0);
            adj_matched_cnt = 0;

            // get valid neighbors
            selected_adj_idx = 0;
            for (std::vector<VtxType>::iterator it2=adj.begin();\
                it2!=adj.end();\
                ++it2)
            {
                
                // only those matched val is false could be selected
                if ( !matched.at( (*it2) ) )
                {

                    selected_adj.push_back( (*it2) );
                    selected_adj_wgts.push_back( adj_wgts.at(selected_adj_idx) );
                    ++adj_matched_cnt;
                    
                }
                ++selected_adj_idx;

            }
            if (adj_matched_cnt != 0)  // make sure there's a neighbor for 
                                       // matching
            {
                select_vertex_with_heaviest_edge(selected_adj,\
                                            selected_adj_wgts,\
                                            max_wgt_vtx);
                cg.add_m_node( (*it1), max_wgt_vtx);

                matched.at( (*it1) ) = true;
                matched.at( max_wgt_vtx ) = true;

                matching.at( (*it1) ) = m_node_idx;
                matching.at( max_wgt_vtx ) = m_node_idx;
                ++m_node_idx;
            }
        }
    }

    // dealing those unmatched node
    int unmatched_idx = 0;
    for (std::vector<bool>::iterator bit=matched.begin();\
        bit!=matched.end();\
        ++bit)
    {
        if ( !(*bit) ) // if unmatched, add as simple one to m_node collection
        {
            cg.add_m_node(unmatched_idx, -1);
            matching.at( unmatched_idx ) = m_node_idx;
            ++m_node_idx;
        }
        ++unmatched_idx;
    }

    // 3. based on pcg, add edge to cg
    std::vector<VtxType> p_adj;
    std::vector<WgtType> p_adj_wgts;
    int                  edge_idx;
    // use to check whether edge is added
    std::vector< std::vector<bool> > edge_added(cg.get_num_vtxs(), \
        std::vector<bool>(cg.get_num_vtxs(), false));

                                                               
    for (int i=0; i<matching.size(); ++i)
    {
        p_adj = pcg.adj(i);
        p_adj_wgts = pcg.adj_wgts(i);
        
        for (int j=0; j<p_adj.size(); ++j)
        {

            if (matching.at(p_adj.at(j)) != matching.at(i))
            {

                if ( !edge_added.at(matching.at(i)).at( matching.at(p_adj.at(j)) ) )
                {

                    cg.add_edge( matching.at(i), matching.at(p_adj.at(j)) );
                    edge_added.at(matching.at(i)).at( matching.at(p_adj.at(j)) ) = true;

                }
                edge_idx = cg.get_edge_id(matching.at(i), matching.at(p_adj.at(j)));
                cg.add_edge_weight( matching.at(i), edge_idx, p_adj_wgts.at(j) );
                
                
            }
        }
        
    }

    return SUCCESS_MATCHING;
}


bool is_coarsening_terminate(CGraph::CGraph& pcg, CGraph::CGraph& cg, int partNum)
{
    if (cg.get_num_vtxs() < THRESHOLD_VTX * partNum)
    {
        return true;
    }
    else if ( ((double)(pcg.get_num_vtxs() - cg.get_num_vtxs()))/pcg.get_num_vtxs() \
        < THRESHOLD_RATIO)
    {
        return true;
    }
    return false;
}


int coarsening_phase (Graph::Graph& g, std::vector<CGraph::CGraph>& cgs, int partNum)
{
    // create g_0 in cgs
    CGraph cg;
    for (int i=0; i<g.get_num_vtxs(); ++i)
    {
        cg.add_m_node(NullVtx, NullVtx);
    }
    // adding the edges
    std::vector<VtxType> adj;
    std::vector<WgtType> adj_wgts;
    int                  edge_id;
    for (int i=0; i<g.get_num_vtxs(); ++i)
    {
        adj = g.adj(i);
        adj_wgts = g.adj_wgts(i);
        for (int j=0; j<g.adj(i).size(); ++j)
        {
            cg.add_edge( i, adj.at(j) );
            edge_id = cg.get_edge_id( i, adj.at(j) );
            cg.add_edge_weight(i, edge_id, adj_wgts.at(j) );
        }
    }
    cg.print_graph();
    cgs.push_back(cg);

    // start the coarsening process
    // avoiding blob
    CGraph pcg = cg;
    cg = CGraph();
    heavy_edge_matching(pcg, cg);
    cg.print_graph();
    cgs.push_back(cg);

std::cout << "cg num = " << cg.get_num_vtxs() << std::endl;
        std::cout << "ratio = " << ((double)(pcg.get_num_vtxs() - cg.get_num_vtxs()))/pcg.get_num_vtxs() << std::endl;
    while( !is_coarsening_terminate(pcg, cg, partNum) )
    {
        std::cout << "cg num = " << cg.get_num_vtxs() << std::endl;
        std::cout << "ratio = " << ((double)(pcg.get_num_vtxs() - cg.get_num_vtxs()))/(double)pcg.get_num_vtxs() << std::endl;
        // clear cg and record pcg for this round
        pcg = cg;
        cg = CGraph();

        // matching process
        heavy_edge_matching(pcg, cg);
        cgs.push_back(cg);

    }

    return SUCCESS_COARSENING;

}



int partition_graph (Graph::Graph& g, std::vector<int>& partition, int partNum)
{
    std::vector<CGraph::CGraph> cgs;
    coarsening_phase(g, cgs, partNum);
    // partitioning_phase();
    // uncoarsening_phase();

    return SUCCESS_PARTITION;
}