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


/******************************************************************************
 *                           Coarsening Phase                                 *
 ******************************************************************************/
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


/******************************************************************************
 *                           Partitioning Phase                               *
 ******************************************************************************/

bool is_partition_num_exceed(std::vector<int>& partition, int partID, int partNum)
{
    int part_id_num = 0;
    for (std::vector<int>::iterator it=partition.begin(); it!=partition.end(); ++it)
    {
        if ( (*it) == partID) ++part_id_num;
    }
    if (part_id_num > partition.size()/partNum+1)
    {
        std::cout << "partition threshold = " << partition.size()/partNum+1 << std::endl;
        std::cout << "partition id number = " << part_id_num << std::endl;
        return true;
    }
    else
    {
        return false;
    }
}

bool in_partition(VtxType vtx, std::vector<VtxType>& partSet)
{
    for (std::vector<VtxType>::iterator it=partSet.begin();\
        it!=partSet.end();\
        ++it)
    {
        if ( vtx == (*it) ) return true;
    }
    return false;
}


bool is_adj_lock(std::vector<VtxType>& part_order, \
                std::vector<bool>& lock, \
                VtxType adjVtx)
{
    int lock_id = 0;
    for (int i=0; i<part_order.size(); ++i)
    {
        if (adjVtx == part_order.at(i))
        {
            lock_id = i;
            break;
        }
    }
    if (lock.at(lock_id)) return true;
    else return false;
}


void compute_gain(CGraph::CGraph& coarsetCG, \
    std::vector<VtxType>& part_order, \
    std::vector<WgtType>& gain, \
    std::vector<bool>& lock)
{
    // reset gain with lock is false
    for (int i=0; i<gain.size(); ++i)
    {
        if (!lock.at(i)) gain.at(i) = 0;
    }


    // record the partition vertices
    std::vector<VtxType> part_a;
    for (int a=0; a<part_order.size()/2; ++a) part_a.push_back(part_order.at(a));
    std::vector<VtxType> part_b;
    for (int b=part_order.size()/2; b<part_order.size(); ++b)
        part_b.push_back(part_order.at(b));

    // Compute the gain
    std::vector<VtxType> adj;
    std::vector<WgtType> adj_wgts;
    int                  adj_idx;

    // 1) Partition A gain
    for (int a=0; a<part_a.size(); ++a)
    {
        if (!lock.at(a))
        {
            adj = coarsetCG.adj(part_a.at(a));
            adj_wgts = coarsetCG.adj_wgts(part_a.at(a));
            adj_idx = 0;
            for (std::vector<VtxType>::iterator it=adj.begin(); it!=adj.end(); ++it)
            {
                if ( !is_adj_lock(part_order, lock, (*it)) )
                {
                    if ( in_partition((*it), part_b) )
                    {
                        gain.at(a) += adj_wgts.at(adj_idx);
                    }
                    if ( in_partition((*it), part_a) )
                    {
                        gain.at(a) -= adj_wgts.at(adj_idx);
                    }
                }
                ++adj_idx;
            }
        }
    }
    // 2) Partition B gain
    for (int b=0; b<part_b.size(); ++b)
    {
        if (!lock.at(b+gain.size()/2))
        {
            adj = coarsetCG.adj(part_b.at(b));
            adj_wgts = coarsetCG.adj_wgts(part_b.at(b));
            adj_idx = 0;
            for (std::vector<VtxType>::iterator it=adj.begin(); it!=adj.end(); ++it)
            {
                if ( !is_adj_lock(part_order, lock, (*it)) )
                {
                    if ( in_partition((*it), part_a) )
                    {
                        gain.at(b+part_a.size()) += adj_wgts.at(adj_idx);
                    }
                    if ( in_partition((*it), part_b) )
                    {
                        gain.at(b+part_a.size()) -= adj_wgts.at(adj_idx);
                    }
                }
                ++adj_idx;
            }
        }
    }
}



void find_max_gain(CGraph::CGraph& coarsetCG, \
                    std::vector<VtxType>& partOrder, \
                    std::vector<WgtType>& gain, \
                    std::vector<bool>& lock, \
                    WgtType& d, \
                    VtxType& a, \
                    VtxType& b)
{
    // Select the gain that is maximum in first half in gain
    // and second half in gain
    // 1) part A
    int max_a_gain = MAX_GAIN_INIT_VAL;
    int max_a_idx = 0;
    for (int i=0; i<gain.size()/2; ++i)
    {
        if ((gain.at(i) > max_a_gain) && (!lock.at(i))) // if not lock, then select
        {
            max_a_gain = gain.at(i);
            max_a_idx = i;
        }
    }

    // 2) part B
    int max_b_gain = MAX_GAIN_INIT_VAL;
    int max_b_idx = 0;
    for (int i=gain.size()/2; i<gain.size(); ++i)
    {
        if ((gain.at(i) > max_b_gain) && (!lock.at(i)))
        {
            max_b_gain = gain.at(i);
            max_b_idx = i;
        }
    }

    // assign a, b to corresponding vertex id in cg
    a = partOrder.at(max_a_idx);
    b = partOrder.at(max_b_idx);
    // calculate difference and assign to d
    if (coarsetCG.is_edge(a, b)){
        d = max_a_gain + max_b_gain - 2*coarsetCG.get_pwgt(a, b);
    }
    else
    {
        d = max_a_gain + max_b_gain;
    }
    
    // set both vertices idx in lock true
    lock.at(max_a_idx) = true;
    lock.at(max_b_idx) = true;
}


void find_k_max_diff(std::vector<WgtType>& diff, \
                    WgtType& diff_max, \
                    int& k)
{
    WgtType diff_curr_sum = 0;
    diff_max = 0;
    int i = 0;
    for (std::vector<WgtType>::iterator it=diff.begin(); it!=diff.end(); ++it)
    {
        diff_curr_sum += (*it);
        if (diff_curr_sum > diff_max)
        {
            diff_max = diff_curr_sum;
            k = i;
        }
        ++i;
    }

}


void exchange_av_with_bv(std::vector<VtxType>& partOrder, \
                        std::vector<VtxType>& av, \
                        std::vector<VtxType>& bv, \
                        int k)
{
    int a_idx;
    int b_idx;
    int tmp_vtx;

    for (int i=0; i<k+1; ++i) // because k is a index, starting from 0
    {
        // get a, b index in partOrder
        for (int pi=0; pi<partOrder.size(); ++pi)
        {
            if (av.at(i) == partOrder.at(pi)) a_idx = pi;
            if (bv.at(i) == partOrder.at(pi)) b_idx = pi;
        }

        // swap order
        tmp_vtx = partOrder.at(a_idx);
        partOrder.at(a_idx) = partOrder.at(b_idx);
        partOrder.at(b_idx) = tmp_vtx;
    }
}

/*
 * The variables naming are based on KK99
 */
void kl_bisection(CGraph::CGraph& coarsetCG, \
    std::vector<int>& coarsetPartition, \
    int partID, int nextPartID)
{
    // [Initial step] randomly giving the partition id
    // vertices id in part_order <= |part_order|/2 -> first partition
    // else -> second parition
    std::vector<VtxType> part_order;
    for (int i=0; i<coarsetCG.get_num_vtxs(); ++i)
    {
        // put the vertex inside the order if it's in the expected partition id
        if (coarsetPartition.at(i) == partID) part_order.push_back(i);
    }
    srand(1);
    std::random_shuffle(part_order.begin(), part_order.end());

    // [Iteration step] kl iteration
    std::vector<WgtType> gain(part_order.size());
    std::vector<bool>    lock(part_order.size());
    std::vector<VtxType> av;
    std::vector<VtxType> bv;
    std::vector<WgtType> diff; // term to decide diff is clean or not
    WgtType diff_max = 0;
    int            k = 0;
    WgtType        d;
    VtxType        a;
    VtxType        b;

    do
    {
        std::fill(lock.begin(), lock.end(), false);
        compute_gain(coarsetCG, part_order, gain, lock);
        av.resize(0);
        bv.resize(0);
        diff.resize(0);

        for (int i=0; i<part_order.size()/2; ++i)
        {

            find_max_gain(coarsetCG, part_order, gain, lock, d, a, b);

            diff.push_back(d);
            av.push_back(a);
            bv.push_back(b);

            compute_gain(coarsetCG, part_order, gain, lock);
            

        }

        find_k_max_diff(diff, diff_max, k);

        if (diff_max > 0)
        {
            exchange_av_with_bv(part_order, av, bv, k);
        }
    } while (diff_max > 0);

    // assign those belong to next partition
    for (int i=part_order.size()/2; i<part_order.size(); ++i)
    {
        for (int vi=0; vi<coarsetPartition.size(); ++vi)
        {
            if (part_order.at(i) == vi)
            {
                coarsetPartition.at(vi) = nextPartID;
            }
        }
    }

}


int partitioning_phase(CGraph::CGraph& coarsetCG, \
    std::vector<int>& coarsetPartition, \
    int partNum)
{
    coarsetPartition.resize(coarsetCG.get_num_vtxs());
    std::fill(coarsetPartition.begin(), coarsetPartition.end(), 0);
    int next_part_idx = 0;

    // examine whether the number of vertices in current partition exceeding the
    // expected number (= n/k, where k is partition number)
    for (int p=0; p<partNum; ++p)
    {
        while (is_partition_num_exceed(coarsetPartition, p, partNum))
        {
            
            kl_bisection(coarsetCG, coarsetPartition, p, ++next_part_idx);
            std::cout << "work on partition " << p << " to " << next_part_idx << std::endl;
            std::cout << "current coarset partition" << std::endl;
            for (std::vector<VtxType>::iterator it=coarsetPartition.begin(); \
                it!=coarsetPartition.end();\
                ++it)
            {
                std::cout << (*it) << " ";
            }
            std::cout << std::endl;

        }

    }
    // kl_bisection(coarsetCG, coarsetPartition, 0, 1);
    // std::cout << "current coarset partition" << std::endl;
    // for (std::vector<VtxType>::iterator it=coarsetPartition.begin(); \
    //     it!=coarsetPartition.end();\
    //     ++it)
    // {
    //     std::cout << (*it) << " ";
    // }
    // std::cout << std::endl;


    return SUCCESS_INIT_PARTITION;

}



int partition_graph (Graph::Graph& g, std::vector<int>& partition, int partNum)
{

    std::vector<CGraph::CGraph> cgs;
    coarsening_phase(g, cgs, partNum);

    std::vector<int> coarset_partition;
    partitioning_phase(cgs.at(cgs.size()-1), coarset_partition, partNum);
    // uncoarsening_phase();

    return SUCCESS_PARTITION;
}