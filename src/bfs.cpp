#include <queue>
#include "bfs.h"


int bfs(Graph::Graph& g, int gSize, VtxType s, std::vector<WgtType>& dist)
{

    // preprocessing
    std::queue<VtxType> Q;                                  // queue init

    std::vector<bool> explored(gSize, false);               // explored flag
                                                            // init

    std::fill(dist.begin(), dist.end(), VTX_NOT_CONNECTED); // dist init

    // initial step
    Q.push(s);
    explored.at(s) = true;
    dist.at(s) = 0;

    // iterative step
    VtxType v;                       // processing vertex
    std::vector<VtxType> nbors;      // neighbors of processing vertex
    std::vector<WgtType> nbors_wgts; // neighbors' weights of processing vertex
    WgtType farthest_dist;           // current farthest dist from source

    while(!Q.empty())
    {
        v = Q.front();
        Q.pop();

        farthest_dist = dist.at(v);
        nbors = g.adj(v);
        nbors_wgts = g.adj_wgts(v);
        for (int i=0; i<nbors.size(); ++i)
        {
            if ( !explored.at(nbors.at(i)) )
            {
                // [modify!] no weight edge assign
                dist.at(nbors.at(i)) = farthest_dist + nbors_wgts.at(i);  

                Q.push(nbors.at(i));
                explored.at(nbors.at(i)) = true;
            }
        }
    }


    // for non-connected components
    for (int i=0; i<dist.size(); ++i)
    {
        // [todo]: explained 1?
        if (!explored.at(i)) dist.at(i) = farthest_dist + DISCONNECTED_OFFSET;

    }

    return SUCCESS_BFS;
}


int bfs_create_clusters_graph(Graph::Graph& g, int gSize, VtxType s,\
    std::vector< WgtType >& radii,\
    std::vector<int>& clusters, int nCluster, Graph::Graph& cg)
{

    // preprocessing
    std::queue<VtxType> Q;                                  // queue init

    std::vector<bool> explored(gSize, false);               // explored flag
                                                            // init

    int edge_cnt = 0; // dist init
    const double COMPLETE_GRAPH_N_EDGES = (double)(nCluster*(nCluster-1)) / 2;

    std::vector< std::vector<VtxType> > adj_mat(nCluster,\
                                                std::vector<VtxType>(nCluster, 0));

    // initial step
    Q.push(s);
    explored.at(s) = true;


    // iterative step
    VtxType v;                       // processing vertex
    std::vector<VtxType> nbors;      // neighbors of processing vertex
    std::vector<WgtType> nbors_wgts; // neighbors' weights of processing vertex
    VtxType cur_c;                   // current vtx's cluster
    VtxType adj_c;                   // adj vtx's cluster
    WgtType wgt;                     // added edge's weight

    int cnt=0;
    while(!Q.empty() && (edge_cnt <= COMPLETE_GRAPH_N_EDGES))
    {
        v = Q.front();
        Q.pop();

        nbors = g.adj(v);
        // std::cout << "node = " << v << std::endl;
        // std::cout << "neighbors" << std::endl;
        for (std::vector<VtxType>::iterator it2=nbors.begin();\
        it2!=nbors.end();
        ++it2)
        {
            std::cout << *it2 << " ";
        }
        std::cout << std::endl;
        for (int i=0; i<nbors.size(); ++i)
        {
            if ( !explored.at(nbors.at(i)) )
            {
                cur_c = clusters.at( v );
                adj_c = clusters.at( nbors.at(i) );

                if ( (adj_mat.at(cur_c).at(adj_c) == 0) &&\
                     (adj_mat.at(adj_c).at(cur_c) == 0) &&\
                     cur_c != adj_c)
                {
                    wgt = (radii.at(cur_c) + radii.at(adj_c)) * 1.1; // 1.1 is the offset
                    cg.add_edge(cur_c, adj_c, wgt);
                    adj_mat.at(cur_c).at(adj_c) = 1;
                    adj_mat.at(adj_c).at(cur_c) = 1;
                    edge_cnt++;
                }

                Q.push(nbors.at(i));
                // std::cout << "current node cluster = " << cur_c << std::endl;
                // std::cout << "neighbor node cluster = " << adj_c << std::endl;
                explored.at(nbors.at(i)) = true;
                
            }
        }
        ++cnt;
    }

    // std::cout << "loop cnt = " << cnt << std::endl;
    // std::cout << "COMPLETE_GRAPH_N_EDGES = " << COMPLETE_GRAPH_N_EDGES << std::endl;
    // cg.print_graph();


    // std::cout << "center adj mat" << std::endl;
    // for (std::vector< std::vector<VtxType> >::iterator it1=adj_mat.begin();\
    //     it1!=adj_mat.end();
    //     ++it1)
    // {
    //     for (std::vector<VtxType>::iterator it2=(*it1).begin();\
    //     it2!=(*it1).end();
    //     ++it2)
    //     {
    //         std::cout << *it2 << " ";
    //     }
    //     std::cout << std::endl;
    // }

    return SUCCESS_BFS;
}