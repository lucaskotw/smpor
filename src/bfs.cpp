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
    VtxType v;                  // processing vertex
    std::vector<VtxType> nbors; // neighbors of processing vertex
    WgtType farthest_dist;      // current farthest dist from source

    while(!Q.empty())
    {
        v = Q.front();
        Q.pop();

        farthest_dist = dist.at(v);
        nbors = g.adj(v);
        for (int i=0; i<nbors.size(); ++i)
        {
            if ( !explored.at(nbors.at(i)) )
            {
                // [modify!] no weight edge assign
                dist.at(nbors.at(i)) = farthest_dist + 1;  

                Q.push(nbors.at(i));
                explored.at(nbors.at(i)) = true;
            }
        }
    }


    // for non-connected components
    for (int i=0; i<dist.size(); ++i)
    {
        // [todo]: explained 10?
        if (!explored.at(i)) dist.at(i) = farthest_dist + 10;

    }

    return SUCCESS_BFS;    
}