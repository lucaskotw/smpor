#include "graph.h"



/***************
 * Constructor *
 ***************/
Graph::Graph(int nNodes)
{
    vtxs.resize(nNodes);
    // fill all vtx as 0
    // std::fill(v.begin(), v.end(), 0);
}


/**************
 * Destructor *
 **************/
Graph::~Graph()
{
    vtxs.clear();
}


/************
 * add edge *
 ************/
void Graph::add_edge(int u, int v, WgtType wgt)
{
    vtxs.at(u).edges.push_back(v);
    vtxs.at(u).pWgts.push_back(wgt);
    vtxs.at(u).nEdges = vtxs.at(u).edges.size();

    vtxs.at(v).edges.push_back(u);
    vtxs.at(v).pWgts.push_back(wgt);
    vtxs.at(v).nEdges = vtxs.at(v).edges.size();
}

/******************************
 * get the number of vertices *
 ******************************/
int Graph::get_num_vtxs()
{
    return vtxs.size();
}


/*****************************************
 * get the neighbors of the given vertex *
 *****************************************/
std::vector<VtxType> Graph::adj(VtxType s)
{
    return vtxs.at(s).edges;
}


/*************************************************
 * get the neighbors' weight of the given vertex *
 *************************************************/
std::vector<WgtType> Graph::adj_wgts(VtxType s)
{
    return vtxs.at(s).pWgts;
}



/******************
 * show the graph *
 ******************/
void Graph::print_graph()
{

    std::cout << "incident edges:" << std::endl;
    int idx = 0;
    for (std::vector<v_dat>::iterator it1=vtxs.begin(); it1!=vtxs.end(); ++it1)
    {
        std::cout << idx << " <--> ";
        for (std::vector<VtxType>::iterator it2=(*it1).edges.begin(); \
            it2!=(*it1).edges.end(); \
            ++it2)
        {
            std::cout << *it2 << " ";
        }
        std::cout << std::endl;
        ++idx;
    }

}