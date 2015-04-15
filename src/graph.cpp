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
std::vector<DistType> Graph::adj(VtxType s)
{
    return vtxs.at(s).edges;
}