#include "pgraph.h"



/***************
 * Constructor *
 ***************/
PGraph::PGraph()
{
}


/**************
 * Destructor *
 **************/
PGraph::~PGraph()
{
    vtxs.clear();
}


/**************
 * add p-node *
 **************/
void PGraph::add_node(VtxType pCenter, WgtType pRadius)
{
    p_node vtx;
    vtx.pCenter = pCenter;
    vtx.pRadius = pRadius;
    vtxs.push_back(vtx);
}


/************
 * add edge *
 ************/
void PGraph::add_edge(VtxType u, VtxType v, WgtType pWgt)
{
    vtxs.at(u).edges.push_back(v);
    vtxs.at(u).pWgts.push_back(pWgt);
}


/**********************************************
 * get the original graph vertex id of p-node *
 **********************************************/
VtxType PGraph::get_center_id(VtxType pID)
{
    return vtxs.at(pID).pCenter;
}


/****************************
 * get the radius of p-node *
 ****************************/
WgtType PGraph::get_radius(VtxType pID)
{
    return vtxs.at(pID).pRadius;
}



/***************************************
 * get the number of vertices in graph *
 ***************************************/
int PGraph::get_num_vtxs()
{
    return vtxs.size();
}


/*****************************************
 * get the neighbors of the given vertex *
 *****************************************/
std::vector<VtxType> PGraph::adj(VtxType s)
{
    return vtxs.at(s).edges;
}


/*************************************************
 * get the neighbors' weight of the given vertex *
 *************************************************/
std::vector<WgtType> PGraph::adj_wgts(VtxType s)
{
    return vtxs.at(s).pWgts;
}


/*******************
 * show the result *
 *******************/
void PGraph::print_graph()
{

    std::cout << "corresponding original vertex:" << std::endl;

    for (int i=0; i<vtxs.size(); ++i)
    {
        std::cout << i << " <--> " << vtxs.at(i).pCenter << std::endl;
    }

    std::cout << "corresponding p-radius:" << std::endl;

    for (int i=0; i<vtxs.size(); ++i)
    {
        std::cout << i << " <--> " << vtxs.at(i).pRadius << std::endl;
    }

    std::cout << "incident edges:" << std::endl;
    int idx = 0;
    for (std::vector<p_node>::iterator it1=vtxs.begin(); it1!=vtxs.end(); ++it1)
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

    std::cout << "corresponding weights:" << std::endl;
    idx = 0;
    for (std::vector<p_node>::iterator it1=vtxs.begin(); it1!=vtxs.end(); ++it1)
    {
        std::cout << idx << " <--> ";
        for (std::vector<WgtType>::iterator it2=(*it1).pWgts.begin(); \
            it2!=(*it1).pWgts.end(); \
            ++it2)
        {
            std::cout << *it2 << " ";
        }
        std::cout << std::endl;
        ++idx;
    }
}