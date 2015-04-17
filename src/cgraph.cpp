#include "cgraph.h"



/***************
 * Constructor *
 ***************/
CGraph::CGraph()
{
}


/**************
 * Destructor *
 **************/
CGraph::~CGraph()
{
    vtxs.clear();
}



/*****************
 * add multinode *
 *****************/
void CGraph::add_m_node(VtxType pu, VtxType pv)
{
    m_node vtx;
    vtx.pu = pu;
    vtx.pv = pv;
    vtxs.push_back(vtx);
}


/************
 * add edge *
 ************/
void CGraph::add_edge(int u, int v)
{
    vtxs.at(u).edges.push_back(v);
    vtxs.at(u).pWgts.push_back(0); // seperate weight adding with edge adding
    vtxs.at(u).nEdges = vtxs.at(u).edges.size();

    // vtxs.at(v).edges.push_back(u);
    // vtxs.at(v).pWgts.push_back(wgt);
    // vtxs.at(v).nEdges = vtxs.at(v).edges.size();
}


void CGraph::add_edge_weight(int u, int edgeId, WgtType wgt)
{
    vtxs.at(u).pWgts.at(edgeId) += wgt;
}



VtxType CGraph::get_edge_id(VtxType vtx, VtxType corrVtx)
{
    int edge_id = 0;
    for (std::vector<VtxType>::iterator it=vtxs.at(vtx).edges.begin();\
        it!=vtxs.at(vtx).edges.end();\
        ++it)
    {
        if (corrVtx == (*it))
        {
            break;
        }
        else
        {
            ++edge_id;
        }
        
    }

    return edge_id;
}
/*****************************************
 * get the neighbors of the given vertex *
 *****************************************/
std::vector<VtxType> CGraph::adj(VtxType s)
{
    return vtxs.at(s).edges;
}


/*************************************************
 * get the neighbors' weight of the given vertex *
 *************************************************/
std::vector<WgtType> CGraph::adj_wgts(VtxType s)
{
    return vtxs.at(s).pWgts;
}


/******************************
 * get the number of vertices *
 ******************************/
int CGraph::get_num_vtxs()
{
    return vtxs.size();
}


/**********************************************************
 * check whether the given vertices pair is an edge in CG *
 **********************************************************/
bool CGraph::is_edge(VtxType a, VtxType b)
{
    for (int i=0; i<vtxs.at(a).edges.size(); ++i)
    {
        if (vtxs.at(a).edges.at(i) == b) return true;
    }

    return false;
}

/**********************************************************
 * check whether the given vertices pair is an edge in CG *
 **********************************************************/
WgtType CGraph::get_pwgt(VtxType a, VtxType b)
{
    if (is_edge(a, b))
    {
        for (int i=0; i<vtxs.at(a).edges.size(); ++i)
        {
            if (vtxs.at(a).edges.at(i) == b) return vtxs.at(a).pWgts.at(i);
        }
    }
    return -1;
}


/******************
 * show the graph *
 ******************/
void CGraph::print_graph()
{
    std::cout << "corresponding previous node:" << std::endl;
    std::cout << "curr pu pv" << std::endl;
    std::cout << "----------" << std::endl;
    int idx = 0;
    for (std::vector<m_node>::iterator it1=vtxs.begin(); it1!=vtxs.end(); ++it1)
    {
        std::cout << idx << " " << (*it1).pu << " " << (*it1).pv << std::endl;
        ++idx;
    }

    std::cout << "incident edges:" << std::endl;
    idx = 0;
    for (std::vector<m_node>::iterator it1=vtxs.begin(); it1!=vtxs.end(); ++it1)
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

    std::cout << "corresponding edges weight:" << std::endl;
    idx = 0;
    for (std::vector<m_node>::iterator it1=vtxs.begin(); it1!=vtxs.end(); ++it1)
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