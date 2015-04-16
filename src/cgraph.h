#ifndef CGRAPH_H
#define CGRAPH_H


#include <vector>
#include <iostream>


#define VtxType  int     // type of vertex id
#define WgtType  double  // type of edge weight value

#define NullVtx  -1


struct m_node
{
    /*********************
     * Multinode struct
     *
     * Different cases:
     * 1. G_0: pu and pv point to NullVtx
     * 2. those vertics which are not concatenated by two previous point:
     *    1) pu -> previous id
     *    2) pv -> NullVtx
     */
    VtxType pu;                 // previous one end point
    VtxType pv;                 // previous another end point
    int nEdges;
    std::vector<VtxType> edges; // each edges[0, ..., k] represent the vertex
                                // connects to this one
    std::vector<WgtType> pWgts; // prefered edges weights corresponds to 
                                // edges[0, ..., k]

};


class CGraph
{
    private:
        std::vector<m_node> vtxs;
    public:
        CGraph();
        ~CGraph();
        void add_m_node(VtxType pu, VtxType pv);
        void add_edge(VtxType u, VtxType v);
        void add_edge_weight(int u, int edgeId, WgtType wgt);
        std::vector<VtxType> adj(VtxType s);
        std::vector<WgtType> adj_wgts(VtxType s);
        VtxType get_edge_id(VtxType vtx, VtxType corrVtx);
        int get_num_vtxs();
        void print_graph();
};

#endif