#ifndef PBA_H
#define PBA_H


#include <vector>
#include <math.h>

#include "config.h"
#include "graph.h"
#include "pgraph.h"


void fill_edges(Graph::Graph& g, std::vector< std::vector<VtxType> > & edges);
void define_ports_and_boundary_pts(std::vector< std::vector<VtxType> >& edges,\
    std::vector<int>& partition,\
    std::vector< std::vector<VtxType> >& ports,\
    std::vector< std::vector<VtxType> >& boundary_pts);
void find_boundary_point(PGraph::PGraph & pg,\
    std::vector<int>& partition,\
    std::vector< std::vector<CoordType> >& coord,\
    VtxType corrPt,\
    std::vector<CoordType>& boundary_pt_coord);
void find_port(PGraph::PGraph & pg,\
    std::vector<int>& partition,\
    VtxType endPt1, VtxType endPt2,\
    std::vector< std::vector<CoordType> >& coord,\
    std::vector<CoordType>& port1_coord,\
    std::vector<CoordType>& port2_coord);
void assign_ports_and_boundary_pts(PGraph::PGraph & pg,\
    std::vector< std::vector<VtxType> >& edges,\
    std::vector<int>& partition,\
    std::vector< std::vector<VtxType> >& ports,\
    std::vector< std::vector<VtxType> >& boundary_pts,\
    std::vector< std::vector<CoordType> > & ports_coords,\
    std::vector< std::vector<CoordType> > & boundary_pts_coords);
void calculate_control_points(PGraph::PGraph& pg,\
    std::vector<int>& partition,\
    std::vector< std::vector<CoordType> > & coord,\
    std::vector< std::vector<VtxType> > & edges,\
    std::vector< std::vector<VtxType> > & ports,\
    std::vector< std::vector<VtxType> > & boundary_pts,\
    std::vector< std::vector<CoordType> > & ports_coords,\
    std::vector< std::vector<CoordType> > & boundary_pts_coords,\
    std::vector< std::vector<CoordType> >& ctrl_pts_coords);

int port_and_boundary_assignment(Graph::Graph& g, PGraph::PGraph& pg,\
    std::vector<int>& partition,\
    std::vector< std::vector<CoordType> > & coord,\
    std::vector< std::vector<VtxType> > & edges,\
    std::vector< std::vector<VtxType> > & ports,\
    std::vector< std::vector<VtxType> > & boundary_pts,\
    std::vector< std::vector<CoordType> > & ports_coords,\
    std::vector< std::vector<CoordType> > & boundary_pts_coords,\
    std::vector< std::vector<CoordType> >& ctrl_pts_coords);

#endif