#include "pba.h"




/******************************************************************************
 *                            Sub-Function                                    *
 ******************************************************************************/
void fill_edges(Graph::Graph& g, std::vector< std::vector<VtxType> > & edges)
{
    std::vector<VtxType> adj;
    std::vector<VtxType> pair(2);
    for (int vtx=0; vtx<g.get_num_vtxs(); ++vtx)
    {
        adj = g.adj(vtx);
        for (int nb=0; nb<adj.size(); ++nb)
        {
            if (vtx<adj.at(nb))
            {
                pair.at(0) = vtx;
                pair.at(1) = adj.at(nb);
                edges.push_back(pair);
            }
        }
    }
}


void define_ports_and_boundary_pts(std::vector< std::vector<VtxType> >& edges,\
    std::vector<int>& clusters,\
    std::vector< std::vector<VtxType> >& ports,\
    std::vector< std::vector<VtxType> >& boundary_pts)
{
    VtxType end_pt1;
    VtxType end_pt2;
    int ports_cnt = -1;
    int b_pt_cnt  = -1;

    std::vector<VtxType> port_pair(2);
    std::vector<VtxType> b_pts_pair(2);
    for (int e=0; e<edges.size(); ++e)
    {
        end_pt1 = edges.at(e).at(0);
        end_pt2 = edges.at(e).at(1);

        if (clusters.at(end_pt1) == clusters.at(end_pt2))
        {
            port_pair.at(0)  = NULL_PORT;
            port_pair.at(1)  = NULL_PORT;
            b_pts_pair.at(0) = NULL_BOUNDARY_PT;
            b_pts_pair.at(1) = NULL_BOUNDARY_PT;
        }
        else
        {
            port_pair.at(0) = ++ports_cnt;
            port_pair.at(1) = ++ports_cnt;
            b_pts_pair.at(0) = ++b_pt_cnt;
            b_pts_pair.at(1) = ++b_pt_cnt;
        }

        ports.push_back(port_pair);
        boundary_pts.push_back(b_pts_pair);
    }

}


void find_boundary_point(PGraph::PGraph & pg,\
    std::vector<int>& clusters,\
    std::vector< std::vector<CoordType> >& coord,\
    VtxType corrPt,\
    std::vector<CoordType>& boundary_pt_coord)
{
    VtxType center_id;
    center_id = pg.get_center_id( clusters.at(corrPt) );

    std::vector<CoordType> center_coord(2);
    center_coord = coord.at(center_id);

    std::vector<CoordType> corr_pt_coord(2);
    corr_pt_coord = coord.at(corrPt);

    // calculate the angle refer to center of x
    double dist_to_center;
    dist_to_center = sqrt(\
        pow( (corr_pt_coord.at(0)-center_coord.at(0)), 2) + \
        pow( (corr_pt_coord.at(1)-center_coord.at(1)), 2)\
        );
    double cos_theta;
    cos_theta = (corr_pt_coord.at(0)-center_coord.at(0)) / dist_to_center;
    double sin_theta;
    sin_theta = (corr_pt_coord.at(1)-center_coord.at(1)) / dist_to_center;

    // boundary_point assignment
    double radius = pg.get_radius( clusters.at(corrPt) );
    boundary_pt_coord.at(0) = center_coord.at(0) + radius*cos_theta;
    boundary_pt_coord.at(1) = center_coord.at(1) + radius*sin_theta;



}


void find_port(PGraph::PGraph & pg,\
    std::vector<int>& clusters,\
    VtxType endPt1, VtxType endPt2,\
    std::vector< std::vector<CoordType> >& coord,\
    std::vector<CoordType>& port1_coord,\
    std::vector<CoordType>& port2_coord)
{
    std::vector<CoordType> end_pt1_coord;
    end_pt1_coord = coord.at(endPt1);
    std::vector<CoordType> end_pt2_coord;
    end_pt2_coord = coord.at(endPt2);

    // start from endpoint 1
    // calculate the angle refer to endpoint 1
    double dist_to_pt1;
    dist_to_pt1 = sqrt(\
        pow( (end_pt1_coord.at(0)-end_pt2_coord.at(0)), 2) + \
        pow( (end_pt1_coord.at(1)-end_pt2_coord.at(1)), 2)\
        );
    double cos_theta;
    cos_theta = (end_pt1_coord.at(0)-end_pt2_coord.at(0)) / dist_to_pt1;
    double sin_theta;
    sin_theta = (end_pt1_coord.at(1)-end_pt2_coord.at(1)) / dist_to_pt1;

    // assign the port based on radius and offset
    WgtType end_pt1_radius = pg.get_radius( clusters.at(endPt1) );
    port1_coord.at(0) = end_pt1_coord.at(0) +\
        (end_pt1_radius+PORT_BOUNDARY_OFFSET)*cos_theta;
    port1_coord.at(1) = end_pt1_coord.at(1) +\
        (end_pt1_radius+PORT_BOUNDARY_OFFSET)*sin_theta;

    WgtType end_pt2_radius = pg.get_radius( clusters.at(endPt2) );
    port2_coord.at(0) = end_pt2_coord.at(0) -\
        (end_pt2_radius+PORT_BOUNDARY_OFFSET)*cos_theta;
    port2_coord.at(1) = end_pt2_coord.at(1) -\
        (end_pt2_radius+PORT_BOUNDARY_OFFSET)*sin_theta;


}


void assign_ports_and_boundary_pts(PGraph::PGraph & pg,\
    std::vector< std::vector<VtxType> >& edges,\
    std::vector<int>& clusters,\
    std::vector< std::vector<CoordType> >& coord,\
    std::vector< std::vector<VtxType> >& ports,\
    std::vector< std::vector<VtxType> >& boundary_pts,\
    std::vector< std::vector<CoordType> > & ports_coords,\
    std::vector< std::vector<CoordType> > & boundary_pts_coords)
{
    // assign boundary points
    std::vector<CoordType> boundary_pt_coord(2);
    for (int bp=0; bp<boundary_pts.size(); ++bp)
    {
        if (boundary_pts.at(bp).at(0) != NULL_BOUNDARY_PT)
        {
            find_boundary_point(pg, clusters, coord, edges.at(bp).at(0), boundary_pt_coord);
            boundary_pts_coords.push_back(boundary_pt_coord);

            find_boundary_point(pg, clusters, coord, edges.at(bp).at(1), boundary_pt_coord);
            boundary_pts_coords.push_back(boundary_pt_coord);
        }
            
    }

    // assign port points
    std::vector<CoordType> port1_coord(2);
    std::vector<CoordType> port2_coord(2);
    for (int p=0; p<ports.size(); ++p)
    {
        if (ports.at(p).at(0) != NULL_PORT)
        {
            find_port(pg, clusters, edges.at(p).at(0), edges.at(p).at(1),\
                coord, port1_coord, port2_coord);
            ports_coords.push_back(port1_coord);
            ports_coords.push_back(port2_coord);
        }
            
    }


}



/******************************************************************************
 *                       Calculate Contrl Points                              *
 ******************************************************************************/
void calculate_control_points(PGraph::PGraph& pg,\
    std::vector<int>& clusters,\
    std::vector< std::vector<CoordType> > & coord,\
    std::vector< std::vector<VtxType> > & edges,\
    std::vector< std::vector<VtxType> > & ports,\
    std::vector< std::vector<VtxType> > & boundary_pts,\
    std::vector< std::vector<CoordType> > & ports_coords,\
    std::vector< std::vector<CoordType> > & boundary_pts_coords,\
    std::vector< std::vector<CoordType> >& ctrl_pts_coords)
{
    VtxType center_id;
    std::vector<CoordType> center_coord(2);
    double dist_to_center;
    double cos_theta;
    double sin_theta;
    double radius;
    std::vector<CoordType> middle_pt(2);
    std::vector<CoordType> ctrl_pt(2);

    int p_cnt = -1;
    int b_cnt = -1;

    for (int i=0; i<ports.size(); ++i)
    {
        if (ports.at(i).at(0) != NULL_PORT)
        {
            // calculate first middle points
            center_id = pg.get_center_id( clusters.at(edges.at(i).at(0)) );
            center_coord = coord.at(center_id);

            ++p_cnt;
            ++b_cnt;
            middle_pt.at(0) = \
                (ports_coords.at(p_cnt).at(0)+boundary_pts_coords.at(b_cnt).at(0))/2;
            middle_pt.at(1) = \
                (ports_coords.at(p_cnt).at(1)+boundary_pts_coords.at(b_cnt).at(1))/2;
            // start from center of clusters
            // calculate the angle refer to center of clusters
            dist_to_center  = sqrt(\
                pow( (middle_pt.at(0)-center_coord.at(0)), 2) + \
                pow( (middle_pt.at(1)-center_coord.at(1)), 2)\
            );        
            cos_theta = (middle_pt.at(0)-center_coord.at(0)) / dist_to_center;
            sin_theta = (middle_pt.at(1)-center_coord.at(1)) / dist_to_center;
            radius = pg.get_radius( clusters.at(edges.at(i).at(0)) );
            ctrl_pt.at(0) = center_coord.at(0) +\
                (radius+PORT_BOUNDARY_OFFSET)*cos_theta;
            ctrl_pt.at(1) = center_coord.at(1) +\
                (radius+PORT_BOUNDARY_OFFSET)*cos_theta;
            ctrl_pts_coords.push_back(ctrl_pt);



            // calculate second middle points
            center_id = pg.get_center_id( clusters.at(edges.at(i).at(1)) );
            center_coord = coord.at(center_id);
            ++p_cnt;
            ++b_cnt;
            middle_pt.at(0) = \
                (ports_coords.at(p_cnt).at(0)+boundary_pts_coords.at(b_cnt).at(0))/2;
            middle_pt.at(1) = \
                (ports_coords.at(p_cnt).at(1)+boundary_pts_coords.at(b_cnt).at(1))/2;
            // start from center of clusters
            // calculate the angle refer to center of clusters
            dist_to_center  = sqrt(\
                pow( (middle_pt.at(0)-center_coord.at(0)), 2) + \
                pow( (middle_pt.at(1)-center_coord.at(1)), 2)\
            );        
            cos_theta = (middle_pt.at(0)-center_coord.at(0)) / dist_to_center;
            sin_theta = (middle_pt.at(1)-center_coord.at(1)) / dist_to_center;
            radius = pg.get_radius( clusters.at(edges.at(i).at(1)) );
            ctrl_pt.at(0) = center_coord.at(0) +\
                (radius+PORT_BOUNDARY_OFFSET)*cos_theta;
            ctrl_pt.at(1) = center_coord.at(1) +\
                (radius+PORT_BOUNDARY_OFFSET)*cos_theta;
            ctrl_pts_coords.push_back(ctrl_pt);

        }


    }
}


/******************************************************************************
 *                            Main Process                                    *
 ******************************************************************************/
int port_and_boundary_assignment(Graph::Graph& g, PGraph::PGraph& pg,\
    std::vector<int>& clusters,\
    std::vector< std::vector<CoordType> > & coord,\
    std::vector< std::vector<VtxType> > & edges,\
    std::vector< std::vector<VtxType> > & ports,\
    std::vector< std::vector<VtxType> > & boundary_pts,\
    std::vector< std::vector<CoordType> > & ports_coords,\
    std::vector< std::vector<CoordType> > & boundary_pts_coords,\
    std::vector< std::vector<CoordType> >& ctrl_pts_coords)
{
    // ports and boundary_pts
    // 1) id: edge's row id
    // 2) value: ports and boundary_pts id
    // ports_coords and boundary_pts_coords
    // 1) id: ports and boundary_pts id
    // 2) value: x, y coordinates of ports and boundary points


    // fill edges
    fill_edges(g, edges);

    // seperate the process of
    // 1) define ports and boundary-points id
    // 2) coordinates assignments
    define_ports_and_boundary_pts(edges, clusters, ports, boundary_pts);
    assign_ports_and_boundary_pts(pg, edges, clusters, coord,\
        ports, boundary_pts, ports_coords, boundary_pts_coords);

    calculate_control_points(pg, clusters, coord, edges, ports, boundary_pts,\
        ports_coords, boundary_pts_coords, ctrl_pts_coords);


    return SUCCESS_PBA;
}