/**
 * Declare layout drawing header
 */


#ifndef DRAW_LAYOUT_H
#define DRAW_LAYOUT_H




#include <stdlib.h>  // for exit()
#include <math.h>

#include <IL/ilut.h> // for output png


#include "config.h"
#include <vector>
#include <iostream>
#include <GLFW/glfw3.h>


/* Operation */
#define MAX_WIDTH_INIT_VAL  -10000000
#define BEZIER_INCREMENT    0.1

/* Property */
#define WINDOW_WIDTH  800
#define WINDOW_HEIGHT 800

void set_points_attributes();
void draw_vertices(std::vector< std::vector<CoordType> >& coord, \
    std::vector<PartType>& partition);
void draw_line(CoordType p1_x, CoordType p1_y,\
               CoordType p2_x, CoordType p2_y);
void draw_bezier(CoordType p1_x, CoordType p1_y,\
                 CoordType p2_x, CoordType p2_y,\
                 CoordType ctrl_x, CoordType ctrl_y);
void draw_edges(std::vector< std::vector<VtxType> >& edges, \
    std::vector< std::vector<CoordType> >& coord,\
    std::vector< std::vector<VtxType> >& ports,\
    std::vector< std::vector<VtxType> >& boundary_pts,\
    std::vector< std::vector<CoordType> >& ports_coords,\
    std::vector< std::vector<CoordType> >& boundary_pts_coords);
void draw_radius(std::vector< std::vector<CoordType> >& coord,\
                 std::vector< std::vector<CoordType> >& centers,\
                 std::vector<WgtType>& radius);
void draw_layout(std::vector< std::vector<VtxType> >& edges,\
    std::vector< std::vector<CoordType> >& coord,
    std::vector<PartType>& partition,\
    std::vector< std::vector<VtxType> >& ports,\
    std::vector< std::vector<VtxType> >& boundary_pts,\
    std::vector< std::vector<CoordType> >& ports_coords,\
    std::vector< std::vector<CoordType> >& boundary_pts_coords,\
    std::vector< std::vector<CoordType> >& ctrl_pts_coords,\
    std::vector< std::vector<CoordType> >& centers,\
    std::vector<WgtType>& radius);


#endif