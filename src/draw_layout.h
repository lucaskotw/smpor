/**
 * Declare layout drawing header
 */


#ifndef DRAW_LAYOUT_H
#define DRAW_LAYOUT_H


extern "C"
{
    #include <stdlib.h>  // for exit()
}


#include "config.h"
#include <vector>
#include <GLFW/glfw3.h>


/* Operation */
#define MAX_WIDTH_INIT_VAL  -10000000

/* Property */
#define WINDOW_WIDTH  800
#define WINDOW_HEIGHT 800

void set_points_attributes();
void draw_vertices(std::vector< std::vector<CoordType> >& coord, \
    std::vector<PartType>& partition);
void draw_edges(std::vector< std::vector<VtxType> >& edges, \
    std::vector< std::vector<CoordType> >& coord);
void draw_layout(std::vector< std::vector<VtxType> >& edges, \
    std::vector< std::vector<CoordType> >& coord,
    std::vector<PartType>& partition);


#endif