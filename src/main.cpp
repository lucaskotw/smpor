/**
 * Source code used to build the stress majorization with p-node overlap removal
 * 
 * Expected Drawing Dimension = 2
 */
#include <iostream>
#include <ctime> 


#define CLUSTER_NUM 3
#define DRAWING_DIM   2


#include "load_graph.h"
#include "distance.h"
#include "sm.h"
#include "smpor.h"
#include "lr.h"
#include "draw_layout.h"

int main(int argc, char** argv)
{
    /* Start Timer */
    using namespace std;
    clock_t begin = clock();




    /* Code to measure the elapsed time */
    /* Load the Graph*/
    Graph g(0);
    // load_graph_from_gml(argv[1], g);
    load_graph_from_mm(argv[1], g);
    // g.print_graph();

    // create distance matrix
    DenseMat dist_mat(g.get_num_vtxs(), g.get_num_vtxs());
    distance_matrix(g, dist_mat);
    std::cout << dist_mat << std::endl;

    // /* Clusters the Graph */

    // NP04 sample
    std::vector<int> clusters;
    for (int i=0; i<6;i++)
    {
        clusters.push_back(0);
    }
    for (int i=0; i<7;i++)
    {
        clusters.push_back(1);
    }
    for (int i=0; i<8;i++)
    {
        clusters.push_back(2);
    }


    /* SMPOR */
    std::vector< std::vector<CoordType> > coord(g.get_num_vtxs(),\
                                                std::vector<CoordType>(DRAWING_DIM));

    std::vector< std::vector<CoordType> > center_coord(CLUSTER_NUM);

    std::vector< WgtType > radii(CLUSTER_NUM);
    smpor(g, g.get_num_vtxs(), dist_mat, coord, center_coord, radii, clusters, 3);

    // Onriginal Stress Majorization

    std::cout << "coord" << std::endl;
    for (std::vector< std::vector<CoordType> >::iterator it1=coord.begin();\
        it1!=coord.end();
        ++it1)
    {
        for (std::vector<CoordType>::iterator it2=(*it1).begin();\
        it2!=(*it1).end();
        ++it2)
        {
            std::cout << *it2 << " ";
        }
        std::cout << std::endl;
    }

    // std::cout << "center coord" << std::endl;
    // for (std::vector< std::vector<CoordType> >::iterator it1=center_coord.begin();\
    //     it1!=center_coord.end();
    //     ++it1)
    // {
    //     for (std::vector<CoordType>::iterator it2=(*it1).begin();\
    //     it2!=(*it1).end();
    //     ++it2)
    //     {
    //         std::cout << *it2 << " ";
    //     }
    //     std::cout << std::endl;
    // }

    // std::cout << "radius" << std::endl;
    // for (std::vector<CoordType>::iterator it=radii.begin();\
    // it!=radii.end();
    // ++it)
    // {
    //     std::cout << *it << " ";
    // }
    // std::cout << std::endl;

    // /* port and boundary assignment - create edges */
    std::vector< std::vector<VtxType> > edges;
    std::vector< std::vector<VtxType> > ports;
    std::vector< std::vector<VtxType> > boundary_pts;
    std::vector< std::vector<CoordType> > ports_coords;
    std::vector< std::vector<CoordType> > boundary_pts_coords;
    std::vector< std::vector<CoordType> > ctrl_pts_coords;
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

    // // port_and_boundary_assignment(g, pg, clusters, coord, \
    // //     edges, ports, boundary_pts, ports_coords, boundary_pts_coords,\
    // //     ctrl_pts_coords);
    layout_refinement(g, dist_mat, clusters, CLUSTER_NUM, coord);


    /* Ended of the elapsed time measure */
    clock_t end = clock();


    /* Draw the Layout*/

    draw_layout(edges, coord, clusters, ports, boundary_pts,\
        ports_coords, boundary_pts_coords, ctrl_pts_coords,\
        center_coord, radii);

    /* Show the elapsed time */    
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    std::cout << "elapsed time = " << elapsed_secs << " secs" << std::endl;

    return 0;

}