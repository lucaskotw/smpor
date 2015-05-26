/**
 * Source code used to build the stress majorization with p-node overlap removal
 * 
 * Expected Drawing Dimension = 2
 */
#include <iostream>
#include <ctime> 


#define CLUSTER_NUM 4
#define DRAWING_DIM   2


#include "load_graph.h"
#include "distance.h"
#include "smpor.h"
#include "pba.h"
#include "draw_layout.h"

int main(int argc, char** argv)
{
    /* Start Timer */
    using namespace std;
    clock_t begin = clock();




    /* Code to measure the elapsed time */
    /* Load the Graph*/
    Graph g(0);
    load_graph_from_mm(argv[1], g);
    // g.print_graph();

    // create distance matrix
    DenseMat dist_mat(g.get_num_vtxs(), g.get_num_vtxs());
    distance_matrix(g, dist_mat);
    std::cout << dist_mat << std::endl;
    // std::cout << dist_mat.row(0) << std::endl;
    // std::cout << dist_mat.rows() << std::endl;
    // std::cout << dist_mat.cols() << std::endl;
    // // check pointer
    // double* dist_pt;
    // dist_pt = dist_mat.data();
    // double** dist_ppt;
    // dist_ppt[0] = *dist_pt;
    // std::cout << dist_pt << std::endl;
    // std::cout << dist_ppt << std::endl;


    


    // std::cout << *mat_pt << std::endl;
    // std::cout << dist_mat.col(0) << std::endl;
    // std::cout << dist_mat.col(1) << std::endl;


    // // run hierarchical clustering
    // double** data;
    // int** mask;
    // double weight [];
    // int transpose = 0;
    // char dist;
    // char method = 'm';
    // treecluster(dist_mat.rows(), dist_mat.cols(), data, mask,\
    //             weight, transpose, dist, method, dist_pt);

    // std::cout << "number of vertices = " << g.get_num_vtxs() << std::endl;
    // /* Clusters the Graph */
    // std::vector<int> clusters;
    // clusters_graph(g, clusters, clusters_NUM);

    std::vector<int> clusters;
    for (int i=0; i<5;i++)
    {
        clusters.push_back(0);
    }
    for (int i=0; i<7;i++)
    {
        clusters.push_back(1);
    }
    for (int i=0; i<4;i++)
    {
        clusters.push_back(2);
    }
    



    /* Stress Majorization */
    std::vector< std::vector<CoordType> > coord(g.get_num_vtxs(),\
                                                std::vector<CoordType>(DRAWING_DIM));

    std::vector< std::vector<CoordType> > center_coord(3); // 3 stands for
                                                           // clusters number
    std::vector< WgtType > radii(3);
    smpor(g.get_num_vtxs(), dist_mat, coord, center_coord, radii, clusters, 3);

    std::cout << "smpor coord" << std::endl;
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

    std::cout << "center coord" << std::endl;
    for (std::vector< std::vector<CoordType> >::iterator it1=center_coord.begin();\
        it1!=center_coord.end();
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

    std::cout << "radius" << std::endl;
    for (std::vector<CoordType>::iterator it=radii.begin();\
    it!=radii.end();
    ++it)
    {
        std::cout << *it << " ";
    }
    std::cout << std::endl;

    // /* port and boundary assignment - create edges */
    // std::vector< std::vector<VtxType> > edges;
    // std::vector< std::vector<VtxType> > ports;
    // std::vector< std::vector<VtxType> > boundary_pts;
    // std::vector< std::vector<CoordType> > ports_coords;
    // std::vector< std::vector<CoordType> > boundary_pts_coords;
    // std::vector< std::vector<CoordType> > ctrl_pts_coords;
    // port_and_boundary_assignment(g, pg, clusters, coord, \
    //     edges, ports, boundary_pts, ports_coords, boundary_pts_coords,\
    //     ctrl_pts_coords);
    

    // // create center and radius
    // std::vector< std::vector<CoordType> > centers;
    // std::vector<WgtType> radius;
    // std::vector<CoordType> center_pair(2);
    // for (int i=0; i<pg.get_num_vtxs(); ++i)
    // {
    //     center_pair.at(0) = coord.at(pg.get_center_id(i)).at(0);
    //     center_pair.at(1) = coord.at(pg.get_center_id(i)).at(1);
    //     centers.push_back(center_pair);
    //     radius.push_back(pg.get_radius(i));
    // }


    /* Ended of the elapsed time measure */
    clock_t end = clock();


    /* Draw the Layout*/

    // draw_layout(edges, coord, clusters, ports, boundary_pts,\
    //     ports_coords, boundary_pts_coords, ctrl_pts_coords,\
    //     centers, radius);

    /* Show the elapsed time */    
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    std::cout << "elapsed time = " << elapsed_secs << " secs" << std::endl;

    return 0;

}