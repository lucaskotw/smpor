/**
 * Source code used to build the stress majorization with p-node overlap removal
 * 
 * Expected Drawing Dimension = 2
 */
#include <iostream>
#include <ctime> 
#include <fstream>
#include <string>
#include <regex>
#include <time.h>


#define CLUSTER_NUM 3
#define DRAWING_DIM   2


#include "load_graph.h"
#include "load_clusters.h"
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
    g.print_graph();

    // create distance matrix
    DenseMat dist_mat(g.get_num_vtxs(), g.get_num_vtxs());
    distance_matrix(g, dist_mat);
    cout << "dist matrix created" << std::endl;

    /* Clusters the Graph */

    // NP04 sample
    std::vector<int> clusters(g.get_num_vtxs());
    int n_cls;
    load_clusters_from_group(argv[2], clusters, n_cls);
    cout << "cluster size = " << n_cls << endl;
    std::cout << "cluster" << std::endl;
    for (std::vector<int>::iterator it=clusters.begin();\
    it!=clusters.end();
    ++it)
    {
        std::cout << *it << " ";
    }
    std::cout << std::endl;


    /* SMPOR */
    std::vector< std::vector<CoordType> > coord(g.get_num_vtxs(),\
                                                std::vector<CoordType>(DRAWING_DIM));

    std::vector< std::vector<CoordType> > center_coord(n_cls);

    std::vector< WgtType > radii(n_cls);
    double interpolation = atof(argv[3]);

    smpor(g, g.get_num_vtxs(), dist_mat, coord, center_coord, radii, clusters, n_cls, interpolation);
    std::cout << "smpor finish" << std::endl;

    // Original Stress Majorization

    // std::cout << "coord" << std::endl;
    // for (std::vector< std::vector<CoordType> >::iterator it1=coord.begin();\
    //     it1!=coord.end();
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

    // layout_refinement(g, dist_mat, clusters, n_cls, center_coord, radii, coord);


    /* Ended of the elapsed time measure */
    clock_t end = clock();


    /* Show the elapsed time */    
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    std::cout << "elapsed time = " << elapsed_secs << " secs" << std::endl;

    /* Output the information and coordinates */
    char curr_time[100];
    time_t t;
    struct tm * timeinfo;
    time (&t);
    timeinfo = localtime (&t);
    strftime(curr_time, 100, "%Y%m%d%H%M%S", timeinfo);
    regex rgx(".*/(\\w+)/.*");
    smatch match;
    string data(argv[1]);
    regex_search(data, match, rgx);
    string outfile = string(curr_time) + "_" + string(match[1]) + ".out";
    string outdir = "output/";
    cout << outfile << endl;

    fstream fo;
    fo.open(outdir+outfile, fstream::app);
    strftime(curr_time, 100, "%Y%m%d%H%M%S", timeinfo);
    fo << "elapsed_secs:" << elapsed_secs << endl;
    fo << "interpolation:" << interpolation << endl;
    fo << "coordination:" << endl;
    for (int i=0; i<coord.size(); ++i)
    {
        fo << coord.at(i).at(0) << " " << coord.at(i).at(1) << endl;
    }
    fo.close();


    /* Draw the Layout*/

    draw_layout(edges, coord, clusters, ports, boundary_pts,\
        ports_coords, boundary_pts_coords, ctrl_pts_coords,\
        center_coord, radii);

    return 0;

}