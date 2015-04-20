/**
 * Source code used to build the stress majorization with p-node overlap removal
 * 
 * Expected Drawing Dimension = 2
 */
#include <iostream>
#include <ctime> 


#define PARTITION_NUM 4


#include "load_graph.h"
#include "partition.h"
#include "smpor.h"
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
    g.print_graph();

    std::cout << "number of vertices = " << g.get_num_vtxs() << std::endl;
    /* Partition the Graph */
    std::vector<int> partition;
    partition_graph(g, partition, PARTITION_NUM);


    /* Stress Majorization */
    std::vector< std::vector<CoordType> > coord;

    smpor(g, coord, partition, PARTITION_NUM);


    /* Ended of the elapsed time measure */
    clock_t end = clock();


    /* Draw the Layout*/

    // draw_layout(edges, coord, partition);

    /* Show the elapsed time */    
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    std::cout << "elapsed time = " << elapsed_secs << " secs" << std::endl;

}