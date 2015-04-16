/**
 * Source code used to build the main program
 */
#include <iostream>
#include <ctime> 


#include "load_graph.h"
#include "partition.h"

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
    partition_graph(g, partition, 4);




    /* Ended of the elapsed time measure */
    clock_t end = clock();

    /* Show the elapsed time */    
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    std::cout << "elapsed time = " << elapsed_secs << " secs" << std::endl;

}