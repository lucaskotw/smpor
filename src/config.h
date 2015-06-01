#ifndef CONFIG_H
#define CONFIG_H


/**********
 * OFFSET *
 **********/
#define PORT_BOUNDARY_OFFSET 1

#define NULL_PORT        -1
#define NULL_BOUNDARY_PT -1


/*******************
 * Data Type Alias *
 *******************/
#define VtxType   int     // type of vertex id
#define WgtType   double  // type of edge weight value
#define CoordType double  // type of coordinates value
#define PartType  int     // type of partition value
#define DistType  int     // type of distance value


/**********************
 * Status Assignement *
 **********************/

/* partition stage */
#define SUCCESS_MATCHING       0
#define SUCCESS_COARSENING     0
#define SUCCESS_INIT_PARTITION 0
#define SUCCESS_PARTITION      0
#define SUCCESS_UNCOARSENING   0
#define FAIL_PARTITION        -1

/* smpor stage */
#define SUCCESS_CREATE_SMALL_GRAPH_LIST         0
#define SUCCESS_BFS                   0
#define SUCCESS_CLUSTER_PLACEMENT   0
#define SUCCESS_INTRA_STRESS_MAJORIZATION 0
#define SUCCESS_GET_RADII_CENTERS 0
#define SUCCESS_INTER_STRESS_MAJORIZATION 0
#define SUCCESS_SHIFT_INTRA_CLUSTER       0
#define SUCCESS_SM_PG                 0
#define SUCCESS_SM                    0
#define SUCCESS_SMPOR                 0


/* ports and boundary assignment stage */
#define SUCCESS_PBA  0
#define SUCCESS_CTRL 0
#define SUCCESS_LAYOUT_REFINEMENT 0;

#endif