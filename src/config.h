#ifndef CONFIG_H
#define CONFIG_H


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
#define SUCCESS_CREATE_PGRAPH  0
#define SUCCESS_BFS            0
#define SUCCESS_SMPOR          0

#endif