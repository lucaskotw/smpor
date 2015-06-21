#ifndef LOAD_CLUSTERS_H
#define LOAD_CLUSTERS_H

#include <vector>


/* global var for */
#define FILE_BUFFER_SIZE       256
#define SUCCESS_LOAD_CLUSTERS  0


int load_clusters_from_group (char* filePath, std::vector<int>& clusters, int& nCls);


#endif