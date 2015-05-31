#ifndef LOAD_GRAPH_H
#define LOAD_GRAPH_H

extern "C"
{
    #include <stdio.h>
}

#include <fstream>
#include "graph.h"

#include <igraph/igraph.h>


/* move to main.cpp */
#define FAIL_OPEN_FILE         2
#define FAIL_READ_FILE_HEADER  3
#define FAIL_COUNT_SPACE       4
#define SUCCESS_CREATE_GRAPH   0
#define SUCCESS_READ_EDGE      0

/* global var for */
#define FILE_BUFFER_SIZE       256


int load_graph_from_mm (char* filePath, Graph::Graph& g);
int load_graph_from_gml(char* filePath, Graph::Graph& g);


#endif