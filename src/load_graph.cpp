#include "load_graph.h"


int read_edge_from_input(char * buff, int buffSize, \
    VtxType & u, VtxType & v, WgtType & pWgt)
{
    pWgt = 0;      // initialize prefer edge weight
    int space = 0; // calculate the space in buffer
    for (int i=0; i<buffSize; ++i) {
        if (buff[i] == '\n') { // check the end of buffer
            break;
        }
        if (buff[i] == ' ') { // record if a space appear
            ++space;
        }
    }

    // Deal with the space condition
    if (space == 1)
    {
        sscanf(buff, "%d %d", &u, &v);
    } 
    else if (space == 2)
    {
        sscanf(buff, "%d %d %lg", &u, &v, &pWgt);
    }
    else
    {
        return FAIL_COUNT_SPACE;
    }

    return SUCCESS_READ_EDGE;
}


int load_graph_from_mm(char* filePath, Graph::Graph& g)
{
    // FILE * f;
    // MM_typecode matcode; // this var will record the type of the matrix

    // error handling while fopen
    // if ((f = fopen(filePath, "r")) == NULL)
    // {
    //     return FAIL_OPEN_FILE;
    // }

    // read in the type of the matrix
    // if (mm_read_banner(f, &matcode) != 0)
    // {
    //     return FAIL_READ_FILE_HEADER;
    // }

    std::ifstream fin(filePath);

    // Ignore headers and comments:
    while (fin.peek() == '%') fin.ignore(2048, '\n');

    // get the size and basic information of the mm file
    int rr;
    int cc;
    int nZero;
    fin >> rr >> cc >> nZero;
    // if ((ret_code = mm_read_mtx_crd_size(f, &Rows, &Cols, &nZero)) !=0)
    // {
    //     exit(1);
    // }

    // initialize graph g
    g = Graph(rr);

    // read the matrix content
    VtxType u;
    VtxType v;
    WgtType pWgt;
    char buff[FILE_BUFFER_SIZE]; // buffer to read each line
    fin.getline(buff, FILE_BUFFER_SIZE); // (CHECK) why this be null?
    while (fin.getline(buff, FILE_BUFFER_SIZE))
    {
        read_edge_from_input(buff, FILE_BUFFER_SIZE, u, v, pWgt);
        std::cout << buff << std::endl;
        std::cout << u << std::endl;
        std::cout << v << std::endl;
        std::cout << pWgt << std::endl;
        if (pWgt == 0) {
            pWgt = 1;
            g.add_edge(u-1, v-1, pWgt);
        } else {
            g.add_edge(u-1, v-1, pWgt);
        }   
    }

    return SUCCESS_CREATE_GRAPH;
    
}