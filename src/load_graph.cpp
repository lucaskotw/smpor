#include "load_graph.h"


static
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
        if (pWgt == 0) {
            pWgt = 1;
            g.add_edge(u-1, v-1, pWgt);
        } else {
            g.add_edge(u-1, v-1, pWgt);
        }   
    }

    return SUCCESS_CREATE_GRAPH;
    
}


/******************************************************************************
 *     Load the data from Graph Modeling Language to Graph data structure     *
 ******************************************************************************/


int load_graph_from_gml(char* filePath, Graph::Graph& g)
{
    FILE * f;

    // error handling while fopen
    if ((f = fopen(filePath, "r")) == NULL)
    {
        return FAIL_OPEN_FILE;
    }

    std::cout << f << std::endl;

    igraph_t ig;
    igraph_read_graph_graphml(&ig, f, 0);

    fclose(f);

    return SUCCESS_CREATE_GRAPH;
    
}



/******************************************************************************
 *                Load the cluster data from cluster file                     *
 ******************************************************************************/
int load_cluster_from_cls(char* filePath, std::vector<int>& clusters, int& nCls)
{

    // Steps
    // 1. get file file path
    // 2. get the header -> indicate cluster number
    // 3. each line will be the vtxs in cluster, separate by space
    using namespace std;
    stringstream ss;
    // Step 1
    ifstream fin(filePath);

    // Step 2
    string n_cls;
    getline(fin, n_cls);
    cout << n_cls << endl;
    ss << n_cls;
    ss >> nCls;


    // Step 3
    string line;
    int vtx;
    int c=0;
    while ( getline(fin, line) ) {
        if (line.empty())
        {
            continue; // be careful: an empty line might be read
        }
        else
        {
            ss.str("");
            ss.clear();
            ss << line;
            while (ss >> vtx)
            {
                clusters.at(vtx-1) = c;
            }
            ++c;
        }
            
    }

    return SUCCESS_LOAD_CLUSTERS;
    
}