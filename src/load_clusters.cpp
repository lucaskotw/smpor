#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "load_clusters.h"




int load_clusters_from_group (char* filePath, std::vector<int>& clusters, int& nCls)
{
    // Steps
    // 1. get file file path
    // 2. Read each line
    //    if line start with "GROUP", skip
    //    else read the node id
    using namespace std;

    stringstream ss;
    // Step 1
    ifstream fin(filePath);

    // Step 2
    string line;
    int vtx;
    int c=0;
    string grp_note = "GROUP";
    std::size_t found;
    while ( getline(fin, line) ) {
        found = line.find(grp_note);
        if (found!=string::npos)
        {
            // cout << "Cluster id#" << c << endl;
            ++c;
        }
        else
        {
            ss.str("");
            ss.clear();
            ss << line;
            ss >> vtx;
            clusters.at(vtx) = c-1;
        }
            
    }
    nCls = c;
    
    return SUCCESS_LOAD_CLUSTERS;
}