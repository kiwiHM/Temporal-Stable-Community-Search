#ifndef TEMPORALGRAPH
#define TEMPORALGRAPH

#include "commonfunctions.h"

class TemporalGraph {
       
    public:

        // pair <int, int> --> to, insertion_time

        // n: the number of vertices; m: the number of edges.
        int n = 0, m = 0;

        // tmax: the maximum time of all temporal edges.
        int tmax = 0;

        // edges[vertex] --> the edges from this vertex.
        // pii = (timestamp, v)
        vector<vector<pii>> edges;

        // temporal_edge[t] --> the edge set at time [t][u].
        vector<unordered_map<int, vector<int>>> temporal_edge;

        // tedges --> the edge set at time at time [t]
        vector<vector<pii>> tedges;

        // queries
        vector<piiii> queries;

        void read_graph(string fp1, string fp2 = "");
        void read_graph_only(string fp);

        vector<int> keynodes(int ts, int te);

        // vector <int> tonodes(int u, int ts, int te){

        //     assert(ts >= 1 && ts <= tmax);
        //     assert(te >= 1 && te <= tmax);

        //     vector <int> ret;
        //     for (int i = ts; i <= te; i++){
        //         for (auto item : temporal_edge[i][u])
        //             ret.push_back
        //     }

        //     return ret;

        // }
};

#endif