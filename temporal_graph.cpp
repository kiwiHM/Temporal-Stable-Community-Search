#include "temporal_graph.h"

void TemporalGraph::read_graph_only(string fp){
    puts("Start reading the graph...");
    FILE *File_graph = freopen(fp.c_str(), "r", stdin);
    
    read(n), read(m), read(tmax);
    printf("n = %d, m = %d, tmax = %d\n", n, m, tmax);

    edges.resize(n + 1); // indexes of nodes ranging in [1, n]
    temporal_edge.resize(tmax + 1);
    tedges.resize(tmax + 1);

    int cnt_ignored = 0;
    for (int i = 1; i <= m; i++){
        int u, v, t; read(u), read(v), read(t);
        assert(0 <= u && u < n);
        assert(0 <= v && v < n);
        assert(0 <= t && t < tmax);
        if (u == v){
            cnt_ignored++;
            continue;
        }
        edges[u].push_back(make_pair(t, v));
        edges[v].push_back(make_pair(t, u));
        temporal_edge[t][u].push_back(v);
        temporal_edge[t][v].push_back(u);
        tedges[t].push_back(make_pair(u, v));
    }
    m -= cnt_ignored;
    puts("Finish reading the graph...");
}

void TemporalGraph::read_graph(string fp1, string fp2){
    puts("Start reading the graph...");
    FILE *File_graph = freopen(fp1.c_str(), "r", stdin);
    
    read(n), read(m), read(tmax);
    printf("n = %d, m = %d, tmax = %d\n", n, m, tmax);

    edges.resize(n + 1); // indexes of nodes ranging in [1, n]
    temporal_edge.resize(tmax + 1);
    tedges.resize(tmax + 1);

    int cnt_ignored = 0;
    for (int i = 1; i <= m; i++){
        int u, v, t; read(u), read(v), read(t);
        assert(0 <= u && u < n);
        assert(0 <= v && v < n);
        assert(0 <= t && t < tmax);
        if (u == v){
            cnt_ignored++;
            continue;
        }
        edges[u].push_back(make_pair(t, v));
        edges[v].push_back(make_pair(t, u));
        temporal_edge[t][u].push_back(v);
        temporal_edge[t][v].push_back(u);
        tedges[t].push_back(make_pair(u, v));
    }
    m -= cnt_ignored;

    if (fp2 != ""){
        FILE *File_query = freopen(fp2.c_str(), "r", stdin);    
        puts("Start reading the queries...");
        puts(fp2.c_str());
    }

    queries.clear();
    int q; read(q);
    printf("number of queries = %d\n", q);
    for (int i = 1; i <= q; i++){
        int ts, te, u, k;
        read(ts), assert(0 <= ts && ts < tmax);
        read(te), assert(ts <= te && te < tmax);
        read(u), assert(0 <= u && u < n);
        read(k), assert(1 <= k && k < n);
        queries.push_back(make_pair(
                    make_pair(ts, te), make_pair(u, k)));
    }
    puts("Finish reading the graph...");
}

vector <int> TemporalGraph::keynodes(int ts, int te){
    unordered_map <int, bool> apr;
    vector <int> ret;
    for (int i = ts; i <= te; i++){
        for (auto e : tedges[i]){
            int u = e.first, v = e.second;
            if (apr.find(u) == apr.end())
                ret.push_back(u);
            if (apr.find(v) == apr.end())
                ret.push_back(v);
        }
    }
    return ret;
}