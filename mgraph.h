#include "commonfunctions.h"
#include "temporal_graph.h"
#include "ds.h"
#include "Graph.h"


const string core_idx_path = "/new_home/yige/core-idx.txt";
const string core_tree_path = "/new_home/yige/core-tree.txt";
const string core_log_path = "core_log.txt";
const string compressed_index_path = "/new_home/yige/com-idx.txt";
const string log_path = "log.txt";
const string main_output_path = "output.txt";
const string advanced_output_path = "output.txt";


class MFU{
    int n;
    vector <int> p, siz;
public:
    MFU(int n_);
    ~MFU();
    void init(int n_);
    int find(int u);
    bool merge(int u, int v);
};

class mgraph{

private:

    int n_;
    int t_max_;
    int k_max_;
    int min_k_; 
    
    unsigned long long core_idx_size;
    unsigned long long time_part1, time_part2, time_part3;

    vector <vector<pii>>* core_t_;

    TemporalGraph G;

    vector <vector<vector<pii> > >* treeedge;

    // (acttm, u, v)
    vector < vector<piii> >* added_edges;

    // adding time and deleting time
    vector < vector<piii> >* idx_add;
    vector < vector<piii> >* idx_del; 

public:

    mgraph();

    void load_idx(string idx_path = core_idx_path);
    void build_coreness_index(string graph_path, string idx_path = core_idx_path, string log_path = core_log_path);

    // query for the time, where the given node u is engaged into k-core.
    int query_coretime(int u, int k, int ts);

    void load_graph(string fp1, string fp2 = "");

    // ASF-index
    int build_index_tree_time_anchored(int k, int t_st);
    void build_index_trees();
    void build_index_trees_pruned();

    void output_index_trees();

    void load_index_trees(string idx_path);
    void save_index_trees(string idx_path);

    piii query_time_anchored(int t_st, int t_ed, int u, int k);
    piii query_time_LR_fixed(int L, int R, int t_ed, int u, int k);
    piii query(int ts, int te, int u, int k);
    piii query_plain(int ts, int te, int u, int k);
    vector <bool> query_cs(int ts, int te, int u, int k);
    vector <bool> query_cs_cckcore(int ts, int te, int u, int k);
    void handle_queries(string output_filepath);
    void handle_queries_cs(string output_filepath);
    void handle_queries_basic(string output_filepath);
    void handle_queries_kcore();
    void handle_queries_connected_kcore();

    // AG-index part
    void perpare_for_edges(); // record the label-changed edges at each ts. 
    void build_compressed_index_old(); // The old AG-build function without optimized
    void build_compressed_index(); // The main building function

    void load_compressed_index(string idx_path);
    void save_compressed_index(string idx_path);

    piii query_time_anchored_idxgiven(int ts, int te, int u, int k, 
            const vector< unordered_set<pii, pair_hash> > &idx_tr);
    piii query_time_anchored_idxgiven_opt(int t_st, int t_ed, int u, int k, 
                const vector <unordered_set<pii, pair_hash> > &idx_tr);
    piii query_time_LR_fixed_compressed(int L, int R, int t_ed, int u, int k,
                const vector <unordered_set<pii, pair_hash> > &idx_tr);
    piii query_compressed(int ts, int te, int u, int k);
    // piii query_compressed_opt(int ts, int te, int u, int k);
    piii query_compressed_basic(int ts, int te, int u, int k);

    vector <int> cs_search(int ts, int tr, int u, int k,
        const vector <unordered_set<pii, pair_hash> > &idx_tr);
    vector<int> query_cs_compressed(int ts, int te, int u, int k, piii &ret_info);

    void handle_queries_compressed(string output_filepath);
    void handle_queries_compressed_cs(string output_filepath);
    void handle_queries_compressed_basic(string output_filepath);

    long long hash_set(const set<int> &S);
    int query_distinct_core_time_anchored(int ts, int te, int u, int k, 
                const vector <unordered_set<pii, pair_hash> > &idx_tr,
                const vector <long long> &randval,
                unordered_set <long long> &kcore);
    double query_distinct_core(int ts, int te, int u, int k);
    void handle_queries_distinct_core(string output_filepath);

    void specific_index(int k, int ts);

    // For expriment
    void output_topnodes(int topk);

    void debug();
};
