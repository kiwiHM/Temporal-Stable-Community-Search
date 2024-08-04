#include "temporal_graph.h"
using namespace std;

namespace Community_search{

    TemporalGraph G;

    void init(string fp1, string fp2 = ""){
        G.read_graph(fp1, fp2);
    }

    // return the edgeset remained after reduction into k-core.
    vector <pii> core_decomposition(int L, int R, int k){
        // cout << "core_decomp in ts = " << L << ", te = " << R << ", k = " << k << endl;
        vector <pii> ret;
        vector <int> keys = G.keynodes(L, R);
        unordered_map <int, bool> dlted;
        unordered_map <int, int> deg;
        unordered_map <int, vector<int> > tonode;
        // undirected edges remained.
        unordered_set <pii, pair_hash> edgeset; 
        vector < unordered_set<int> > tax(G.n);
        for (int t = L; t <= R; t++)
            for (auto e : G.tedges[t]){
                pii ne = make_pair(min(e.first, e.second), max(e.first, e.second));
                edgeset.insert(ne);
            }
        for (auto e : edgeset){
            int u = e.first, v = e.second;
            deg[u]++, deg[v]++;
            tonode[u].push_back(v);
            tonode[v].push_back(u);
        }
        int dmax = 0;
        for (auto u : keys){
            assert(0 <= deg[u] && deg[u] < G.n);
            tax[deg[u]].insert(u);
            dmax = max(dmax, deg[u]);
        }

        // the queue records the nodes being deleted.
        queue <int> que;
        // inque markes if the node has been added to the queue.
        unordered_map <int, bool> inque;
        for (int cur = 0; cur < k; cur++){
            // the corenesses of the remained nodes are \ge cur;
            for (auto u : tax[cur])
                if (inque[u] == false)
                    que.push(u), inque[u] = true;
            while (que.size()){
                int u = que.front(); que.pop();
                dlted[u] = true;
                for (auto v : tonode[u]){
                    tax[deg[v]].erase(v);
                    deg[v]--;
                    tax[deg[v]].insert(v);
                    if (inque[v] == false && deg[v] <= cur)
                        que.push(v), inque[v] = true;
                }
            }
        }
        for (auto e : edgeset)
            if ((dlted[e.first] || dlted[e.second]) == false)
                ret.push_back(e);
        return ret;
    }

    vector <int> find_component(int u, vector <pii> edgeset){
        vector <int> ret;
        unordered_map <int, bool> vis;
        unordered_map <int, vector<int> > tonode;
        for (auto e : edgeset){
            int u = e.first, v = e.second;
            tonode[u].push_back(v);
            tonode[v].push_back(u);
        }
        if (tonode[u].size() == 0){
            // in this case, coreness[u] < K
            return ret;
        }
        queue <int> que;
        que.push(u), vis[u] = true, ret.push_back(u);
        while (que.size()){
            int u = que.front(); que.pop();
            for (auto v : tonode[u]){
                if (vis.find(v) == vis.end()){
                    vis[v] = true;
                    que.push(v);
                    ret.push_back(v);
                }
            }
        }
        return ret;
    }

    vector <int> community_query(int L, int R, int u, int k){
        vector <pii> edgeset = core_decomposition(L, R, k);
        // cout << "Finish core decomposition" << endl;
        return find_component(u, edgeset);
    }

    // return (lasting_time, (left_timepoint, right_timepoint))
    // if two have the same lasting_time, we try minimizing L, then R.
    piii query(int ts, int te, int u, int k){
        int lst_tm = 0;
        pii ret_itv = make_pair(-1, -1);
        for (int L = ts; L <= te; L++){
            for (int R = L; R <= te; R++){
                vector <int> com = community_query(L, R, u, k);
                if (com.size() == 0) continue;
                // print_vec(com, "The initial com is: ", L == 6 && R == 8);
                // Try using binary search to find the lasting time.
                int lft = R, rgt = te, ret = -1;
                while (lft <= rgt){
                    int mid = (lft + rgt) >> 1;
                    vector <int> ncom = community_query(L, mid, u, k);
                    // print_vec(ncom, "The update com is: ", L == 6 && R == 8);
                    if (checksame_vector(com, ncom))
                        lft = mid + 1, ret = mid;
                    else rgt = mid - 1;
                }

                int cur_tm = ret == -1 ? 0 : ret - R + 1;
                if (cur_tm > lst_tm){
                    lst_tm = cur_tm;
                    ret_itv = make_pair(L, R);
                }
                // printf("ts = %d, te = %d, te' = %d\n", L, R, ret);
            }
        }
        return make_pair(lst_tm, ret_itv);
    }

    void handle_queries(string output_filepath){
        cout << "Handling queries..." << endl;
        unsigned long long start_time_whole = currentTime();
        ofstream outputFile(output_filepath);
        assert(outputFile.is_open());
        for (auto q : G.queries){
            int ts = q.first.first, te = q.first.second;
            int u = q.second.first, k = q.second.second;
            auto ans = query(ts, te, u, k);
            // printf("The community: L = %d, R = %d, with lasting time = %d\n", 
            //         ans.second.first, ans.second.second, ans.first);
            outputFile << "Query: " << ts << ' ' << te << ' ' << u << ' ' << k << endl; 
            outputFile << "  lasting time = " << ans.first << "  ";
            outputFile << "Com_ts = " << ans.second.first << ", Com_te = " << ans.second.second << endl;
        }
        unsigned long long end_time_whole = currentTime();
        cout << "Time Cost in online-query " << timeFormatting((end_time_whole - start_time_whole) / 1e6).str() << endl;
        outputFile << "Time Cost in online-query " << timeFormatting((end_time_whole - start_time_whole) / 1e6).str() << endl;
        outputFile.close();
    }
};

int main(int argc, char *argv[]){
    if (argc > 2) Community_search::init(argv[1], argv[2]);
    else Community_search::init(argv[1]);
    Community_search::handle_queries(argv[3]);
    return 0;
}