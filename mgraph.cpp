#include "mgraph.h"

MFU::MFU(int n_){
    n = n_;
    p.resize(n);
    siz.resize(n);
    for (int i = 0; i < n; i++)
        p[i] = i, siz[i] = 1;
}
MFU::~MFU(){
    vector <int>().swap(p);
    vector <int>().swap(siz);
}
int MFU::find(int u){
    if (u == this->p[u]) return u;
    return this->p[u] = this->find(this->p[u]);
}
bool MFU::merge(int u, int v){
    int fu = this->find(u);
    int fv = this->find(v);
    if (fu == fv) return false;
    if (siz[fu] > siz[fv]){
        this->p[fv] = fu;
        this->siz[fu] += this->siz[fv];
    } else {
        this->p[fu] = fv;
        this->siz[fv] += this->siz[fu];
    }
    return true;
}


mgraph::mgraph(){
    min_k_ = 2;
    core_t_ = nullptr;
}

void mgraph::load_idx(string idx_path) {
    auto fp = fopen(idx_path.c_str(),"rb");
    fread(&n_,sizeof(unsigned int),1,fp);
    k_max_ = 0;
    t_max_ = 0;

    if (core_t_ == nullptr) core_t_ = new vector<vector<pii>>[n_];

    core_idx_size = 0;
    for (int u = 0; u < n_; ++u) {
        int cs;
        fread(&cs,sizeof(int),1,fp);
        core_idx_size += sizeof(int);
        core_t_[u].resize(cs);
        if(k_max_ < cs-1) k_max_ = cs-1;

        for (int k = min_k_; k < cs; ++k) {
            int ccs;
            fread(&ccs,sizeof(int),1,fp);
            core_idx_size += sizeof(int);

            for (int i = 0; i < ccs; ++i) {
                int a,b;
                fread(&a,sizeof(int),1,fp);
                fread(&b,sizeof(int),1,fp);
                core_idx_size += sizeof(int) * 2;
                core_t_[u][k].emplace_back(make_pair(a,b));
            }
            if (t_max_ < core_t_[u][k].back().second) t_max_ = core_t_[u][k].back().second;
        }
    }
    fclose(fp);
}

void mgraph::build_coreness_index(string graph_path, string idx_path, string log_path){
    auto *g = new Graph();
    g->init_log(log_path);
    g->load(graph_path, true);
    g->index();
    g->write_idx(idx_path);
    load_idx(idx_path);
    cout << "Finish building the Core-index." << endl;

    // for (int i = 0; i < n_; i++){
    //     cout << endl << "Node = " << i << ' ' << "Coresize = " << core_t_[i].size() << endl;
    //     for (int k = 1; k <= k_max_; k++){
    //         if (k >= core_t_[i].size()) continue;

    //         cout << "K = " << k << "   siz = " << core_t_[i][k].size() << endl;
    //         for (int j = 0; j < core_t_[i][k].size(); j++){
    //             cout << "[" << core_t_[i][k][j].first << ", ";
    //             cout << core_t_[i][k][j].second << "], ";
    //         }
    //         cout << endl;
    //     }
    // }
}

int mgraph::query_coretime(int u, int k, int ts){
    // printf("Query_coretime : u = %d, k = %d, ts = %d\n", u, k, ts); 
    // printf("   core_t_[%u].size() = %u\n", u, core_t_[u].size());
    assert(0 <= ts && ts < G.tmax);
    if (k == 1) return ts;
    if (k >= core_t_[u].size())
        return G.tmax + 1;
    assert(core_t_[u][k].size() != 0);
    pii tmp = make_pair(ts, G.tmax + 5); 
    // Find the first it, where it.first greater than ts
    auto it = upper_bound(core_t_[u][k].begin(), core_t_[u][k].end(), tmp);
    // if (it == core_t_[u][k].end())
    //     return G.tmax + 1;
    assert(it != core_t_[u][k].begin());
    if (it != core_t_[u][k].begin())
        --it;
    int Its = it->first, Ite = it->second;
    // cout << "      !!! " << Its << ' ' << Ite << endl;
    assert(Its <= ts); // TODO
    if (Ite >= G.tmax)
        return G.tmax + 1;
    return Ite;
}

void mgraph::load_graph(string fp1, string fp2){
    G.read_graph(fp1, fp2);

    // cout << "Tmax = " << t_max_ << endl;
    // for (int k = 1; k <= k_max_; k++){
    //     cout << "\n k = " << k << endl;
    //     for (int i = 0; i < n_; i++){
    //         cout << endl << "Node = " << i << endl;
    //         // cout << "K = " << k << "   siz = " << core_t_[i][k].size() << endl;
    //         for (int ts = 0; ts < t_max_; ts++)
    //             cout << query_coretime(i, k, ts) << ' '; 
    //         cout << endl;
    //     }
    // }
}

int mgraph::build_index_tree_time_anchored(int k, int t_st){
    // cout << "   On k = " << k << ", start_time = " << t_st << endl;
    treeedge[k][t_st].resize(G.n); // Adj-linked-list
    vector <piii> edges;
    
    unsigned long long start_time = currentTime();
    for (int t = t_st; t < G.tmax; t++){
        for (pii e : G.tedges[t]){
            int tsp = t;
            tsp = max(tsp, query_coretime(e.first, k, t_st));
            tsp = max(tsp, query_coretime(e.second, k, t_st));
            assert(t_st <= tsp);
            if (tsp < G.tmax)
                edges.push_back(make_pair(tsp, e));
        }
    }
    unsigned long long time_point1 = currentTime();
    time_part1 += time_point1 - start_time;

    sort(edges.begin(), edges.end());
    unsigned long long time_point2 = currentTime();
    time_part2 += time_point2 - time_point1;
    
    MFU c(G.n + 1);
    int cnt = 0;
    for (auto e : edges){
        int t = e.first, u = e.second.first, v = e.second.second;
        if (c.merge(u, v)){
            cnt++;
            treeedge[k][t_st][u].push_back(make_pair(t, v));
            treeedge[k][t_st][v].push_back(make_pair(t, u));
        }
        if (cnt == G.n - 1) 
        break;
    }
    unsigned long long time_point3 = currentTime();
    time_part3 += time_point3 - time_point2;
    return edges.size();
}

void mgraph::build_index_trees(){
    cout << "Start building the ASF-index" << endl;
    treeedge = new vector<vector<vector<pii>>>[k_max_ + 1];
    unsigned long long start_time_whole = currentTime();
    for (int k = 1; k <= k_max_; k++){
        unsigned long long start_time = currentTime();
        treeedge[k].resize(G.tmax);
        long long sum_edges = 0;
        time_part1 = time_part2 = time_part3 = 0;
        for (int t = 0; t < G.tmax; t++)
            sum_edges += build_index_tree_time_anchored(k, t);
        unsigned long long end_time = currentTime();
        cout << "Time in Iteration k = " << k << ": " << timeFormatting((end_time - start_time) / 1e6).str() << endl;
        // cout << "Average edges.size() = " << ((float) sum_edges / G.tmax) << ' ';
        // cout << "Average time_part {1, 2, 3} = {";
        // cout << timeFormatting(time_part1 / 1e6).str() << ", ";
        // cout << timeFormatting(time_part2 / 1e6).str() << ", ";
        // cout << timeFormatting(time_part3 / 1e6).str() << "}." << endl;
    }
    unsigned long long end_time_whole = currentTime();
    cout << "Finish building the ASF-index" << endl;
    ofstream logFile(log_path.c_str(), ios::app);
    logFile << "Time Cost in the ASF-index Construction " << timeFormatting((end_time_whole - start_time_whole) / 1e6).str() << endl;
}

void mgraph::build_index_trees_pruned(){
    cout << "Start building the ASF-index (optimized)" << endl;
    // not utilized.
    cout << "Finish building the ASF-index" << endl;
}

void mgraph::output_index_trees(){
    for (int k = 1; k <= k_max_; k++){
        for (int t = 0; t < G.tmax; t++){
            cout << "k = " << k << ", t_s = " << t << endl;
            for (int i = 0; i < G.n; i++){
                for (auto e : treeedge[k][t][i])
                    if (e.second > i)
                        cout << "{" << i << ' ' << e.second << ' ' << e.first << "}, ";
            }
            cout << endl;
        }
    }
}

void mgraph::load_index_trees(string idx_path){
    cout << "Start loading the ASF-index" << endl;
    auto fp = fopen(idx_path.c_str(), "rb");
    treeedge = new vector<vector<vector<pii>>>[k_max_ + 1];
    for (int k = 1; k <= k_max_; k++){
        treeedge[k].resize(G.tmax);
        for (int t = 0; t < G.tmax; t++){
            treeedge[k][t].resize(G.n); 
            int cnt_e;
            fread(&cnt_e, sizeof(int), 1, fp);
            for (int i = 0; i < cnt_e; i++){
                int u, v, timestamp;
                fread(&u, sizeof(int), 1, fp);
                fread(&v, sizeof(int), 1, fp);
                fread(&timestamp, sizeof(int), 1, fp);
                treeedge[k][t][u].push_back(make_pair(timestamp, v));
                treeedge[k][t][v].push_back(make_pair(timestamp, u));
            }
        }
    }
    cout << "Finish loading the ASF-index" << endl;
}

void mgraph::save_index_trees(string idx_path){
    cout << "Start saving the ASF-index" << endl;
    auto fp = fopen(idx_path.c_str(), "wb");
    unsigned long long idx_size = 0;
    for (int k = 1; k <= k_max_; k++){
        for (int t = 0; t < G.tmax; t++){
            int cnt_e = 0;
            for (int i = 0; i < G.n; i++){
                for (auto e : treeedge[k][t][i])
                    if (e.second > i) cnt_e++;
            }
            fwrite(&cnt_e, sizeof(int), 1, fp);
            idx_size += sizeof(int);
            for (int i = 0; i < G.n; i++){
                for (auto e : treeedge[k][t][i])
                    if (e.second > i){
                        fwrite(&i, sizeof(int), 1, fp);
                        fwrite(&e.second, sizeof(int), 1, fp);
                        fwrite(&e.first, sizeof(int), 1, fp);
                        idx_size += sizeof(int) * 3;
                    }
            }
        }
    }
    cout << "Core-index size: " << ((double) core_idx_size / 1024 / 1024) << "MB" << endl;
    cout << "ASF-index size: " << ((double) idx_size / 1024 / 1024) << "MB" << endl;
    cout << "Total index size: " << ((double) (core_idx_size + idx_size) / 1024 / 1024) << "MB" << endl; 
    cout << "Finish saving the ASF-index" << endl;
}

piii mgraph::query_time_anchored(int t_st, int t_ed, int u, int k){
    // Do BFS to get udpate the spanning-tree
    int core_t_u = query_coretime(u, k, t_st);
    if (core_t_u >= G.tmax || core_t_u > t_ed)
        return make_pair(0, make_pair(-1, -1));
    queue <int> que;
    vector <int> ctm(G.n);
    vector <bool> vis(G.n);
    for (int i = 0; i < G.n; i++) 
        ctm[i] = vis[i] = 0;
    que.push(u), ctm[u] = core_t_u, vis[u] = true;
    while (que.size()){
        int cu = que.front(); que.pop();
        for (auto e_out : treeedge[k][t_st][cu]){
            int t = e_out.first, v = e_out.second;
            assert(0 <= t && t < G.tmax);
            assert(0 <= v && v < G.n);
            if (vis[v]) continue; // visited
            que.push(v), ctm[v] = max(ctm[cu], t), vis[v] = true;
        }
    }
    vector <int> vec;
    for (int i = 0; i < G.n; i++)
        if (vis[i] && ctm[i] <= t_ed) vec.push_back(ctm[i]);
    vec.push_back(t_ed + 1);
    sort(vec.begin(), vec.end());
    // print_vec(vec, "timing_vector :");
    int lst_tm = 0, R = -1;
    for (int i = 1, si = vec.size(); i < si; i++)
        if (vec[i] - vec[i - 1] > lst_tm){
            lst_tm = vec[i] - vec[i - 1];
            R = vec[i - 1];
        }
    return make_pair(lst_tm, make_pair(t_st, R));
}


piii mgraph::query(int ts, int te, int u, int k){
    piii ret = make_pair(0, make_pair(-1, -1));
    for (int L = ts; L <= te; L++){
        piii tmp = query_time_anchored(L, te, u, k);
        if (tmp.first > ret.first)
            ret = tmp;
    }
    return ret;
}

vector <bool> mgraph::query_cs(int ts, int te, int u, int k){
    vector <bool> ret(G.n);
    for (int i = 0; i < G.n; i++)
        ret[i] = 0;
    piii com_info = query(ts, te, u, k);
    if (com_info.first == 0) 
        return ret;
    int tl = com_info.second.first, tr = com_info.second.second;
    queue <int> que;
    que.push(u), ret[u] = true;
    while (que.size()){
        int cu = que.front(); que.pop();
        for (auto e_out : treeedge[k][tl][cu]){
            int t = e_out.first, v = e_out.second;
            if (ret[v]) continue;
            if (t <= tr) que.push(v), ret[v] = true;
        }
    }
    return ret;
}

vector <bool> mgraph::query_cs_cckcore(int ts, int te, int u, int k){
    vector <bool> ret(G.n);
    for (int i = 0; i < G.n; i++)
        ret[i] = 0;
    if (query_coretime(u, k, ts) >= G.tmax)
        return ret;
    queue <int> que;
    que.push(u), ret[u] = true;
    while (que.size()){
        int cu = que.front(); que.pop();
        for (auto e_out : treeedge[k][ts][cu]){
            int t = e_out.first, v = e_out.second;
            if (ret[v] == false)
                que.push(v), ret[v] = true;
        }
    }
    return ret;
}

piii mgraph::query_time_LR_fixed(int L, int R, int t_ed, int u, int k){
    queue <int> que;
    vector <bool> vis(G.n);
    for (int i = 0; i < G.n; i++)
        vis[i] = false;
    que.push(u), vis[u] = true;
    int upd_tm = t_ed + 1;
    while (que.size()){
        int cu = que.front(); que.pop();
        for (auto e_out : treeedge[k][L][cu]){
            int t = e_out.first, v = e_out.second;
            assert(0 <= t && t < G.tmax);
            assert(0 <= v && v < G.n);
            if (vis[v]) continue; // visited
            if (t <= R) que.push(v), vis[v] = true;
            else upd_tm = min(upd_tm, t);
        }
    }
    return make_pair(upd_tm - R, make_pair(L, R));
}

piii mgraph::query_plain(int ts, int te, int u, int k){
    piii ret = make_pair(0, make_pair(-1, -1));
    for (int L = ts; L <= te; L++){
        for (int R = query_coretime(u, k, L); R <= te; R++){
            piii tmp = query_time_LR_fixed(L, R, te, u, k);
            if (tmp.first > ret.first)
                ret = tmp;
        }
    }
    return ret;
}

void mgraph::handle_queries(string output_filepath){
    cout << "Solving by AITE on the ASF-index" << endl;
    unsigned long long start_time = currentTime();
    ofstream outputFile(output_filepath);
    assert(outputFile.is_open());
    int Q = G.queries.size(), Q_5 = Q / 5, cnt = 0, percent = 0;
    for (auto q : G.queries){
        int ts = q.first.first, te = q.first.second;
        int u = q.second.first, k = q.second.second;
        auto ans = this->query(ts, te, u, k);
        // printf("The community: L = %d, R = %d, with lasting time = %d\n", 
        //         ans.second.first, ans.second.second, ans.first);
        outputFile << "Query: " << ts << ' ' << te << ' ' << u << ' ' << k << endl; 
        outputFile << "  lasting time = " << ans.first << "  ";
        outputFile << "Com_ts = " << ans.second.first << ", Com_te = " << ans.second.second << endl;

        cnt++;
        if (Q_5 && cnt % Q_5 == 0){
            percent += 20;
            cout << "Finished " << percent << "% queries." << endl;
        }
    }
    unsigned long long end_time = currentTime();
    ofstream logFile(log_path.c_str(), ios::app);
    logFile << "Time used in queries: " << timeFormatting((end_time - start_time) / 1e6).str() << endl;
    cout << "Time used in queries: " << timeFormatting((end_time - start_time) / 1e6).str() << endl;
    outputFile << "Time used in queries: " << timeFormatting((end_time - start_time) / 1e6).str() << endl;
    outputFile.close();
}

void mgraph::handle_queries_cs(string output_filepath){
    cout << "Solving by community search" << endl;
    unsigned long long start_time = currentTime();
    ofstream outputFile(output_filepath);
    assert(outputFile.is_open());
    int Q = G.queries.size(), Q_20 = Q / 20, cnt = 0, percent = 0;
    for (auto q : G.queries){
        int ts = q.first.first, te = q.first.second;
        int u = q.second.first, k = q.second.second;
        auto ans = this->query_cs(ts, te, u, k);
        // The first line is the query node
        // The second line is the vector Y_pre.
        if (query_coretime(u, k, ts) < G.tmax){
            outputFile << u << endl;
            for (int i = 0; i < G.n; i++)
                outputFile << int(ans[i]) << (i == G.n - 1 ? '\n' : ' ');
        }
        cnt++;
        if (Q_20 && cnt % Q_20 == 0){
            percent += 5;
            cout << "Finished " << percent << "% queries." << endl;
        }
    }
    outputFile.close();
    unsigned long long end_time = currentTime();
    cout << "Time used in queries: " << timeFormatting((end_time - start_time) / 1e6).str() << endl;
}

void mgraph::handle_queries_kcore(){
    cout << "Solving by kcore()" << endl;
    unsigned long long start_time = currentTime();
    ofstream outputFile("main_kcore_output.txt");
    assert(outputFile.is_open());
    int Q = G.queries.size(), Q_20 = Q / 20, cnt = 0, percent = 0;
    for (auto q : G.queries){
        int ts = q.first.first, te = q.first.second;
        int u = q.second.first, k = q.second.second;
        // The first line is the query node
        // The second line is the vector Y_pre.
        if (query_coretime(u, k, ts) < G.tmax){
            outputFile << u << endl;
            for (int i = 0; i < G.n; i++)
                outputFile << (query_coretime(i, k, ts) < G.tmax) << (i == G.n - 1 ? '\n' : ' ');
        }
        cnt++;
        if (Q_20 && cnt % Q_20 == 0){
            percent += 5;
            cout << "Finished " << percent << "% queries." << endl;
        }
    }
    outputFile.close();
    unsigned long long end_time = currentTime();
    cout << "Time used in queries: " << timeFormatting((end_time - start_time) / 1e6).str() << endl;
}


void mgraph::handle_queries_basic(string output_filepath){
    cout << "Solving by BITE on ASF-index" << endl;
    unsigned long long start_time = currentTime();
    ofstream outputFile(output_filepath);
    assert(outputFile.is_open());
    int Q = G.queries.size(), Q_20 = Q / 20, cnt = 0, percent = 0;
    for (auto q : G.queries){
        int ts = q.first.first, te = q.first.second;
        int u = q.second.first, k = q.second.second;
        auto ans = this->query_plain(ts, te, u, k);
        // printf("The community: L = %d, R = %d, with lasting time = %d\n", 
        //         ans.second.first, ans.second.second, ans.first);
        outputFile << "Query: " << ts << ' ' << te << ' ' << u << ' ' << k << endl; 
        outputFile << "  lasting time = " << ans.first << "  ";
        outputFile << "Com_ts = " << ans.second.first << ", Com_te = " << ans.second.second << endl;

        cnt++;
        if (Q_20 && cnt % Q_20 == 0){
            percent += 5;
            cout << "Finished " << percent << "% queries." << endl;
        }
    }
    unsigned long long end_time = currentTime();
    ofstream logFile(log_path.c_str(), ios::app);
    logFile << "Time used in queries: " << timeFormatting((end_time - start_time) / 1e6).str() << endl;
    cout << "Time used in queries: " << timeFormatting((end_time - start_time) / 1e6).str() << endl;
    outputFile << "Time used in queries: " << timeFormatting((end_time - start_time) / 1e6).str() << endl;
    outputFile.close();
}

void mgraph::handle_queries_connected_kcore(){
    cout << "Solving by cc-kcore query()" << endl;
    unsigned long long start_time = currentTime();
    ofstream outputFile("main_cckcore_output.txt");
    assert(outputFile.is_open());
    int Q = G.queries.size(), Q_20 = Q / 20, cnt = 0, percent = 0;
    for (auto q : G.queries){
        int ts = q.first.first, te = q.first.second;
        int u = q.second.first, k = q.second.second;
        auto ans = this->query_cs_cckcore(ts, te, u, k);
        // printf("The community: L = %d, R = %d, with lasting time = %d\n", 
        //         ans.second.first, ans.second.second, ans.first);
        if (query_coretime(u, k, ts) < G.tmax){
            outputFile << u << endl;
            for (int i = 0; i < G.n; i++)
                outputFile << int(ans[i]) << (i == G.n - 1 ? '\n' : ' ');
        }
        cnt++;
        if (Q_20 && cnt % Q_20 == 0){
            percent += 5;
            cout << "Finished " << percent << "% queries." << endl;
        }
    }
    unsigned long long end_time = currentTime();
    ofstream logFile(log_path.c_str(), ios::app);
    logFile << "Time used in queries: " << timeFormatting((end_time - start_time) / 1e6).str() << endl;
    cout << "Time used in queries: " << timeFormatting((end_time - start_time) / 1e6).str() << endl;
    outputFile << "Time used in queries: " << timeFormatting((end_time - start_time) / 1e6).str() << endl;
    outputFile.close();
}

// For AG-index

void mgraph::perpare_for_edges(){
    cout << "Start labeling the edges" << endl;
    unsigned long long start_time_whole = currentTime();
    added_edges = new vector< vector<piii> > [k_max_ + 1];
    for (int k = 2; k <= k_max_; k++){
        cout << "Labeling on k = " << k << endl;
        // if (k != 12) continue; // TODO
        added_edges[k].resize(G.tmax);
        // record (left_value, node, id)
        priority_queue <piii, vector<piii>, less<piii> > que;
        vector <int> coretm(G.n), upd_nodes;
        for (int i = 0; i < G.n; i++)
            coretm[i] = G.tmax + 5; // inf
        for (int i = 0; i < G.n; i++){
            if (k >= core_t_[i].size() || core_t_[i][k].size() == 0)
                continue;
            int bt = core_t_[i][k].size() - 1;
            // cout << "Queue inserted: " << core_t_[i][k][bt].first << ' ' << i << endl;
            que.push(make_pair(core_t_[i][k][bt].first, make_pair(i, bt)));
            coretm[i] = core_t_[i][k][bt].second;
        }
        unordered_map <pii, int, pair_hash> cur_acttm;
        for (int curtm = G.tmax - 1; curtm >= 0; --curtm){
            // cout << "time = " << curtm << endl;
            // for (int i = 0; i < G.n; i++){
            //     cout << coretm[i] << (i == G.n - 1 ? '\n' : ' ');
            // }
            // add edges with timestamp = current_time.
            for (auto e : G.tedges[curtm]){
                int u = e.first, v = e.second;
                if (coretm[u] >= G.tmax || coretm[v] >= G.tmax)
                    continue;
                int acttm = max(curtm, max(coretm[u], coretm[v]));
                added_edges[k][curtm].push_back(make_pair(acttm, make_pair(u, v)));
                cur_acttm[make_pair(u, v)] = acttm;
            }
            // mark the nodes whose coretime changed.
            while (que.size() && que.top().first > curtm){
                auto tmp = que.top(); que.pop();
                int ltm = tmp.first, u = tmp.second.first, p = tmp.second.second;
                // cout << "Core updated" << ' ' << u << endl;
                // cout << "Next que top: " << que.top().first << ' ' << curtm << ' ' << que.top().second.first << endl;
                assert(p != 0);
                upd_nodes.push_back(u);
                coretm[u] = core_t_[u][k][p - 1].second;
                que.push(make_pair(core_t_[u][k][p - 1].first, make_pair(u, p - 1)));
            }
            // cout << "Out of while: " << que.top().first << ' ' << curtm << ' ' << que.top().second.first << endl;
            // if (upd_nodes.size()) for (int i = 0; i < upd_nodes.size(); i++) cout << "upd" << upd_nodes[i] << "\n";
            // found all edges whose active time changed.
            for (auto u : upd_nodes){
                for (int i = G.edges[u].size() - 1; i >= 0; i--){
                    int tsp = G.edges[u][i].first, v = G.edges[u][i].second;
                    // if (u == X && curtm > 64){
                        // cout << "Curtime = " << curtm << endl;
                        // cout << "    !!! " << tsp << ' ' << v << " coretime: " << coretm[u] << ' ' << coretm[v] << endl;
                        // cout << "    current act time = " << cur_acttm[make_pair(u, v)] << endl;
                    // }
                    if (tsp < curtm) 
                        break;
                    int nacttm = max(tsp, max(coretm[u], coretm[v]));
                    if (nacttm >= G.tmax) continue;
                    pii e = make_pair(u, v);
                    if (cur_acttm.find(e) == cur_acttm.end() || nacttm < cur_acttm[e]){
                        // active time of e changes.
                        added_edges[k][curtm].push_back(make_pair(nacttm, e));
                        cur_acttm[e] = nacttm;
                    }
                }
            }
            upd_nodes.clear();
        }
    }
    cout << "Finish labeling the edges" << endl;
    unsigned long long end_time_whole = currentTime();
    cout << "Time Cost in labeling " << timeFormatting((end_time_whole - start_time_whole) / 1e6).str() << endl;
}

void mgraph::build_compressed_index_old(){
    perpare_for_edges();
    cout << "Start building AG-index" << endl;
    unsigned long long start_time_whole = currentTime();
    idx_add = new vector < vector<piii> > [k_max_ + 1];
    idx_del = new vector < vector<piii> > [k_max_ + 1];
    for (int k = 2; k <= k_max_; k++){
        idx_add[k].resize(G.tmax);
        idx_del[k].resize(G.tmax);
        // cout << "The current k = " << k << endl;
        // if (k != 12) continue; // TODO
        unsigned long long start_time = currentTime();
        Minimum_Spanning_Tree mst(G.n);
        for (int ts = G.tmax - 1; ts >= 0; ts--){
            for (auto e : added_edges[k][ts]){
                // if (e.second.first == X || e.second.second == X) 
                    // cout << "Try: " << e.first << ' ' << e.second.first << ' ' << e.second.second << " at time = " << ts << endl;
                piii de = make_pair(-1, make_pair(0, 0));
                if (mst.insert(e, de)){
                    // if (e.second.first == X || e.second.second == X){
                        // cout << "Insert: " << e.first << ' ' << e.second.first << ' ' << e.second.second << " at time = " << ts << endl;
                        // cout << "   k = " << k << ' ' << "ts = " << ts << endl;
                    // }
                    idx_add[k][ts].push_back(e);
                    if (de.first != -1){
                        idx_del[k][ts].push_back(de);
                        // cout << "Delete: " << de.first << ' ' << de.second.first << ' ' << de.second.second << endl;
                    }
                }
            }
        }
        unsigned long long end_time = currentTime();
        cout << "Time Cost in k = " << k << ", " << timeFormatting((end_time - start_time) / 1e6).str() << endl;
    }
    cout << "finish building AG-index" << endl;
    unsigned long long end_time_whole = currentTime();
    cout << "Time Cost in AG-index Construction " << timeFormatting((end_time_whole - start_time_whole) / 1e6).str() << endl;
}

void mgraph::build_compressed_index(){
    cout << "Start building AG-index optimized" << endl;
    unsigned long long start_time_whole = currentTime();
    idx_add = new vector < vector<piii> > [k_max_ + 1];
    idx_del = new vector < vector<piii> > [k_max_ + 1];
    long double tot_c = 0;
    for (int k = 2; k <= k_max_; k++){
        idx_add[k].resize(G.tmax);
        idx_del[k].resize(G.tmax);
        // cout << "The current k = " << k << endl;
        // if (k != 12) continue; // TODO
        unsigned long long start_time = currentTime();
        priority_queue <piii, vector<piii>, less<piii> > que;
        vector <int> coretm(G.n), upd_nodes;
        for (int i = 0; i < G.n; i++)
            coretm[i] = G.tmax + 5; // inf
        for (int i = 0; i < G.n; i++){
            if (k >= core_t_[i].size() || core_t_[i][k].size() == 0)
                continue;
            int bt = core_t_[i][k].size() - 1;
            // cout << "Queue inserted: " << core_t_[i][k][bt].first << ' ' << i << endl;
            que.push(make_pair(core_t_[i][k][bt].first, make_pair(i, bt)));
            coretm[i] = core_t_[i][k][bt].second;
        }        
        unordered_map <pii, int, pair_hash> cur_acttm;
        Minimum_Spanning_Tree mst(G.n);
        int tot_edges = 0;
        for (int ts = G.tmax - 1; ts >= 0; ts--){
            vector <piii> tau;
            for (auto e : G.tedges[ts]){
                int u = e.first, v = e.second;
                if (coretm[u] >= G.tmax || coretm[v] >= G.tmax)
                    continue;
                int acttm = max(ts, max(coretm[u], coretm[v]));
                tau.push_back(make_pair(acttm, make_pair(u, v)));
                cur_acttm[make_pair(u, v)] = acttm;
            }
            // mark the nodes whose coretime changed.
            while (que.size() && que.top().first > ts){
                auto tmp = que.top(); que.pop();
                int ltm = tmp.first, u = tmp.second.first, p = tmp.second.second;
                // cout << "Core updated" << ' ' << u << endl;
                // cout << "Next que top: " << que.top().first << ' ' << curtm << ' ' << que.top().second.first << endl;
                assert(p != 0);
                upd_nodes.push_back(u);
                coretm[u] = core_t_[u][k][p - 1].second;
                que.push(make_pair(core_t_[u][k][p - 1].first, make_pair(u, p - 1)));
            }
            // cout << "Out of while: " << que.top().first << ' ' << curtm << ' ' << que.top().second.first << endl;
            // if (upd_nodes.size()) for (int i = 0; i < upd_nodes.size(); i++) cout << "upd" << upd_nodes[i] << "\n";
            // found all edges whose active time changed.
            for (auto u : upd_nodes){
                for (int i = G.edges[u].size() - 1; i >= 0; i--){
                    int tsp = G.edges[u][i].first, v = G.edges[u][i].second;
                    // if (u == X && curtm > 64){
                        // cout << "Curtime = " << curtm << endl;
                        // cout << "    !!! " << tsp << ' ' << v << " coretime: " << coretm[u] << ' ' << coretm[v] << endl;
                        // cout << "    current act time = " << cur_acttm[make_pair(u, v)] << endl;
                    // }
                    if (tsp < ts) 
                        break;
                    int nacttm = max(tsp, max(coretm[u], coretm[v]));
                    if (nacttm >= G.tmax) continue;
                    pii e = make_pair(u, v);
                    if (cur_acttm.find(e) == cur_acttm.end() || nacttm < cur_acttm[e]){
                        // active time of e changes.
                        tau.push_back(make_pair(nacttm, e));
                        cur_acttm[e] = nacttm;
                    }
                }
            }
            upd_nodes.clear();

            tot_edges += tau.size();
            for (auto e : tau){
                // if (e.second.first == X || e.second.second == X) 
                    // cout << "Try: " << e.first << ' ' << e.second.first << ' ' << e.second.second << " at time = " << ts << endl;
                piii de = make_pair(-1, make_pair(0, 0));
                if (mst.insert(e, de)){
                    // if (e.second.first == X || e.second.second == X){
                    // cout << "Insert: " << e.first << ' ' << e.second.first << ' ' << e.second.second << " at time = " << ts << endl;
                        // cout << "   k = " << k << ' ' << "ts = " << ts << endl;
                    // }
                    idx_add[k][ts].push_back(e);
                    if (de.first != -1){
                        idx_del[k][ts].push_back(de);
                        // cout << "Delete: " << de.first << ' ' << de.second.first << ' ' << de.second.second << endl;
                    }
                }
            }
        }
        tot_c += (long double) tot_edges / G.m;
        unsigned long long end_time = currentTime();
        cout << "Time Cost in k = " << k << ", " << timeFormatting((end_time - start_time) / 1e6).str() << endl;
    }
    tot_c /= (k_max_ - 1);
    cout << "Finish building AG-index" << endl;
    cout << "The c_bar = " << tot_c << endl;
    unsigned long long end_time_whole = currentTime();
    cout << "Time Cost in AG-index Construction " << timeFormatting((end_time_whole - start_time_whole) / 1e6).str() << endl;
    ofstream logFile(log_path.c_str(), ios::app);
    logFile << "Time Cost in the AG-index Construction " << timeFormatting((end_time_whole - start_time_whole) / 1e6).str() << endl;
}

piii mgraph::query_compressed(int ts, int te, int u, int k){
    piii ret = make_pair(0, make_pair(-1, -1));
    if (k > k_max_) return ret;
    vector < unordered_set<pii, pair_hash> > idx_tr(G.n);
    for (int curtm = G.tmax - 1; curtm >= ts; curtm--){
        // cout << curtm << ' ' << ts << ' ' << u << ' ' << k << endl;
        // cout << "Current time = " << curtm << endl;
        for (auto e : idx_add[k][curtm]){
            int w = e.first, u = e.second.first, v = e.second.second;
            // cout << "Insert: " << w << ' '  << u << ' ' << v << endl;
            idx_tr[u].insert(make_pair(w, v));
            idx_tr[v].insert(make_pair(w, u));
        }
        for (auto e : idx_del[k][curtm]){
            int w = e.first, u = e.second.first, v = e.second.second;
            // cout << "Erase: " << w << ' ' << u << ' ' << v << endl;
            idx_tr[u].erase(make_pair(w, v));
            idx_tr[v].erase(make_pair(w, u));
        }
        if (curtm <= te){
            piii tmp = query_time_anchored_idxgiven_opt(curtm, te, u, k, idx_tr);
            if (tmp.first >= ret.first)
                ret = tmp;
        }
    }
    for (int i = 0; i < G.n; i++)
        idx_tr[i].clear();
    return ret;
}

piii mgraph::query_time_anchored_idxgiven_opt(int t_st, int t_ed, int u, int k, 
                const vector< unordered_set<pii, pair_hash> > &idx_tr){
    // cout << "In query" << endl;
    // Do BFS to get udpate the spanning-tree
    int core_t_u = query_coretime(u, k, t_st);
    if (core_t_u >= G.tmax || core_t_u > t_ed){
        // cout << "Out of query" << endl;
        return make_pair(0, make_pair(-1, -1));
    }
    queue <int> que;
    vector <int> ctm(G.n);
    vector <bool> vis(G.n);
    vector <int> vislist;
    for (int i = 0; i < G.n; i++) 
        ctm[i] = vis[i] = 0;
    que.push(u), ctm[u] = core_t_u, vis[u] = true, vislist.push_back(u);
    while (que.size()){
        int cu = que.front(); que.pop();
        for (auto e_out : idx_tr[cu]){
            int t = e_out.first, v = e_out.second;
            assert(0 <= t && t < G.tmax);
            assert(0 <= v && v < G.n);
            // assert(t >= query_coretime(v, k, t_st)); // TODO
            if (vis[v]) continue; // visited
            que.push(v), ctm[v] = max(ctm[cu], t), vis[v] = true, vislist.push_back(v);
        }
    }
    // cout << "ts = " << t_st << endl;
    // for (int i = 0; i < G.n; i++)
    //     cout << ctm[i] << (i == G.n - 1 ? '\n' : ' ');
    vector <int> vec;
    if (G.tmax <= G.n){
        vector <int> tax(G.tmax + 20);
        for (auto i : vislist)
            if (vis[i] && ctm[i] <= t_ed)
                tax[ctm[i]]++;
        tax[t_ed + 1]++;
        for (int i = 0; i < G.tmax + 20; i++)
            for (int j = 0; j < tax[i]; j++)
                vec.push_back(i);
    } else {
        for (int i = 0; i < G.n; i++)
            if (vis[i] && ctm[i] <= t_ed) vec.push_back(ctm[i]);
        vec.push_back(t_ed + 1);
        sort(vec.begin(), vec.end());
    }
    // print_vec(vec, "timing_vector :");
    int lst_tm = 0, R = -1;
    for (int i = 1, si = vec.size(); i < si; i++)
        if (vec[i] - vec[i - 1] > lst_tm){
            lst_tm = vec[i] - vec[i - 1];
            R = vec[i - 1];
        }
    // cout << "Out of query" << endl;
    return make_pair(lst_tm, make_pair(t_st, R));
}

piii mgraph::query_time_anchored_idxgiven(int t_st, int t_ed, int u, int k, 
                const vector< unordered_set<pii, pair_hash> > &idx_tr){
    // cout << "In query" << endl;
    // Do BFS to get udpate the spanning-tree
    int core_t_u = query_coretime(u, k, t_st);
    if (core_t_u >= G.tmax || core_t_u > t_ed){
        // cout << "Out of query" << endl;
        return make_pair(0, make_pair(-1, -1));
    }
    queue <int> que;
    vector <int> ctm(G.n);
    vector <bool> vis(G.n);
    for (int i = 0; i < G.n; i++) 
        ctm[i] = vis[i] = 0;
    que.push(u), ctm[u] = core_t_u, vis[u] = true;
    while (que.size()){
        int cu = que.front(); que.pop();
        for (auto e_out : idx_tr[cu]){
            int t = e_out.first, v = e_out.second;
            assert(0 <= t && t < G.tmax);
            assert(0 <= v && v < G.n);
            // assert(t >= query_coretime(v, k, t_st)); // TODO
            if (vis[v]) continue; // visited
            que.push(v), ctm[v] = max(ctm[cu], t), vis[v] = true;
        }
    }
    // cout << "ts = " << t_st << endl;
    // for (int i = 0; i < G.n; i++)
    //     cout << ctm[i] << (i == G.n - 1 ? '\n' : ' ');
    vector <int> vec;
    for (int i = 0; i < G.n; i++)
        if (vis[i] && ctm[i] <= t_ed) vec.push_back(ctm[i]);
    vec.push_back(t_ed + 1);
    sort(vec.begin(), vec.end());
    // print_vec(vec, "timing_vector :");
    int lst_tm = 0, R = -1;
    for (int i = 1, si = vec.size(); i < si; i++)
        if (vec[i] - vec[i - 1] > lst_tm){
            lst_tm = vec[i] - vec[i - 1];
            R = vec[i - 1];
        }
    // cout << "Out of query" << endl;
    return make_pair(lst_tm, make_pair(t_st, R));
}

void mgraph::handle_queries_compressed(string output_filepath){
    cout << "Solving by AITE on AG-index" << endl;
    unsigned long long start_time = currentTime();
    ofstream outputFile(output_filepath);
    assert(outputFile.is_open());
    int Q = G.queries.size(), Q_5 = Q / 5, cnt = 0, percent = 0;
    for (auto q : G.queries){
        int ts = q.first.first, te = q.first.second;
        int u = q.second.first, k = q.second.second;
        auto ans = this->query_compressed(ts, te, u, k);
        // printf("The community: L = %d, R = %d, with lasting time = %d\n", 
        //         ans.second.first, ans.second.second, ans.first);
        outputFile << "Query: " << ts << ' ' << te << ' ' << u << ' ' << k << endl; 
        outputFile << "  lasting time = " << ans.first << "  ";
        outputFile << "Com_ts = " << ans.second.first << ", Com_te = " << ans.second.second << endl;

        cnt++;
        if (Q_5 && cnt % Q_5 == 0){
            percent += 20;
            cout << "Finished " << percent << "% queries." << endl;
        }
    }
    unsigned long long end_time = currentTime();
    ofstream logFile(log_path.c_str(), ios::app);
    logFile << "Time used in queries: " << timeFormatting((end_time - start_time) / 1e6).str() << endl;
    cout << "Time used in queries: " << timeFormatting((end_time - start_time) / 1e6).str() << endl;
    outputFile << "Time used in queries: " << timeFormatting((end_time - start_time) / 1e6).str() << endl;
    outputFile.close();
}

piii mgraph::query_time_LR_fixed_compressed(int L, int R, int t_ed, int u, int k,
                const vector< unordered_set<pii, pair_hash> > &idx_tr){
    queue <int> que;
    vector <bool> vis(G.n);
    for (int i = 0; i < G.n; i++)
        vis[i] = false;
    que.push(u), vis[u] = true;
    int upd_tm = t_ed + 1;
    while (que.size()){
        int cu = que.front(); que.pop();
        for (auto e_out : idx_tr[cu]){
            int t = e_out.first, v = e_out.second;
            assert(0 <= t && t < G.tmax);
            assert(0 <= v && v < G.n);
            if (vis[v]) continue; // visited
            if (t <= R) que.push(v), vis[v] = true;
            else upd_tm = min(upd_tm, t);
        }
    }
    return make_pair(upd_tm - R, make_pair(L, R));
}

piii mgraph::query_compressed_basic(int ts, int te, int u, int k){
    piii ret = make_pair(0, make_pair(-1, -1));
    if (k > k_max_) return ret;
    vector < unordered_set<pii, pair_hash> > idx_tr(G.n);
    for (int curtm = G.tmax - 1; curtm >= ts; curtm--){
        // cout << curtm << ' ' << ts << ' ' << u << ' ' << k << endl;
        // cout << "Current time = " << curtm << endl;
        for (auto e : idx_add[k][curtm]){
            int w = e.first, u = e.second.first, v = e.second.second;
            // cout << "Insert: " << w << ' '  << u << ' ' << v << endl;
            idx_tr[u].insert(make_pair(w, v));
            idx_tr[v].insert(make_pair(w, u));
        }
        for (auto e : idx_del[k][curtm]){
            int w = e.first, u = e.second.first, v = e.second.second;
            // cout << "Erase: " << w << ' ' << u << ' ' << v << endl;
            idx_tr[u].erase(make_pair(w, v));
            idx_tr[v].erase(make_pair(w, u));
        }
        if (curtm <= te){
            for (int R = te, stR = query_coretime(u, k, curtm); R >= stR; R--){
                piii tmp = query_time_LR_fixed_compressed(curtm, R, te, u, k, idx_tr);
                if (tmp.first >= ret.first)
                    ret = tmp;
            }
        }
    }
    for (int i = 0; i < G.n; i++)
        idx_tr[i].clear();
    return ret;
}

void mgraph::handle_queries_compressed_basic(string output_filepath){
    cout << "Solving by basic BITE on AG-index" << endl;
    unsigned long long start_time = currentTime();
    ofstream outputFile(output_filepath);
    assert(outputFile.is_open());
    int Q = G.queries.size(), Q_5 = Q / 5, cnt = 0, percent = 0;
    for (auto q : G.queries){
        int ts = q.first.first, te = q.first.second;
        int u = q.second.first, k = q.second.second;
        auto ans = this->query_compressed_basic(ts, te, u, k);
        // printf("The community: L = %d, R = %d, with lasting time = %d\n", 
        //         ans.second.first, ans.second.second, ans.first);
        outputFile << "Query: " << ts << ' ' << te << ' ' << u << ' ' << k << endl; 
        outputFile << "  lasting time = " << ans.first << "  ";
        outputFile << "Com_ts = " << ans.second.first << ", Com_te = " << ans.second.second << endl;

        cnt++;
        if (Q_5 && cnt % Q_5 == 0){
            percent += 20;
            cout << "Finished " << percent << "% queries." << endl;
        }
    }
    unsigned long long end_time = currentTime();
    ofstream logFile(log_path.c_str(), ios::app);
    logFile << "Time used in queries: " << timeFormatting((end_time - start_time) / 1e6).str() << endl;
    cout << "Time used in queries: " << timeFormatting((end_time - start_time) / 1e6).str() << endl;
    outputFile << "Time used in queries: " << timeFormatting((end_time - start_time) / 1e6).str() << endl;
    outputFile.close();
}

void mgraph::save_compressed_index(string idx_path){
    cout << "Start saving AG-index" << endl;
    // cout << "M = " << G.m << " " << "Tmax = " << G.tmax << endl;
    auto fp = fopen(idx_path.c_str(), "wb");
    unsigned long long idx_size = 0;
    long double mbar = 0;
    for (int k = 2; k <= k_max_; k++){
        for (int t = 0; t < G.tmax; t++){
            int cnt_e = idx_add[k][t].size();
            fwrite(&cnt_e, sizeof(int), 1, fp);
            idx_size += (cnt_e * 3 + 1) * sizeof(int);
            mbar += cnt_e;
            for (auto e : idx_add[k][t]){
                fwrite(&e.first, sizeof(int), 1, fp);
                fwrite(&e.second.first, sizeof(int), 1, fp);
                fwrite(&e.second.second, sizeof(int), 1, fp);
            }
            cnt_e = idx_del[k][t].size();
            idx_size += (cnt_e * 3 + 1) * sizeof(int);
            fwrite(&cnt_e, sizeof(int), 1, fp);
            for (auto e : idx_del[k][t]){
                fwrite(&e.first, sizeof(int), 1, fp);
                fwrite(&e.second.first, sizeof(int), 1, fp);
                fwrite(&e.second.second, sizeof(int), 1, fp);
            }
        }
    }
    mbar /= k_max_;
    cout << "M_bar = " << mbar << ", M = " << G.m << ", N = " << G.n << ", tmax = " << G.tmax << endl;
    cout << "Core-index size: " << ((double) core_idx_size / 1024 / 1024) << "MB" << endl;
    cout << "AG-index size: " << ((double) idx_size / 1024 / 1024) << "MB" << endl;
    cout << "Total index size: " << ((double) (core_idx_size + idx_size) / 1024 / 1024) << "MB" << endl; 
    cout << "Finish saving AG-index" << endl;
}

void mgraph::load_compressed_index(string idx_path){
    cout << "Start loading AG-index" << endl;
    // cout << "M = " << G.m << " " << "Tmax = " << G.tmax << endl;    
    auto fp = fopen(idx_path.c_str(), "rb");
    unsigned long long idx_size = 0;
    long double mbar = 0;
    idx_add = new vector < vector<piii> > [k_max_ + 1];
    idx_del = new vector < vector<piii> > [k_max_ + 1];
    for (int k = 2; k <= k_max_; k++){
        idx_add[k].resize(G.tmax);
        idx_del[k].resize(G.tmax);
        for (int t = 0; t < G.tmax; t++){
            int cnt_e;
            fread(&cnt_e, sizeof(int), 1, fp);
            idx_size += (cnt_e * 3 + 1) * sizeof(int);
            mbar += cnt_e;
            for (int i = 0; i < cnt_e; i++){
                int ts, u, v;
                fread(&ts, sizeof(int), 1, fp);
                fread(&u, sizeof(int), 1, fp);
                fread(&v, sizeof(int), 1, fp);
                idx_add[k][t].push_back(make_pair(ts, make_pair(u, v)));
            }
            fread(&cnt_e, sizeof(int), 1, fp);
            idx_size += (cnt_e * 3 + 1) * sizeof(int);
            for (int i = 0; i < cnt_e; i++){
                int ts, u, v;
                fread(&ts, sizeof(int), 1, fp);
                fread(&u, sizeof(int), 1, fp);
                fread(&v, sizeof(int), 1, fp);
                idx_del[k][t].push_back(make_pair(ts, make_pair(u, v)));
            }
        }
    }
    mbar /= k_max_;
    // cout << "M_bar = " << mbar << ", M = " << G.m << ", N = " << G.n << ", tmax = " << G.tmax << endl;
    // cout << "Core index size: " << ((double) core_idx_size / 1024 / 1024) << "MB" << endl;
    // cout << "Compressed index size: " << ((double) idx_size / 1024 / 1024) << "MB" << endl;
    // cout << "Total index size = " << ((double) (core_idx_size + idx_size) / 1024 / 1024) << "MB" << endl; 
    cout << "Finish loading AG-index" << endl;
}

vector <int> mgraph::cs_search(int ts, int tr, int u, int k,
        const vector< unordered_set<pii, pair_hash> > &idx_tr){
    vector <int> ret, vis(G.n);
    for (int i = 0; i < G.n; i++)
        vis[i] = 0;
    queue <int> que;
    que.push(u), vis[u] = true, ret.push_back(u);
    while (que.size()){
        int cu = que.front(); que.pop();
        for (auto e_out : idx_tr[cu]){
            int t = e_out.first, v = e_out.second;
            if (vis[v]) continue;
            if (t <= tr) que.push(v), vis[v] = true, ret.push_back(v);
        }
    }
    return ret;
}

vector<int> mgraph::query_cs_compressed(int ts, int te, int u, int k, piii &ret_info){
    // cout << "  ! inquery: " << ts << ' ' << te << ' ' << u << ' ' << k << endl;
    ret_info = make_pair(0, make_pair(-1, -1));
    vector <int> ret;
    if (k > k_max_) return ret;
    vector < unordered_set<pii, pair_hash> > idx_tr(G.n);
    for (int curtm = G.tmax - 1; curtm >= ts; curtm--){
        // cout << curtm << ' ' << ts << ' ' << u << ' ' << k << endl;
        // cout << "Current time = " << curtm << endl;
        for (auto e : idx_add[k][curtm]){
            int w = e.first, u = e.second.first, v = e.second.second;
            // cout << "Insert: " << w << ' '  << u << ' ' << v << endl;
            idx_tr[u].insert(make_pair(w, v));
            idx_tr[v].insert(make_pair(w, u));
            // if (u == X || v == X) cout << "inserted, " << idx_tr[X].size() << ' ' << w << ' ' << u << ' ' << v << endl;
        }
        for (auto e : idx_del[k][curtm]){
            int w = e.first, u = e.second.first, v = e.second.second;
            // cout << "Erase: " << w << ' ' << u << ' ' << v << endl;
            idx_tr[u].erase(make_pair(w, v));
            idx_tr[v].erase(make_pair(w, u));
            // if (u == X || v == X) cout << "erased, " << idx_tr[X].size() << ' ' << w << ' ' << u << ' ' << v << endl;
        }
        // cout << "after modified, ts = " << curtm << ", size = " << idx_tr[u].size() << endl;
        if (curtm <= te){
            piii tmp = query_time_anchored_idxgiven_opt(curtm, te, u, k, idx_tr);
            if (tmp.first >= ret_info.first && tmp.first){
                ret_info = tmp;
                ret = cs_search(ret_info.second.first, ret_info.second.second, u, k, idx_tr);
                // if (ret.size() == 1){
                //     cout << "Out edges[" << u << "], (size = " << idx_tr[u].size() << ")" << endl;
                //     cout << u << ' ' << k << ' ' << curtm << endl;
                //     for (auto item : idx_tr[u])
                //         cout << "    ! " << item.first << ' ' << item.second << endl;
                //     cout << "Coretime: " << query_coretime(u, k, ret_info.second.first) << endl;
                //     cout << "Duration = " << ret_info.first << ' ' << ret_info.second.first << ' ' << ret_info.second.second << endl;
                // }
            }
        }
    }
    for (int i = 0; i < G.n; i++)
        idx_tr[i].clear();
    return ret;
}


void mgraph::handle_queries_compressed_cs(string output_filepath){
    cout << "Solving by advanced query()-cs on advanced index" << endl;
    unsigned long long start_time = currentTime();
    ofstream outputFile(output_filepath);
    assert(outputFile.is_open());
    int Q = G.queries.size(), Q_5 = Q / 5, cnt = 0, percent = 0;
    for (auto q : G.queries){
        int ts = q.first.first, te = q.first.second;
        int u = q.second.first, k = q.second.second;
        piii ans_info; vector <int> ans;
        ans = this->query_cs_compressed(ts, te, u, k, ans_info);
        // printf("The community: L = %d, R = %d, with lasting time = %d\n", 
        //         ans.second.first, ans.second.second, ans.first);
        outputFile << "Query: " << ts << ' ' << te << ' ' << u << ' ' << k << endl; 
        outputFile << "  lasting time = " << ans_info.first << "  ";
        outputFile << "Com_ts = " << ans_info.second.first << ", Com_te = " << ans_info.second.second << endl;
        outputFile << "Community size = " << ans.size() << endl;
        for (int i = 0, si = ans.size(); i < min(si, 100) && ans_info.first; i++){
            outputFile << ans[i] << (i == si - 1 ? '\n' : ' ');
            // cout << "coretime = " << query_coretime(ans[i], k, ans_info.second.first) << ' ' << ans_info.second.second << endl;
            assert(query_coretime(ans[i], k, ans_info.second.first) <= ans_info.second.second);
        }

        cnt++;
        if (Q_5 && cnt % Q_5 == 0){
            percent += 20;
            cout << "Finished " << percent << "% queries." << endl;
        }
    }
    outputFile.close();
    unsigned long long end_time = currentTime();
    ofstream logFile(log_path.c_str(), ios::app);
    logFile << "Time used in queries: " << timeFormatting((end_time - start_time) / 1e6).str() << endl;
    cout << "Time used in queries: " << timeFormatting((end_time - start_time) / 1e6).str() << endl;
}

long long mgraph::hash_set(const set <int> &S){
    // long long mod = 998244353, C1 = 133337, C2 = 1337, C3 = 178329197;
    // long long hashval = 0;
    // long long cnt = 0;
    // for (auto i : S){
    //     cnt++;
    //     (hashval += (cnt * cnt % mod * C1 % mod + cnt * C2 % mod + C3) * i % mod) %= mod;
    // }
    // return hashval;
}

int mgraph::query_distinct_core_time_anchored(int t_st, int t_ed, int u, int k,
                const vector <unordered_set<pii, pair_hash> > &idx_tr,
                const vector <long long> &randval,
                unordered_set <long long> &kcore){
    int core_t_u = query_coretime(u, k, t_st);
    if (core_t_u >= G.tmax || core_t_u > t_ed){
        // cout << "Out of query" << endl;
        return 0;
    }
    
    long long hsh_cur = 0;

    queue <int> que;
    vector <int> ctm(G.n);
    vector <bool> vis(G.n);
    vector <int> vislist;
    for (int i = 0; i < G.n; i++) 
        ctm[i] = vis[i] = 0;
    que.push(u), ctm[u] = core_t_u, vis[u] = true, vislist.push_back(u);
    while (que.size()){
        int cu = que.front(); que.pop();
        for (auto e_out : idx_tr[cu]){
            int t = e_out.first, v = e_out.second;
            assert(0 <= t && t < G.tmax);
            assert(0 <= v && v < G.n);
            // assert(t >= query_coretime(v, k, t_st)); // TODO
            if (vis[v]) continue; // visited
            que.push(v), ctm[v] = max(ctm[cu], t), vis[v] = true, vislist.push_back(v);
        }
    }
    vector <pii> vec;
    for (int i = 0; i < G.n; i++)
        if (vis[i] && ctm[i] <= t_ed) vec.push_back(make_pair(ctm[i], i));
    vec.push_back(make_pair(t_ed + 1, -1));
    sort(vec.begin(), vec.end());
    int ret = 0;
    long long mod = 3e18 + 7;
    for (int i = 0, si = vec.size(); i < si - 1; i++){
        (hsh_cur += (__int128_t) randval[vec[i].second] * vec[i].second % mod) %= mod;
        int duration = vec[i + 1].first - vec[i].first;
        if (duration == 0) continue;
        if (kcore.find(hsh_cur) == kcore.end()){
            kcore.insert(hsh_cur);
            ret += duration;
        }
    }
    return ret;
}

double mgraph::query_distinct_core(int ts, int te, int u, int k){
    if (k > k_max_) return 0;
    unordered_set <long long> kcore;
    long long total_duration = 0;
    vector < unordered_set<pii, pair_hash> > idx_tr(G.n);
    vector <long long> randval(G.n);
    for (int i = 0; i < G.n; i++)
        randval[i] = max(1ll * rand() * rand() - rand(), 0ll);
    for (int curtm = G.tmax - 1; curtm >= ts; curtm--){
        // cout << curtm << ' ' << ts << ' ' << u << ' ' << k << endl;
        // cout << "Current time = " << curtm << endl;
        for (auto e : idx_add[k][curtm]){
            int w = e.first, u = e.second.first, v = e.second.second;
            // cout << "Insert: " << w << ' '  << u << ' ' << v << endl;
            idx_tr[u].insert(make_pair(w, v));
            idx_tr[v].insert(make_pair(w, u));
        }
        for (auto e : idx_del[k][curtm]){
            int w = e.first, u = e.second.first, v = e.second.second;
            // cout << "Erase: " << w << ' ' << u << ' ' << v << endl;
            idx_tr[u].erase(make_pair(w, v));
            idx_tr[v].erase(make_pair(w, u));
        }
        if (curtm <= te){
            total_duration += query_distinct_core_time_anchored(
                                    curtm, te, u, k, idx_tr, randval, kcore);
        }
    }
    return (double) total_duration / kcore.size();
}

void mgraph::handle_queries_distinct_core(string output_filepath){
    cout << "Solving by disctinct core on advanced index" << endl;
    unsigned long long start_time = currentTime();
    ofstream outputFile(output_filepath);
    assert(outputFile.is_open());
    int Q = G.queries.size(), Q_5 = Q / 5, cnt = 0, percent = 0;
    for (auto q : G.queries){
        int ts = q.first.first, te = q.first.second;
        int u = q.second.first, k = q.second.second;
        double ans = query_distinct_core(ts, te, u, k);
        outputFile << "On query interval: [" << ts << ", " << te << "]" << endl;
        outputFile << "Average duration = " << ans << endl;
        cnt++;
        if (Q_5 && cnt % Q_5 == 0){
            percent += 20;
            cout << "Finished " << percent << "% queries." << endl;
        }
    }
    outputFile.close();
    unsigned long long end_time = currentTime();
    ofstream logFile(log_path.c_str(), ios::app);
    logFile << "Time used in queries: " << timeFormatting((end_time - start_time) / 1e6).str() << endl;
    cout << "Time used in queries: " << timeFormatting((end_time - start_time) / 1e6).str() << endl;
}

// For debuging and expriment

void mgraph::specific_index(int k, int ts){
    cout << "Start building the ASF-index" << endl;
    treeedge = new vector<vector<vector<pii>>>[k_max_ + 1];
    treeedge[k].resize(G.tmax);
    build_index_tree_time_anchored(k, ts);
    cout << "Finish building the ASF-index" << endl;
}

// It's only used for debuging. 
void mgraph::debug(){
    // int u = 169433;
    // bool flag = false;
    // for (int i = 0; i < G.n; i++)
    //     flag |= treeedge[12][65][i].size();
    // if (flag == false) cout << "non-builded" << endl;
    // cout << "coretime: " << query_coretime(169433, 12, 65) << endl;
    // cout << "!!" << treeedge[12][65][169433].size() << endl;
    // cout << "!!! " << treeedge[12][65][169433][0].first << ' ' << treeedge[12][65][169433][0].second << endl;
    // vector <int> ret, vis(G.n);
    // for (int i = 0; i < G.n; i++)
    //     vis[i] = 0;
    // queue <int> que;
    // que.push(u), vis[u] = true, ret.push_back(u);
    // while (que.size()){
    //     int cu = que.front(); que.pop();
    //     for (auto e_out : treeedge[12][65][cu]){
    //         int t = e_out.first, v = e_out.second;
    //         if (vis[v]) continue;
    //         if (t <= 68) que.push(v), vis[v] = true, ret.push_back(v);
    //     }
    // }
    // cout << ret.size() << endl;
    // // for (int i = 0; i < ret.size(); i++)
    //     // cout << ret[i] << (i == ret.size() - 1 ? '\n' : ' ');
    
    // cout << "-------" << endl;
    // for (int i = 0; i < core_t_[u][12].size(); i++)
    //     cout << core_t_[u][12][i].first << ' ' << core_t_[u][12][i].second << endl;
}

void mgraph::output_topnodes(int topk){
    vector <pii> vec;
    for (int i = 0; i < G.n; i++)
        vec.push_back(make_pair(core_t_[i].size(), i));
    sort(vec.begin(), vec.end()), reverse(vec.begin(), vec.end());
    for (int i = 0; i < topk; i++)
        cout << vec[i].second << ", ";
    cout << endl;
}