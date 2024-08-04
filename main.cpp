#include "mgraph.h"

int main(int argc, char *argv[]){
    assert(argc > 3);

    string graph_path = (string) argv[1];
    string query_path = (string) argv[2];
    string answer_path = (string) argv[3];

    bool flag_ASF = false;
    bool flag_AG = false;
    bool flag_BITE = false;
    bool flag_AITE = false;

    bool flag_load_core_index = false;
    string load_core_index_path;
    bool flag_save_core_index = false;
    string save_core_index_path;
    bool flag_load_ASForAG_index = false;
    string load_ASF_AG_index_path;
    bool flag_save_ASForAG_index = false;
    string save_ASF_AG_index_path;

    for (int i = 4; i < argc; i++){
        if (argv[i][0] != '-') continue;
        string op = (string) argv[i];
        if (op == "-ASF") 
            flag_ASF = true;
        if (op == "-AG")
            flag_AG = true;
        if (op == "-BITE")
            flag_BITE = true;
        if (op == "-AITE")
            flag_AITE = true;
        if (op == "-sc"){
            assert(i + 1 < argc);
            flag_save_core_index = true;
            save_core_index_path = (string) argv[++i];
        }
        if (op == "-lc"){
            assert(i + 1 < argc);
            flag_load_core_index = true;
            load_core_index_path = (string) argv[++i];
        }
        if (op == "-st"){
            assert(i + 1 < argc);
            flag_save_ASForAG_index = true;
            save_ASF_AG_index_path = (string) argv[++i];
        }
        if (op == "-lt"){
            assert(i + 1 < argc);
            flag_load_ASForAG_index = true;
            load_ASF_AG_index_path = (string) argv[++i];
        }
    }
    assert(flag_AG ^ flag_ASF);
    assert(flag_BITE ^ flag_AITE);
    assert(!(flag_load_core_index && flag_save_core_index));
    assert(!(flag_load_ASForAG_index && flag_save_ASForAG_index));

    mgraph mg;
    unsigned long long start_time_index_construction = currentTime();
    // build or load(save) core-index
    if (flag_load_core_index){
        mg.load_idx(load_core_index_path);
    } else if (flag_save_core_index){
        mg.build_coreness_index(graph_path, save_core_index_path);
        mg.load_idx(save_core_index_path);
    } else {
        mg.build_coreness_index(graph_path, core_idx_path);
        mg.load_idx(core_idx_path);
    }
    mg.load_graph(graph_path, query_path);
    // build or load(save) ASF/AG-index
    if (flag_ASF){
        if (flag_load_ASForAG_index){
            mg.load_index_trees(load_ASF_AG_index_path);
        } else {
            mg.build_index_trees();
        }
        if (flag_save_ASForAG_index){
            mg.save_index_trees(save_ASF_AG_index_path);
        }
    }
    if (flag_AG){
        if (flag_load_ASForAG_index){
            mg.load_compressed_index(load_ASF_AG_index_path);
        } else {
            mg.build_compressed_index();
        }
        if (flag_save_ASForAG_index){
            mg.save_compressed_index(save_ASF_AG_index_path);
        }
    }
    unsigned long long end_time_index_construction = currentTime();

    if (flag_AG && flag_AITE){
        mg.handle_queries_compressed(answer_path);
    } else if (flag_AG && flag_BITE){
        mg.handle_queries_compressed_basic(answer_path);
    } else if (flag_ASF && flag_AITE){
        mg.handle_queries(answer_path);
    } else if (flag_ASF && flag_BITE){
        mg.handle_queries_basic(answer_path);
    } else assert(false);

    unsigned long long period = end_time_index_construction - start_time_index_construction;
    cout << "Time on index construction(loading): " << timeFormatting(period / 1e6).str() << endl;
    return 0;
}

