#ifndef LCT
#define LCT

#include "commonfunctions.h"

class Link_Cut_Tree{

private:

    int *son[2];
    int *fa;
    int *tag;
    int *val;
    int *maxval;
    int *maxid;

    bool Rc(int x);
    bool nroot(int x);
    void push_up(int x);
    void push_down(int x);
    void rotate(int x);
    void splay(int x);
    void access(int x);
    void find(int x);

    void makeroot(int x);
    int findroot(int x);
    int split(int x, int y);
    void link(int x, int y);
    void cut(int x, int y);

	bool notconnected(int u, int v);
	void init_node(int id, int v);

public:
    void init(int maxn);
    void release();

	bool try_insert(piii e, int &de_id);
    void insert(piii e, int eid); // suppose that (u, v) is not connected.
    void modify(int id, int v);
};

class Minimum_Spanning_Tree{

private:
    Link_Cut_Tree T;
	queue <int> freeids;
	unordered_map <pii, int, pair_hash> e2id;
	unordered_map <int, pii> id2e;
	unordered_map <int, int> id2val;
public:
    
	Minimum_Spanning_Tree(int n);
    ~Minimum_Spanning_Tree();
    bool insert(piii ad, piii &de);

};

#endif