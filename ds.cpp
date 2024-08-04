#include "ds.h"

// For MST maintainness

Minimum_Spanning_Tree::Minimum_Spanning_Tree(int n){
    int siz = n * 2 + 100;
    T.init(siz);
    for (int i = 1; i <= n; i++)
        freeids.push(i + n);
}

bool Minimum_Spanning_Tree::insert(piii ad, piii &de){
    piii tr_ad = make_pair(ad.first, 
        make_pair(ad.second.first + 1, ad.second.second + 1));

    // cout << "Try insert (transed): " << tr_ad.first << ", (" << tr_ad.second.first << ", " << tr_ad.second.second << ")" << endl;

    // check if the edge exists.
    int fid = e2id[tr_ad.second];
    if (fid != 0){
        if (tr_ad.first >= id2val[fid])
            return false;
        de = make_pair(id2val[fid], ad.second);
        id2val[fid] = tr_ad.first;
        T.modify(fid, tr_ad.first);
        return true;
    }

    // try insert the edge.
    int did = 0;
    if (T.try_insert(tr_ad, did) == false)
        return false;
    // assign id for the new inserted edge.
    assert(freeids.size());
    int nid = freeids.front();
    freeids.pop();
    // cout << "Used, deleted: " << nid << ' ' << did << endl;
    e2id[tr_ad.second] = nid;
    id2e[nid] = tr_ad.second;
    id2val[nid] = tr_ad.first;
    T.insert(tr_ad, nid);
    // get the deleted edge
    if (did != -1){
        assert(e2id[id2e[did]]);
        e2id[id2e[did]] = 0; // mark the deleted edge as non-exist.
        de = make_pair(id2val[did], id2e[did]);
        de.second.first -= 1, de.second.second -= 1;
        freeids.push(did);
    } else de = make_pair(-1, make_pair(-1, -1));
    return true;
}

Minimum_Spanning_Tree::~Minimum_Spanning_Tree(){
    T.release();
}


// LCT part


bool Link_Cut_Tree::Rc(int x){
    return son[1][fa[x]] == x;
}
bool Link_Cut_Tree::nroot(int x){
    return son[0][fa[x]] == x || Rc(x);
}
void Link_Cut_Tree::push_up(int x){
    maxval[x] = val[x], maxid[x] = x;
    int son1 = son[0][x], son2 = son[1][x];
    if (maxval[son1] > maxval[x])
        maxval[x] = maxval[son1], maxid[x] = maxid[son1];
    if (maxval[son2] > maxval[x])
        maxval[x] = maxval[son2], maxid[x] = maxid[son2];
}

void Link_Cut_Tree::push_down(int x){
    if (!tag[x]) return;
    int son1 = son[0][x], son2 = son[1][x];
    tag[son1] ^= 1, tag[son2] ^= 1;
    swap(son[0][x], son[1][x]), tag[x] = 0;
}

void Link_Cut_Tree::rotate(int x){
    int y = fa[x], z = fa[y];
    push_down(y), push_down(x);
    int tg = Rc(x), w = son[tg ^ 1][x];
    son[tg][y] = w, fa[w] = y;
    if (nroot(y)) son[Rc(y)][z] = x;
    fa[x] = z;
    son[tg ^ 1][x] = y, fa[y] = x;
    push_up(y), push_up(x);
}

void Link_Cut_Tree::splay(int x){
    // cout << "In splay: " << x << endl;
    push_down(x);
    for (int y = fa[x]; nroot(x); rotate(x), y = fa[x]){
        // cout << "Splaying " << x << ' ' << fa[x] << ' ' << "...\n";
        if (nroot(y)) rotate(Rc(y) == Rc(x) ? y : x);
        // cout << "After Splaying " << x << ' ' << fa[x] << ' ' << "...\n";
    }
}

void Link_Cut_Tree::access(int x){
    // cout << "In access: " << x << endl;
	for (int y = 0; x; y = x, x = fa[x]){
		splay(x), son[1][x] = y, push_up(x);
        // cout << "Accessing " << x << "...\n";
    }
    // cout << "Finish access." << endl;
}

void Link_Cut_Tree::find(int x){
    access(x); 
    splay(x);
}

void Link_Cut_Tree::makeroot(int x){
    find(x);
    tag[x] ^= 1; 
}

int Link_Cut_Tree::findroot(int x){
    find(x);
    push_down(x);
    while (son[0][x])
        x = son[0][x], push_down(x);
    splay(x); 
    return x;
}

int Link_Cut_Tree::split(int x, int y){
    makeroot(x);
    find(y); 
    return y;
}

void Link_Cut_Tree::link(int x, int y){
    makeroot(x); 
    if (findroot(y) != x) fa[x] = y;
}

void Link_Cut_Tree::cut(int x, int y){
    makeroot(x);
    if (findroot(y) != x || fa[y] != x || son[0][y]) return;
    fa[y] = 0, son[1][x] = 0, push_up(x);
}

bool Link_Cut_Tree::notconnected(int u, int v){
    makeroot(u);
    bool ret = findroot(v) != u;
    return ret;
}

void Link_Cut_Tree::init_node(int id, int v){
    son[0][id] = son[1][id] = fa[id] = tag[id] = 0;
    maxval[id] = val[id] = v;
    maxid[id] = id;
}

void Link_Cut_Tree::init(int maxn){
    maxn += 5;
    // cout << "LCT init, size = " << maxn << endl;
    son[0] = new int[maxn];
    son[1] = new int[maxn];
    fa = new int[maxn];
    tag = new int[maxn];
    val = new int[maxn];
    maxval = new int[maxn];
    maxid = new int[maxn];
    memset(son[0], 0, maxn * sizeof(int));
    memset(son[1], 0, maxn * sizeof(int));
    memset(fa, 0, maxn * sizeof(int));
    memset(tag, 0, maxn * sizeof(int));
    memset(val, 0, maxn * sizeof(int));
    memset(maxval, 0, maxn * sizeof(int));
    memset(maxid, 0, maxn * sizeof(int));
}

bool Link_Cut_Tree::try_insert(piii e, int &did){
    int w = e.first, u = e.second.first, v = e.second.second;
    if (notconnected(u, v)){
        did = -1;
        return true;
    }
    split(u, v);
    int mxid = maxid[v];
    int mxval = maxval[mxid];
    if (mxval <= w) 
        return false;
    splay(mxid);
    fa[son[0][mxid]] = fa[son[1][mxid]] = 0;
    did = mxid;
    return true;
}

void Link_Cut_Tree::insert(piii e, int nid){
    int w = e.first, u = e.second.first, v = e.second.second;
    init_node(nid, w);
    link(u, nid);
    link(nid, v);
}

void Link_Cut_Tree::modify(int u, int w){
    find(u);
    val[u] = w;
    push_up(u);
}

void Link_Cut_Tree::release(){
    delete [] son[0];
    delete [] son[1];
    delete [] fa;
    delete [] tag;
    delete [] val;
    delete [] maxval;
    delete [] maxid;
}
