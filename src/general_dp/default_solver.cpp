// tree dp for every budget j <= K at once, no lagrange, so it's exact for every
// K (cost(k) isn't convex on trees). same chains as lagrange_solver, but cur is
// an array over j. fine for a few thousand nodes and K up to a few hundred.
// cur[j](x) = best cost of the suffix with at most j cuts, carrying x
// ./default_solver [K] [ALLK=0] [ROOT=2] < tree
// ALLK=1 also prints "allk: v0 v1 ... vK"

#include <bits/stdc++.h>
#include "slope_min_helper.h"
#define int long long
#define rep(i,a,b) for(int i=(a);i<(b);++i)
#define pb push_back
#define ff first
#define ss second
using namespace std;
typedef pair<int,int> pii;

const int N = 400005;
int val[N];
int p[N];
int T = 0;
int tin[N], tout[N];
int sz[N];    // subtree weight
int cnt_[N];  // subtree node count
vector<int> g[N];
int here[N];
int ROOT = 1;
int n, K;

bool g_preferMoreCuts = false; // cut counts aren't used here

// heap peak in KiB, -1 without mem_cap.h
static long long peak_heap_kb(){
#ifdef MEM_CAP_H
    return (long long)(mem_cap_peak_bytes() >> 10);
#else
    return -1;
#endif
}

/****** fast input ******/
static char ibuf[1 << 16];
static size_t ilen = 0, ipos = 0;
static inline int igc() {
    if (ipos == ilen) { ilen = fread(ibuf, 1, sizeof(ibuf), stdin); ipos = 0; if (!ilen) return -1; }
    return ibuf[ipos++];
}
static long long readInt() {
    int c = igc(); while (c != '-' && (c < '0' || c > '9')) { if (c == -1) return 0; c = igc(); }
    int sg = 1; if (c == '-') { sg = -1; c = igc(); }
    long long x = 0;
    while (c >= '0' && c <= '9') { x = x * 10 + (c - '0'); c = igc(); }
    return x * sg;
}

void prep(int root){
    static vector<pair<int,int>> st;
    st.clear();
    st.pb({root, 0});
    T = 0;
    tin[root] = ++T; here[T] = root; sz[root] = val[root]; cnt_[root] = 1;
    while(!st.empty()){
        auto &top = st.back();
        int v = top.first;
        if(top.second < (int)g[v].size()){
            int to = g[v][top.second++];
            if(to == p[v]) continue;
            tin[to] = ++T; here[T] = to; sz[to] = val[to]; cnt_[to] = 1;
            st.pb({to, 0});
        } else {
            tout[v] = T;
            st.pop_back();
            if(!st.empty()){ sz[st.back().first] += sz[v]; cnt_[st.back().first] += cnt_[v]; }
        }
    }
}

int prefV[N];
long long Wtot;
static inline long long xmax_at(int t, int h){
    return Wtot - (prefV[tout[h]] - prefV[t - 1]);
}

vector<int> chainHead;
vector<int> leafPref;

void heavy_last_reorder(){
    for(int v = 1; v <= T; v++){
        sort(g[v].begin(), g[v].end(), [](int a, int b){ return cnt_[a] < cnt_[b]; });
    }
    prep(ROOT);
}

void init_dp(){
    int nn = T;
    chainHead.assign(nn + 1, 0);
    leafPref.assign(nn + 1, 0);
    for(int t = 1; t <= nn; t++){
        int v = here[t];
        chainHead[v] = (v == ROOT || tout[v] != tout[p[v]]) ? v : chainHead[p[v]];
        leafPref[t] = leafPref[t - 1] + (tout[v] == t ? 1 : 0);
    }
}

// dpv[v][j] = best cost of v's subtree with at most j cuts, j <= min(K, cnt_[v]-1)
// effK[v] = last j where dpv[v] still decreases
vector<vector<long long>> dpv;
vector<int> effK;
vector<char> done_;

const vector<long long>& solve(int v); // fwd

void build_chain(int h){
    int sub = tout[h] - tin[h] + 1;
    int leaf0 = leafPref[tin[h] - 1];
    int nsnap = leafPref[tout[h] - 1] - leaf0;

    int capMax = (int)min((long long)K, (long long)sub);
    vector<DPFunc> cur(capMax + 1); // cur[j] for j <= cap
    int cap = 0;                    // past the end every j is the same

    struct Snap { vector<DPFunc> f; int cap = 0; };
    vector<Snap> snaps(nsnap);
    vector<pii> relStack;           // (release off, snap idx)

    for(int off = sub - 1; off >= 0; off--){
        int t  = tin[h] + off;
        int at = here[t];
        long long w = val[at];
        long long xm = xmax_at(t, h);

        // j >= #nodes of the suffix are all equal, so cap grows by 1 at most
        int capNew = (int)min((long long)K, (long long)(sub - off));
        if(capNew > cap){ cur[capNew] = cur[cap]; cap = capNew; }

        if(tout[at] == tout[h]){
            // chain node: the suffix is at's subtree, so dpv[at] is read here
            for(int j = 0; j <= cap; j++) cur[j].noCut(w);
            int mk = (int)min((long long)K, cnt_[at] - 1);
            auto &d = dpv[at];
            d.resize(mk + 1);
            for(int j = 0; j <= mk; j++) d[j] = cur[j].at0().first;
            int e = 0; for(int j = 1; j <= mk; j++) if(d[j] < d[j-1]) e = j;
            effK[at] = e;
            done_[at] = 1;
            if(at != h){
                // cut above at: resumes past the end, a constant
                for(int j = cap; j >= 1; j--)
                    cur[j].mergeMinConst(d[min((long long)(j - 1), (long long)mk)], 0, xm);
            }
        } else if(tout[at] == t){
            // light leaf: cutting reads cur[j-1] before its update, so go down in j
            const auto &d = solve(at);
            for(int j = cap; j >= 1; j--){
                cur[j].noCut(w);
                cur[j].mergeMinWithOffsets(cur[j - 1], d[0], 0, xm);
            }
            cur[0].noCut(w);
        } else {
            // light subtree: j1 cuts inside at, j-1-j1 after it
            const auto &d = solve(at);
            int q = tout[at] + 1;
            Snap &S = snaps[leafPref[q - 1] - leaf0 - 1];
            for(int j = 0; j <= cap; j++) cur[j].noCut(w);
            int mk = effK[at]; // larger j1 can't help
            for(int j = cap; j >= 1; j--){
                int j1lo = (int)max(0LL, (long long)(j - 1) - S.cap);
                int j1hi = (int)min((long long)mk, (long long)(j - 1));
                if(j1lo > j1hi && j - 1 > S.cap){
                    // both sides already saturated
                    j1lo = j1hi = mk;
                }
                for(int j1 = j1lo; j1 <= j1hi; j1++){
                    int j2 = (int)min((long long)(j - 1 - j1), (long long)S.cap);
                    cur[j].mergeMinWithOffsets(S.f[j2], d[j1], 0, xm);
                }
            }
        }

        // copy the whole array after a leaf (single-node chains don't need it)
        if(off >= 1 && tout[here[t - 1]] == t - 1){
            int lf = here[t - 1];
            if(chainHead[lf] != lf){
                int sidx = leafPref[t - 1] - leaf0 - 1;
                snaps[sidx].f.assign(cur.begin(), cur.begin() + cap + 1);
                snaps[sidx].cap = cap;
                relStack.pb({tin[chainHead[lf]] - tin[h], sidx});
            }
        }
        while(!relStack.empty() && relStack.back().ff == off){
            snaps[relStack.back().ss].f.clear();
            snaps[relStack.back().ss].f.shrink_to_fit();
            relStack.pop_back();
        }
    }
}

const vector<long long>& solve(int v){
    if(!done_[v]) build_chain(chainHead[v]);
    return dpv[v];
}

/****** RECOVERY ******/
// same as build_chain but keeps f[off][j] for every position
void build_full(int h, vector<vector<DPFunc>>& f, vector<int>& capAt){
    int sub = tout[h] - tin[h] + 1;
    f.assign(sub + 1, {});
    capAt.assign(sub + 1, 0);
    f[sub].assign(1, DPFunc());     // past the end

    for(int off = sub - 1; off >= 0; off--){
        int t  = tin[h] + off;
        int at = here[t];
        long long w = val[at];
        long long xm = xmax_at(t, h);
        int capOld = capAt[off + 1];
        int cap = (int)min((long long)K, (long long)(sub - off));
        capAt[off] = cap;

        f[off].assign(cap + 1, DPFunc());
        for(int j = 0; j <= cap; j++) f[off][j] = f[off + 1][min((long long)j, (long long)capOld)];
        for(int j = 0; j <= cap; j++) f[off][j].noCut(w);

        if(tout[at] == tout[h]){
            if(at != h){
                const auto &d = dpv[at];
                int mk = (int)d.size() - 1;
                for(int j = cap; j >= 1; j--)
                    f[off][j].mergeMinConst(d[min((long long)(j - 1), (long long)mk)], 0, xm);
            }
        } else if(tout[at] == t){
            const auto &d = dpv[at];
            for(int j = cap; j >= 1; j--)
                f[off][j].mergeMinWithOffsets(f[off + 1][min((long long)(j - 1), (long long)capOld)],
                                              d[0], 0, xm);
        } else {
            const auto &d = dpv[at];
            int qoff = (int)(tout[at] + 1 - tin[h]);
            int Scap = capAt[qoff];
            int mk = effK[at];
            for(int j = cap; j >= 1; j--){
                int j1lo = (int)max(0LL, (long long)(j - 1) - Scap);
                int j1hi = (int)min((long long)mk, (long long)(j - 1));
                if(j1lo > j1hi && j - 1 > Scap) j1lo = j1hi = mk;
                for(int j1 = j1lo; j1 <= j1hi; j1++){
                    int j2 = (int)min((long long)(j - 1 - j1), (long long)Scap);
                    f[off][j].mergeMinWithOffsets(f[qoff][j2], d[j1], 0, xm);
                }
            }
        }
    }
}

vector<int> recoverCutVert;
// worklist of (component root, budget of its subtree)
void recover(int root, int budget){
    vector<pii> q;
    q.pb({root, budget});
    while(!q.empty()){
        int v = q.back().ff, jbud = q.back().ss; q.pop_back();
        vector<vector<DPFunc>> f;
        vector<int> capAt;
        build_full(v, f, capAt);

        long long w = 0;
        int j = jbud;
        int t = tin[v];
        while(t != -1 && t <= tout[v]){
            int off = t - tin[v];
            int at  = here[t];
            long long cost = f[off][min((long long)j, (long long)capAt[off])].eval(w).first;

            // option 1: no cut above at
            int jn = (int)min((long long)j, (long long)capAt[off + 1]);
            long long costNoCut = f[off + 1][jn].eval(w + val[at]).first
                                + val[at] * val[at] + 2 * w * val[at];
            if(at == v || costNoCut == cost){
                w += val[at];
                t++;
                continue;
            }

            // option 2: cut above at, j1 cuts inside its subtree
            int nv   = (tout[at] + 1 <= tout[v]) ? (int)(tout[at] + 1) : -1;
            int qoff = (nv == -1) ? -1 : (int)(nv - tin[v]);
            int Scap = (nv == -1) ? 0 : capAt[qoff];
            const auto &d = dpv[at];
            int mk = (int)d.size() - 1;
            bool found = false;
            for(int j1 = 0; j1 <= min((long long)mk, (long long)(j - 1)); j1++){
                int j2 = (int)min((long long)(j - 1 - j1), (long long)Scap);
                long long after = (nv == -1) ? 0 : f[qoff][j2].eval(w).first;
                if(d[j1] + after == cost){
                    recoverCutVert.pb(at);
                    q.pb({at, j1});
                    j = j - 1 - j1;
                    t = nv;
                    found = true;
                    break;
                }
            }
            if(!found){
                fprintf(stderr, "recover: no matching option at node %lld\n", (long long)at);
                abort();
            }
        }
    }
}

long long cost_recovered(const vector<int>& cut){
    vector<char> isCut(T + 1, 0);
    for(int v : cut) isCut[v] = 1;
    vector<long long> sub(T + 1, 0);
    long long res = 0;
    for(int t = T; t >= 1; t--){
        int v = here[t];
        sub[v] += val[v];
        if(v == ROOT || isCut[v]) res += sub[v] * sub[v];
        else sub[p[v]] += sub[v];
    }
    return res;
}

/****** rooting, same as lagrange_solver ******/
int pick_centroid(const vector<vector<int>>& adj){
    int nn = n;
    vector<int> order; order.reserve(nn);
    vector<int> par(nn + 1, 0), cnt(nn + 1, 1);
    vector<char> seen(nn + 1, 0);
    order.pb(1); seen[1] = 1;
    rep(i,0,nn){
        int v = order[i];
        for(int to : adj[v]) if(!seen[to]){ seen[to] = 1; par[to] = v; order.pb(to); }
    }
    for(int i = nn - 1; i >= 1; i--) cnt[par[order[i]]] += cnt[order[i]];
    int best = 1; long long bestVal = LLONG_MAX;
    rep(i,0,nn){
        int v = order[i];
        long long mx = nn - cnt[v];
        for(int to : adj[v]) if(to != par[v]) mx = max(mx, (long long)cnt[to]);
        if(mx < bestVal){ bestVal = mx; best = v; }
    }
    return best;
}

void root_tree(int root, const vector<vector<int>>& adj){
    ROOT = root;
    rep(i,1,n+1) g[i].clear();
    vector<int> order; order.reserve(n);
    vector<char> seen(n + 1, 0);
    order.pb(root); seen[root] = 1; p[root] = 0;
    rep(i,0,n){
        int v = order[i];
        for(int to : adj[v]) if(!seen[to]){
            seen[to] = 1; p[to] = v; g[v].pb(to); order.pb(to);
        }
    }
}

#undef int
int main(int argc, char* argv[]){
#define int long long
    int k;
    n = readInt(); k = readInt();
    if(argc > 1) k = atoll(argv[1]);
    int allk = (argc > 2) ? atoll(argv[2]) : 0;
    int rootMode = (argc > 3) ? atoll(argv[3]) : 2; // 2 = auto, 1 = centroid, 0 = node 1
    rep(i,1,n+1) val[i] = readInt();

    K = max(0LL, min(k, n - 1));

    {
        vector<vector<int>> adj(n + 1);
        rep(i,2,n+1){ int pi = readInt(); adj[pi].pb(i); adj[i].pb(pi); }
        auto light_work = [&](int r)->long long {
            root_tree(r, adj);
            prep(ROOT);
            heavy_last_reorder();
            init_dp();
            long long tot = 0;
            for(int t = 1; t <= T; t++){
                int v = here[t];
                if(chainHead[v] == v) tot += tout[v] - tin[v] + 1;
            }
            return tot - T;
        };
        int root = 1;
        if(rootMode == 1) root = pick_centroid(adj);
        else if(rootMode >= 2){
            int c = pick_centroid(adj);
            root = 1;
            if(c != 1){
                long long w1 = light_work(1);
                long long wc = light_work(c);
                root = (wc < w1) ? c : 1;
            }
        }
        root_tree(root, adj);
    }

    prep(ROOT);
    heavy_last_reorder();
    init_dp();
    rep(i,1,T+1) prefV[i] = prefV[i-1] + val[here[i]];
    Wtot = sz[ROOT];

    // dpv would be too big
    {
        long long cells = 0;
        rep(v,1,T+1) cells += min((long long)K, cnt_[v] - 1) + 1;
        if(cells > 150000000LL){
            fprintf(stderr, "default_solver: dp would need %lld cells (> 1.5e8); "
                            "reduce K or use the Lagrangian lagrange_solver\n", cells);
            return 1;
        }
    }

    dpv.assign(T + 1, {});
    effK.assign(T + 1, 0);
    done_.assign(T + 1, 0);

    // recovery is part of the timed run: the answer is the cut set, not the cost
    auto t2 = clock();
    const auto &d = solve(ROOT);
    long long ans = d[min((long long)K, (long long)d.size() - 1)];
    double dtValue = (clock()-t2)/(double)CLOCKS_PER_SEC;
    long long searchPeakKb = peak_heap_kb();   // dp only

    auto t3 = clock();
    recoverCutVert.clear();
    recover(ROOT, K);
    long long real_cost = cost_recovered(recoverCutVert);
    double dtRecover = (clock()-t3)/(double)CLOCKS_PER_SEC;
    double dt = dtValue + dtRecover;

    assert((long long)recoverCutVert.size() <= K);
    assert(real_cost == ans);

    cout << ans << "  - k=" << k << ", time: " << dt << "s"
         << ", value_time: " << dtValue << "s, recover_time: " << dtRecover << "s"
         << ", cuts: " << recoverCutVert.size()
         << ", status: exact-allk"
         << ", root: " << ROOT;
    if(searchPeakKb >= 0)
        cout << ", search_peak_kb: " << searchPeakKb << ", peak_kb: " << peak_heap_kb();
    cout << "\n";
    if(allk){
        cout << "allk:";
        rep(j,0,K+1) cout << " " << d[min((long long)j, (long long)d.size() - 1)];
        cout << "\n";
    }
    return 0;
}
