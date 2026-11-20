// same dp as lagrange_solver_without_hld.cpp, with heavy-light chains
// ./lagrange_solver [K] [SKIP_BRUTE=1] [DUMP_CUTS=0] [ROOT=2] < tree
// ROOT: 2 = auto, 1 = centroid, 0 = node 1 (needed for brute)

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

bool g_preferMoreCuts = true; // tie-break on cut count

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

// iterative, paths can be very deep
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

int Lambda;
int prefV[N];   // prefix sum of val in euler order
long long Wtot;

// x at position t is weight outside [t, tout[h]], so f only matters up to here
static inline long long xmax_at(int t, int h){
    return Wtot - (prefV[tout[h]] - prefV[t - 1]);
}

int curVisDP = 0;
vector<int> visDP;
vector<pii> dp;           // dp[v] = (cost + lambda*cuts, cuts)
vector<int> chainHead;
vector<int> leafPref;     // number of leaves (tout[here[t]] == t) up to t

// biggest child last: then it has the same tout as its parent (same chain)
void heavy_last_reorder(){
    for(int v = 1; v <= T; v++){
        sort(g[v].begin(), g[v].end(), [](int a, int b){ return cnt_[a] < cnt_[b]; });
    }
    prep(ROOT);
}

void init_dp(){
    int nn = T;
    visDP.assign(nn + 1, 0);
    dp.assign(nn + 1, pii(0, 0));
    chainHead.assign(nn + 1, 0);
    leafPref.assign(nn + 1, 0);
    // tin order: parents before children
    for(int t = 1; t <= nn; t++){
        int v = here[t];
        chainHead[v] = (v == ROOT || tout[v] != tout[p[v]]) ? v : chainHead[p[v]];
        leafPref[t] = leafPref[t - 1] + (tout[v] == t ? 1 : 0);
    }
}

pii solve(int v); // fwd

// goes right to left over h's euler range keeping only cur = f[t].
// a cut resumes right after a leaf, so cur is copied only there, and freed
// after its last reader (the head of that leaf's chain). intervals are
// laminar, so the frees are a stack.
void build_chain(int h){
    int sub = tout[h] - tin[h] + 1;
    int leaf0 = leafPref[tin[h] - 1];
    int nsnap = leafPref[tout[h] - 1] - leaf0;

    vector<DPFunc> snaps(nsnap);
    vector<pii> relStack;           // (release off, snap idx)
    relStack.reserve(min((long long)nsnap, 64LL));

    DPFunc cur;                     // past the end
    for(int off = sub - 1; off >= 0; off--){
        int t  = tin[h] + off;
        int at = here[t];
        long long w = val[at];

        if(tout[at] == tout[h]){
            // on h's chain: cutting resumes past the end, so it's a constant
            cur.noCut(w);
            Val b = cur.at0();
            visDP[at] = curVisDP;
            dp[at] = pii(b.first, b.second);
            if(at != h)
                cur.mergeMinConst(Lambda + dp[at].ff, (int)(1 + dp[at].ss), xmax_at(t, h));
        } else if(tout[at] == t){
            // light leaf: cutting resumes at t+1, which is cur itself
            pii cortar = solve(at);
            // w == 0: cutting can't win and noCut does nothing
            if(!(w == 0 && Lambda + cortar.ff > 0))
                cur.cutLeafMerge(w, Lambda + cortar.ff, (int)(1 + cortar.ss), xmax_at(t, h));
        } else {
            // light subtree: cutting resumes at tout[at]+1
            pii cortar = solve(at);
            cur.noCut(w);
            int q = tout[at] + 1;
            int sidx = leafPref[q - 1] - leaf0 - 1;
            cur.mergeMinWithOffsets(snaps[sidx], Lambda + cortar.ff, (int)(1 + cortar.ss),
                                    xmax_at(t, h));
        }

        // a leaf alone in its chain is handled by cutLeafMerge, no copy
        if(off >= 1 && tout[here[t - 1]] == t - 1){
            int lf = here[t - 1];
            if(chainHead[lf] != lf){
                int sidx = leafPref[t - 1] - leaf0 - 1;
                snaps[sidx] = cur;
                relStack.pb({tin[chainHead[lf]] - tin[h], sidx});
            }
        }
        while(!relStack.empty() && relStack.back().ff == off){
            snaps[relStack.back().ss].release();
            relStack.pop_back();
        }
    }
}

pii solve(int v){
    if(visDP[v] == curVisDP) return dp[v];
    build_chain(chainHead[v]);   // fills dp of the whole chain
    return dp[v];
}

void clear(){ curVisDP++; }

pii solve_lamda(int l){
    clear();
    Lambda = l;
    return solve(ROOT);
}

/****** RECOVERY ******/
// same as build_chain but keeps f for every position
void build_full(int h, vector<DPFunc>& f){
    int sub = tout[h] - tin[h] + 1;
    f.assign(sub + 1, DPFunc());
    for(int off = sub - 1; off >= 0; off--){
        int t  = tin[h] + off;
        int at = here[t];
        long long w = val[at];
        if(tout[at] == tout[h]){
            f[off] = f[off + 1];
            f[off].noCut(w);
            Val b = f[off].at0();
            visDP[at] = curVisDP;
            dp[at] = pii(b.first, b.second);
            if(at != h)
                f[off].mergeMinConst(Lambda + dp[at].ff, (int)(1 + dp[at].ss), xmax_at(t, h));
        } else {
            pii cortar = solve(at);
            f[off] = f[off + 1];
            f[off].noCut(w);
            int nv_off = tout[at] + 1 - tin[h];
            f[off].mergeMinWithOffsets(f[nv_off], Lambda + cortar.ff, (int)(1 + cortar.ss),
                                       xmax_at(t, h));
        }
    }
}

vector<int> recoverCutVert;
// one build_full per chain: a cut on the same chain ends at the same tout, so
// f is already its dp and we keep walking it. only cuts at light children
// rebuild, O(n log n) in total (rebuilding per cut was O(nK) on a path).
void recover(int root){
    vector<int> q;
    q.pb(root);
    while(!q.empty()){
        int v = q.back(); q.pop_back();
        vector<DPFunc> f;
        build_full(v, f);
        int croot = v;              // component root, always on v's chain
        long long w = 0;
        int t = tin[v];
        while(t != -1 && t <= tout[croot]){
            int off = t - tin[v];
            int at  = here[t];

            Val best = f[off].eval(w);
            long long cost = best.first;
            int qtd = best.second;

            // option 1: no cut at `at`
            Val bestTo1 = f[off + 1].eval(w + val[at]);
            long long costNoCut = bestTo1.first + val[at] * val[at] + 2 * w * val[at];
            if(at == croot || (costNoCut == cost && bestTo1.second == qtd)){
                w += val[at];
                t++;
                continue;
            }
            // option 2: cut above `at`
            int nv = (tout[at] + 1 <= tout[croot]) ? (int)(tout[at] + 1) : -1;
            pii cortar = solve(at);
            Val after = (nv == -1 ? Val{0, 0} : f[nv - tin[v]].eval(w));
            long long costHere = cortar.ff + after.first + Lambda;
            long long cuts_here = 1 + cortar.ss + (nv == -1 ? 0 : after.second);
            if(costHere == cost && cuts_here == qtd){
                recoverCutVert.pb(at);
                if(tout[at] == tout[croot]){
                    // same chain: keep f, new component starting at at
                    croot = at;
                    w = 0;
                    t = tin[at];
                } else {
                    q.pb(at);      // light child: its own build_full later
                    t = nv;
                }
            } else {
                fprintf(stderr, "recover: no matching option at node %lld\n", (long long)at);
                break;
            }
        }
    }
}

long long cost_recovered(vector<int> cut){
    int n = T;
    vector<int> isCut(n + 1, 0);
    for(int v : cut) isCut[v] = 1;
    long long res = 0;
    vector<long long> sub(n + 1, 0);
    for(int t = T; t >= 1; t--){
        int v = here[t];
        sub[v] += val[v];
        if(v == ROOT || isCut[v]){ res += sub[v] * sub[v]; }
        else sub[p[v]] += sub[v];
    }
    return res;
}

/****** ALIENS ******/
int n;
long long sumW2;      // lambda above this gives 0 cuts

pii eval_lambda(long long l, bool preferMore){
    g_preferMoreCuts = preferMore;
    return solve_lamda(l);
}

int lastBsSteps = 0;
// LAMBDA_SEARCH=bisect: bisection instead of the secant
static const bool g_bisectOnly = [](){ const char* e = getenv("LAMBDA_SEARCH"); return e && string(e) == "bisect"; }();
const char* lastStatus = "exact";
long long lastGapLB = LLONG_MIN;     // lower bound, only for GAP
bool g_recoverPreferMore = true;     // tie-break that gives the reported solution
long long solve_dp_exact(int k){
    lastStatus = "exact";
    lastGapLB = LLONG_MIN;
    g_recoverPreferMore = true;
    if(k > n - 1) k = n - 1;
    if(k <= 0){ Lambda = sumW2 + 1; return (long long)sz[ROOT] * sz[ROOT]; } // 0 cuts

    // secant: keep (l1, c1 > k) and (l2, c2 < k) with U = V - lambda*c,
    // next lambda is the slope of the chord
    int steps = 0;

    steps++;
    pii r0 = eval_lambda(0, true);
    long long l1 = 0, c1 = r0.ss, U1 = r0.ff;
    if(c1 == k){ lastBsSteps = steps; Lambda = 0; g_preferMoreCuts = true; return U1; }
    if(c1 < k){ // fewer than k cuts is already optimal
        lastBsSteps = steps; Lambda = 0; g_preferMoreCuts = true; return U1;
    }

    steps++;
    pii rT = eval_lambda(sumW2 + 1, true); // forces 0 cuts
    long long l2 = sumW2 + 1, c2 = rT.ss, U2 = rT.ff - l2 * c2;
    if(c2 == k){ lastBsSteps = steps; Lambda = l2; g_preferMoreCuts = true; return U2; }

    long long ans = LLONG_MIN;
    while(true){
        long long lm = (!g_bisectOnly && c1 != c2) ? (U2 - U1) / (c1 - c2) : (l1 + l2) / 2;
        if(lm <= l1) lm = l1 + 1;
        if(lm > l2)  lm = l2;
        steps++;
        pii r = eval_lambda(lm, true);
        long long cm = r.ss, Vm = r.ff, Um = Vm - lm * cm;
        if(cm == k){ ans = Um; Lambda = lm; break; }
        bool progress;
        if(cm > k){ progress = (cm < c1) || (lm > l1); l1 = lm; c1 = cm; U1 = Um; }
        else      { progress = (cm > c2) || (lm < l2); l2 = lm; c2 = cm; U2 = Um; }
        if(!progress || l2 - l1 <= 1){
            long long mid = (l1 + l2) / 2;
            if(mid == l1 || mid == l2){
                // l2 = l1 + 1 and still no k. try both tie-breaks at l1 and
                // l2: exactly k cuts means exact. if not, report the best
                // solution with < k cuts (upper bound) and the lower bound,
                // cut counts are not contiguous on trees.
                bool done = false;
                long long ubBest = LLONG_MAX, ubLambda = l1, lbBest = LLONG_MIN;
                bool ubPreferMore = true;
                for(long long lc : {l1, l2}){
                    steps += 2;
                    pii hi = eval_lambda(lc, true);   // max cuts among optima
                    pii lo = eval_lambda(lc, false);  // min cuts among optima
                    if(lo.ss == k || hi.ss == k){
                        ans = hi.ff - lc * (long long)k; // hi.ff == lo.ff (same optimum)
                        Lambda = lc;
                        lastStatus = "exact-collinear";
                        g_recoverPreferMore = (hi.ss == k);
                        done = true;
                        break;
                    }
                    lbBest = max(lbBest, hi.ff - lc * (long long)k);
                    if(hi.ss < k && hi.ff - lc * hi.ss < ubBest){
                        ubBest = hi.ff - lc * hi.ss; ubLambda = lc; ubPreferMore = true;
                    }
                    if(lo.ss < k && lo.ff - lc * lo.ss < ubBest){
                        ubBest = lo.ff - lc * lo.ss; ubLambda = lc; ubPreferMore = false;
                    }
                }
                if(!done && ubBest == lbBest){
                    // upper == lower bound, still exact (common at lambda 0)
                    ans = ubBest;
                    Lambda = ubLambda;
                    g_recoverPreferMore = ubPreferMore;
                    lastStatus = "exact-squeeze";
                    done = true;
                }
                if(!done){
                    ans = ubBest;
                    Lambda = ubLambda;
                    g_recoverPreferMore = ubPreferMore;
                    lastGapLB = lbBest;
                    lastStatus = "GAP_OR_COLLINEAR";
                    fprintf(stderr,
                        "WARNING: lambda search missed k=%lld (possibly an intermediate collinear optimum); "
                        "reported value is a certified <=k-cut solution (upper bound), "
                        "true opt is in [%lld, %lld].\n",
                        (long long)k, (long long)lbBest, (long long)ubBest);
                }
                break;
            }
            steps++;
            pii rb = eval_lambda(mid, true);
            long long cb = rb.ss, Ub = rb.ff - mid * cb;
            if(cb == k){ ans = Ub; Lambda = mid; break; }
            if(cb > k){ l1 = mid; c1 = cb; U1 = Ub; }
            else      { l2 = mid; c2 = cb; U2 = Ub; }
        }
    }
    lastBsSteps = steps;
    g_preferMoreCuts = true;
    return ans;
}

/****** BRUTE (needs root = node 1) ******/
long long check_brute(vector<int>& cut_dad){
    long long tot = 0;
    int nn = T;
    static vector<long long> sub;
    sub.assign(nn + 1, 0);
    for(int t = T; t >= 1; t--){
        int v = here[t];
        sub[v] += val[v];
        if(v == 1 || cut_dad[v]) tot += sub[v] * sub[v];
        else sub[p[v]] += sub[v];
    }
    return tot;
}
long long brute(int k){
    assert(ROOT == 1);
    if(k <= 0) return (long long)sz[1] * sz[1];
    long long best = (long long)sz[1] * sz[1];
    vector<int> cut_dad(n + 1, 0);
    for(long long msk = 0; msk < (1LL << n); msk++){
        if(msk & 1) continue;
        if(__builtin_popcountll(msk) > k) continue;
        rep(i,1,n+1) cut_dad[i] = (msk >> (i-1) & 1);
        cut_dad[1] = 1;
        best = min(best, check_brute(cut_dad));
    }
    return best;
}

int SKIP_BRUTE = 1;

/****** rooting ******/
// centroid by node count
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
    if(argc > 2) SKIP_BRUTE = atoll(argv[2]); // 0 = also run brute
    int rootMode = (argc > 4) ? atoll(argv[4]) : 2; // 2 = auto, 1 = centroid, 0 = node 1
    rep(i,1,n+1) val[i] = readInt();

    {
        vector<vector<int>> adj(n + 1);
        rep(i,2,n+1){ int pi = readInt(); adj[pi].pb(i); adj[i].pb(pi); }
        // the answer doesn't depend on the root, the work does: auto picks
        // node 1 or the centroid, whichever has less light work
        // (sum of chain head subtree sizes - n; 0 on a path rooted at an end)
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
    sumW2 = (long long)sz[ROOT] * sz[ROOT];

    if(!SKIP_BRUTE){
        auto t0 = clock();
        for(int kk=1; kk<n; kk++) cout << brute(kk) << " ";
        cout << "\n[brute time " << (clock()-t0)/(double)CLOCKS_PER_SEC << "s]\n";
    }

    // recovery is part of the timed run: the answer is the cut set, not the cost
    auto t2 = clock();
    long long ans = solve_dp_exact(k);
    double dtValue = (clock()-t2)/(double)CLOCKS_PER_SEC;
    long long searchPeakKb = peak_heap_kb();   // lambda search only

    auto t3 = clock();
    g_preferMoreCuts = g_recoverPreferMore;
    solve_lamda(Lambda);            // refresh memo at lambda*
    recoverCutVert.clear();
    recover(ROOT);
    long long real_cost = cost_recovered(recoverCutVert);
    double dtRecover = (clock()-t3)/(double)CLOCKS_PER_SEC;
    double dt = dtValue + dtRecover;

    assert((long long)recoverCutVert.size() <= max(0LL, (long long)k));
    assert(real_cost == ans);

    cout << ans << "  - k=" << k << ", time: " << dt << "s"
         << ", value_time: " << dtValue << "s, recover_time: " << dtRecover << "s"
         << ", cuts: " << recoverCutVert.size()
         << ", bs_steps: " << lastBsSteps
         << ", lambda*: " << Lambda << ", status: " << lastStatus;
    if(lastGapLB != LLONG_MIN) cout << ", gap_lb: " << lastGapLB;
    cout << ", root: " << ROOT;
    if(searchPeakKb >= 0)
        cout << ", search_peak_kb: " << searchPeakKb << ", peak_kb: " << peak_heap_kb();
    cout << "\n";

    if(argc > 3 && atoll(argv[3]) == 1){
        cout << "cut edge child-vertices:";
        for(int v : recoverCutVert) cout << " " << v;
        cout << "\n";
    }
    return 0;
}
