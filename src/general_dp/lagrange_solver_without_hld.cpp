#include <bits/stdc++.h>
#include "slope_min_helper_without_hld.h"
#define int long long
#include "debug.h"
#define rep(i,a,b) for(int i=(a);i<(b);++i)
#define pb push_back
#define ff first
#define ss second
using namespace std;
typedef pair<int,int> pii;


const int N = 8000;
int val[N];
int p[N];
int T=0;
int tin[N],tout[N];
int sz[N];
vector<int> g[N];
int here[N];

void prep(int v,int p=-1){
    tin[v] = ++T;
    here[T] = v;
    sz[v] = val[v];
    for(auto to : g[v])if(to != p){
        prep(to, v);
        sz[v] += sz[to];
    }
    tout[v] = T;
}

int prefSum[N];
int Lambda;

// per-node dynamic arrays (allocated in init_dp)
vector<vector<int>> nxt_arr;      // nxt_arr[v][tin[to]-tin[v]]
vector<vector<int>> szAfter_arr;  // szAfter_arr[v][tin[to]-tin[v]]
vector<pii> dp;                   // dp[v]
vector<vector<DPFunc>> dp2;       // dp2[v][t-tin[v]]
int curVisDP = 0;
vector<int> visDP;                // visDP[v], lazy clear via curVisDP
vector<vector<int>> vis2;         // vis2[v][t-tin[v]], lazy clear via curVisDP

int in_sub(int v,int to){
    return (tin[v] <= tin[to] && tout[to] <= tout[v]);
}

// access helpers for nxt/szAfter
int& nxt_ref(int v, int to) { return nxt_arr[v][tin[to] - tin[v]]; }
int& szAfter_ref(int v, int to) { return szAfter_arr[v][tin[to] - tin[v]]; }

void calc(int v,int to, int debug = 0){
    int &nv = nxt_ref(v, to);
    if(nv != 0) return;
    nv = -1;
    for(int t = tin[to] + 1; t <= T && t <= tout[v]; t++){
        if(in_sub(to, here[t]))continue;
        nv = t;
        break;
    }
    if(debug){
        cout << "calc for " << v<< ", " << to << " = " << nv << endl;
    }
}

int sz_after(int v,int to){
    calc(v, to);
    int &x = szAfter_ref(v, to);
    if(x != -1) return x;
    int nv = nxt_ref(v, to);
    if(nv == -1) return x = 0;
    x = prefSum[tout[v]] - prefSum[nv - 1];
    return x;
}

void init_dp(){
    int nn = T;
    nxt_arr.assign(nn + 1, {});
    szAfter_arr.assign(nn + 1, {});
    dp.assign(nn + 1, pii(0, 0));
    dp2.assign(nn + 1, {});
    visDP.assign(nn + 1, 0);
    vis2.assign(nn + 1, {});
    for(int v = 1; v <= nn; v++){
        int sub_sz = tout[v] - tin[v] + 1;
        nxt_arr[v].assign(sub_sz, 0);
        szAfter_arr[v].assign(sub_sz, -1);
        dp2[v].resize(sub_sz + 1); // t offsets 0..sub_sz (includes tout[v]+1)
        vis2[v].assign(sub_sz + 1, 0);
    }
}

const DPFunc emptyDP = DPFunc();

pii solve(int v);


// Option 1: dont pick
// dp[v][at][x] = dp[v][at+1][x + w] + (x+w)*(x+w)-x^2 =
// dp[v][at+1][x+w] + 2xw + w^2
// So dp[v][at](x) = dp[x][at+1](x + w) + 2xw + w^2

// Option 2: cut
// dp[v][at][x] = dp[v][nxt[v][at]][x] + dp[at][at] + lambda

// dp[v][at][x] = min(option1, option2)

void solve2(int v,int t) {
    if(t == -1)return ;
    int to = t - tin[v];
    if (vis2[v][to] == curVisDP) return;
    vis2[v][to] = curVisDP;

    dp2[v][to].reset();
    if( t > tout[v]) return ;

    int at = here[t];
    dbg("solve2", v, t, at, val[at]);
    // nao cortar aqui
    solve2(v, t+1);
    int to1 = t + 1 - tin[v];
    dp2[v][to] = dp2[v][to1];
    DPFunc& noCut = dp2[v][to];
    noCut.shift(val[at]); // x_new = x_old - val[at]
    noCut.addAll(val[at] * val[at]); // f(x) <- f(x) + w^2
    noCut.addSlope(2 * val[at]); // f(x) <- f(x) + x*2*w
    dbg("processing no cut at", v, t, at);
    noCut.debug();
    // cortar
    calc(v, at);
    if(at != v){
        auto cortar = solve(at);
        int nv = nxt_ref(v, at);
        if(!(nv > t) && nv != -1) {
            cout << "Error: nxt[" << v << "][" << at << "] = " << nv << " not > " << t << endl;
            exit(1);
        }
        solve2(v, nv);
        dbg("processing cut at", v, t, at, "cortar =", cortar.ff, cortar.ss);
        if(nv == -1) {
            DPFunc tmp;
            noCut.mergeMinWithOffsets(tmp, Lambda + cortar.ff, 1 + cortar.ss);
        } else {
            int nv_off = nv - tin[v];
            noCut.mergeMinWithOffsets(dp2[v][nv_off], Lambda + cortar.ff, 1 + cortar.ss);
        }
    }
    return ;
}

pii solve(int v){
    if(visDP[v] == curVisDP) return dp[v];
    visDP[v] = curVisDP;
    solve2(v, tin[v]);
    auto bestP = dp2[v][0].at0(); // t=tin[v] → offset 0
    dp[v] = pii(bestP.first, bestP.second);
    return dp[v];
}

void clear(){
    // lazy clear :)
    curVisDP++;
}

pii solve_lamda(int l) {
    clear();
    Lambda = l;
    auto r = solve(1);
    return r;
}


vector<int> recoverCutVert; // must cut from v -> p[v]

void recover(int v);

void recover2(int v,int t,long long w){
    dbg(v,t,w);
    if(t == -1) return;
    if(t > tout[v]) return ;
    int to = t - tin[v];
    int at = here[t];
    auto best = dp2[v][to].eval( w );
    long long cost = best.first;
    long long qtd = best.second;
    // test no cut
    auto bestDpTo1 = dp2[v][to+1].eval( w + val[at] );
    long long costNoCut = bestDpTo1.first + val[at] * val[at] + 2 * w * val[at];

    if(costNoCut == cost && bestDpTo1.second == qtd) {
        // didnt cut here then
        return recover2(v, t + 1, w + val[at]);
    }
    int nv = nxt_ref(v, at);
    int nv_off = nv - tin[v];
    {
        auto cortar = solve(at);
\
        // eval at same w because we didnt pick AT so continue to next subtree
        // with the same weight
        auto best = (nv == -1 ? Val{0 , 0} :  dp2[v][nv_off].eval( w ));
        long long costHere = cortar.first + best.first + Lambda;
        long long k_here = 1 + cortar.second + best.second;

        if(costHere == cost && k_here == qtd) {
            // cut here
            recoverCutVert.pb(at);
            recover(at);
            if(nv == -1)return;
            return recover2(v, nv, w);
        }
    }
    // impossible :(
}

void recover(int v){
    recover2(v, tin[v], 0);
}

long long cost_recovered(vector<int> cut){
    int n = T;
    vector<int> isCut(n + 1);
    for(int v : cut)isCut[v] = 1;

    for(int v : cut){
        dbg("cut", v);
    }

    long long res = 0;
    auto dfs = [&](int v,int p,auto self) -> long long {
        long long w = val[v];
        for(auto to : g[v]){
            if(to == p)continue;
            w += self(to, v, self);
        }

        if( p == -1 || isCut[v] ){
            dbg("add", v, p, w);
            res += w * w;
            return 0;
        }
        return w;
    };

    dfs(1, -1, dfs);
    return res;
}

void check_recovered_solution(int k,int lambda){
    recoverCutVert.clear();
    k = min(k, T-1);
    if(Lambda != lambda ){
        clear();
        Lambda = lambda;
    }
    pii cost = solve(1);
    recover(1);
    cost.ff -= k * lambda;

    long long recovered_cost = cost_recovered(recoverCutVert);
    assert(recoverCutVert.size() == k);
    dbg(cost.ff,cost.ss,recovered_cost);
    assert(recovered_cost == cost.first);
}

int maxPointsInDP = 0;
void compute_num_points(){
    int points = 0;
    rep(v,1,T+1){
        int sub_sz = tout[v] - tin[v] + 1;
        rep(ti,0,sub_sz+1){
            if(vis2[v][ti] == curVisDP)
                points += dp2[v][ti].s.size();
        }
    }
    maxPointsInDP = max(maxPointsInDP, points);
}

int n;
int resp_brute[N];
int brute(int k);
long long check_brute(vector<int> cut_dad);
int bruteK1();
int bruteK2();
void testSmallK(int k);

map<int,int> resp_dp;
int lastBsSteps = 0;
int solve_dp(int k){
    if(k > n - 1) k = n - 1; // no need for larger k
    if(resp_dp.count(k))return resp_dp[k];
    int found = 0;
    int bsSteps = 0;
    pii val;
    long long l = 0, r = sz[1] * sz[1]; //r = 2 * sz[1] * sz[1];
    while( l < r ){
        long long mid = (l + r + 1) / 2;
        bsSteps++;
        val = solve_lamda(mid);
        compute_num_points();
        if(val.ss == k){
            l = mid;
            found = 1;
            break;
        }
        if(val.ss >= k){
            l = mid;
        } else {
            r = mid - 1;
        }
    }
    if(!found) {
        return -1;
    }
    int cost = val.ff - Lambda * k;
    if(val.ss != k){
        cout << "Warning: final cuts = " << val.ss << " != k = " << k << endl;
        assert(0);
    }
    check_recovered_solution(k, Lambda);
    compute_num_points();
    lastBsSteps = bsSteps;
    return resp_dp[k] = cost;
}

const int SKIP_BRUTE = 1;
const int BRUTE_SMALL_K = 0;

#undef int
int main(int argc, char* argv[]){
#define int long long
    int k;
    cin >> n >> k;
    if(argc > 1){
        k = atoi(argv[1]);
    }
    rep(i,1,n+1){
        cin >> val[i];
    }
    rep(i,2,n+1){
        cin >> p[i];
        g[p[i]].pb(i);
    }
    prep(1);

    init_dp();
    rep(i,1,n+1){
        prefSum[i] = prefSum[i-1] + val[here[i]];
    }

    auto start = clock();
    if(!SKIP_BRUTE){
        rep(k,1,n)cout << brute(k) << " ";
        cout << endl;
        cout << "Brute time: " << (clock() - start) / (double) CLOCKS_PER_SEC << "s\n";
    }
    start = clock();
    if(BRUTE_SMALL_K || !SKIP_BRUTE){
        cerr << "Processing N = " << n << " with brute force checks\n";
    }
    int kk = k;
    for(int k : {kk}){
        if( BRUTE_SMALL_K ){
            if(k > 2)break;
            testSmallK(k);
            continue;
        }
        auto t2 = clock();
        cout << solve_dp(k) << " ";
        if(SKIP_BRUTE){

            cout << " - k=" << k << ", ";
            cout << "time: " << (clock() - t2) / (double) CLOCKS_PER_SEC << "s";
            cout << ", bs_steps: " << lastBsSteps << "\n";
            cout.flush();
        }
    }
    cout << endl;
    cout << "Max points in DP: " << maxPointsInDP << endl;
    return 0;
    cout << "DP time: " << (clock() - start) / (double) CLOCKS_PER_SEC << "s\n";
    if(!SKIP_BRUTE){
        vector<int> diff;
        for(int k=1;k<n;k++){
            if(brute(k) != solve_dp(k)){
                cout << "Mismatch at k = " << k << ": brute = " << brute(k) << ", dp = " << solve_dp(k) << endl;
                assert(false);
                return 0;
            }
            diff.push_back(brute(k-1) - brute(k));
            cout << diff.back() << " ";
        }
        cout << endl;
        rep(i,1,diff.size()){
            if(diff[i] > diff[i-1]){
                cout << "weird " << i << endl;
                assert(false);
            }
        }
    }
}

/****** BRUTE FORCE STUFF *******/

void testSmallK(int k){
    auto tStart = clock();
    auto ans = (k == 1 ? bruteK1() : bruteK2());
    auto tEnd = clock();
    cerr << "Brute time for k=" << k << " : " << (tEnd - tStart) / (double) CLOCKS_PER_SEC << "s\n";
    tStart = clock();
    auto dpAns = solve_dp(k);
    tEnd = clock();
    cerr << "dpTime = " << (tEnd - tStart) / (double) CLOCKS_PER_SEC << "s for k = " << k << endl;
    if(ans != solve_dp(k)){
        cerr << "Mismatch at k = " << k << ": brute = " << ans << ", dp = " << solve_dp(k) << endl;
        assert(false);
    }
    cerr << solve_dp(k) << " " << ans << endl;
}

int brute(int k){
    if(resp_brute[k] != 0)return resp_brute[k];
    int best = sz[1] * sz[1];
    vector<int> cut_dad(n+1);

    int tot = 0;
    auto paint = [&](int v,int p, auto self) -> int {
        int Sz = val[v];
        for(auto to : g[v])if(to != p){
            int nxt = self(to, v, self);
            if(!cut_dad[to]) Sz += nxt;
        }
        if(cut_dad[v])tot += Sz * Sz;
        return Sz;
    };

    for(int msk = 1; msk < (1<<n);msk++){
        if(msk&1)continue;
        if(__builtin_popcount(msk) > k)continue;
        rep(i,1,n+1)cut_dad[i] = (msk >> (i-1) & 1);
        cut_dad[1] = 1;
        tot = 0;
        paint(1,-1,paint);
        best = min(best, tot);
    }
    resp_brute[k] = best;
    return best;
}
long long check_brute(vector<int> cut_dad){
    int tot = 0;
    auto paint = [&](int v,int p, auto self) -> int {
        int Sz = val[v];
        for(auto to : g[v])if(to != p){
            int nxt = self(to, v, self);
            if(!cut_dad[to]) Sz += nxt;
        }
        if(cut_dad[v])tot += Sz * Sz;
        return Sz;
    };
    paint(1,-1,paint);
    return tot;
}

void testSmallK(int k);

int bruteK1() {
    vector<int> cut_dad(n+1);
    int best = sz[1] * sz[1];
    for(int i=1;i<=n;i++){
        rep(j,1,n+1)cut_dad[j] = 0;
        cut_dad[i] = 1;
        cut_dad[1] = 1;
        int tot = check_brute(cut_dad);
        best = min(best, tot);
    }
    return best;
}
int bruteK2(){
    int best = sz[1] * sz[1];
    rep(i,1,n+1)rep(j,i+1,n+1){
        vector<int> cut_dad(n+1);
        cut_dad[i] = 1;
        cut_dad[j] = 1;
        cut_dad[1] = 1;
        int tot = check_brute(cut_dad);
        best = min(best, tot);
    }
    return best;
}
