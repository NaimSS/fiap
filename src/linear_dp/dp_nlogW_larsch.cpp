#include <bits/stdc++.h>
#include "convert_utils.h"
#include "larsch.h"
#define rep(i,a,b) for(int i=(a);i<(b);++i)
using namespace std;
typedef long long ll;
typedef pair<ll,ll> pll;

ll sq(ll x) { return x*x; }

const ll INF_LAR = (ll)4e18;

// Penalized DP at fixed lambda, solved with LARSCH (O(n) per lambda).
//
// Recurrence (right-to-left, 1-indexed):
//   dp[n+1] = 0
//   dp[i]   = min_{j>i} { (pref[j-1]-pref[i-1])^2 + dp[j] } + lambda
//
// Segment count for dp[1] is recovered by tracing nxt from position 1.
// tie_dir = +1 maximizes segment count among equally-optimal paths (max-tie);
// tie_dir = -1 minimizes it (min-tie). Encoding: minimize
//     SCALE * (real cost) - tie_dir * cnt[j]
// with SCALE = n+2 so the integer tie term cannot flip the primary order.
// Subtracting a function of only J preserves the Monge property, so LARSCH
// still applies. __int128 avoids overflow from the SCALE multiplication.
void solve_lambda(vector<ll>& pref, int n, vector<pll>& dp, ll lambda,
                  vector<int>& nxt, int tie_dir = +1) {
    auto w = [&](int i, int j) -> ll {
        ll d = pref[j-1] - pref[i-1];
        return d * d;
    };

    vector<ll> dp_val(n + 2, 0);
    vector<int> cnt(n + 2, 0);
    nxt.assign(n + 2, -1);

    const __int128 SCALE = (__int128)(n + 2);
    const __int128 INF128 = (__int128)1 << 100;

    larsch<__int128> la(n, [&](int r, int c) -> __int128 {
        if (c > r) return INF128;
        int i = n - r;
        int j = n + 1 - c;
        ll real = w(i, j) + dp_val[j] + lambda;
        return SCALE * (__int128)real - (__int128)tie_dir * (__int128)cnt[j];
    });

    for (int r = 0; r < n; r++) {
        int c = la.get_argmin();
        int i = n - r;
        int j = n + 1 - c;
        dp_val[i] = w(i, j) + dp_val[j] + lambda;
        cnt[i]    = cnt[j] + 1;
        nxt[i]    = j;
    }

    dp[1] = {dp_val[1], cnt[1]};
}

// Reconstruct a k-segment solution at the critical lambda using lemma 2.1
// of https://www.osti.gov/servlets/purl/10146169.
//
// Run solve_lambda at the critical `lambda` with max-tie and min-tie
// respectively, giving two optimal paths with segment counts hi >= k >= lo.
// The lemma takes the smaller-count path as i-side and the larger-count path
// as j-side (j is denser). Find x in [0, lo-1] such that
//   b_lo[x] < b_hi[u] < b_hi[u+1] <= b_lo[x+1]
// with u = x + k - lo (where 0 <= u, u+1 <= hi). Splice: first u edges of the
// hi path, bridge edge b_hi[u] -> b_lo[x+1], then suffix of the lo path.
// Total = u + 1 + (lo - x - 1) = k edges. Cost equals f(k) exactly.
//
// Returns {true optimal cost for k segments, cut positions (0-indexed)}.
pair<ll, vector<int>> solve_final(vector<ll>& pref, int n, ll lambda, int k) {
    auto w = [&](int i, int j) -> ll {
        ll d = pref[j-1] - pref[i-1];
        return d * d;
    };

    vector<pll> dp_s(n+2), dp_t(n+2);
    vector<int> nxt_s, nxt_t;
    solve_lambda(pref, n, dp_s, lambda, nxt_s, +1); // max-tie
    solve_lambda(pref, n, dp_t, lambda, nxt_t, -1); // min-tie

    int hi = (int)dp_s[1].second;
    int lo = (int)dp_t[1].second;
    assert(hi >= k && k >= lo);

    auto recover = [&](vector<int>& nxt, int cnt) {
        vector<int> b; b.reserve(cnt + 1);
        int pos = 1;
        for (int i = 0; i < cnt; i++) { b.push_back(pos); pos = nxt[pos]; }
        b.push_back(n + 1);
        return b;
    };
    vector<int> b_hi = recover(nxt_s, hi);
    vector<int> b_lo = recover(nxt_t, lo);

    auto emit = [&](const vector<int>& b, int segs) -> pair<ll, vector<int>> {
        ll val = 0;
        for (int i = 0; i < segs; i++) val += w(b[i], b[i+1]);
        vector<int> cuts;
        for (int i = 1; i < segs; i++) cuts.push_back(b[i] - 2);
        return {val, cuts};
    };
    if (hi == k) return emit(b_hi, hi);
    if (lo == k) return emit(b_lo, lo);

    int x_found = -1;
    for (int x = 0; x < lo; x++) {
        int u = x + k - lo;
        if (u < 0 || u + 1 > hi) continue;
        if (b_lo[x] < b_hi[u] && b_hi[u+1] <= b_lo[x+1]) { x_found = x; break; }
    }
    assert(x_found >= 0);
    int x = x_found, u = x + k - lo;

    vector<int> bk;
    bk.reserve(k + 1);
    for (int i = 0; i <= u; i++) bk.push_back(b_hi[i]);
    for (int i = x + 1; i <= lo; i++) bk.push_back(b_lo[i]);

    ll val = 0;
    for (int i = 0; i < k; i++) val += w(bk[i], bk[i+1]);

    vector<int> cuts;
    for (int i = 1; i < k; i++) cuts.push_back(bk[i] - 2);

    return {val, cuts};
}

// Binary search over lambda for lagrange optimization
// Binary search: find the largest lambda where the max-tie DP gives >= k segs.
// At that lambda, f(k) = dp[1] - \lambda*k.
pair<ll, vector<int>> solve_lagrange(vector<ll>& a, int k) {
    int n = (int)a.size();

    vector<ll> pref(n+2, 0);
    rep(i, 1, n+1) pref[i] = pref[i-1] + a[i-1];

    vector<pll> dp(n+2);
    vector<int> nxt(n+2);

    ll l = 0, r = sq(pref[n]) + 1;
    int found = 0;
    while (l < r) {
        ll mid = (l + r + 1) / 2;
        solve_lambda(pref, n, dp, mid, nxt);
        if(dp[1].second == k){
            l = mid;
            found = 1;
            break;
        }
        if (dp[1].second >= k){
            l = mid;
        }else{
            r = mid - 1;
        }                   
    }

    if(found){
        // easy reconstruct: exactly K cuts
        vector<int> cuts;
        int pos = 1;
        for (int rem = k; rem > 1; rem--) {
            int j = nxt[pos];
            cuts.push_back(j - 2); // cuts 0-indexed
            pos = j;
        }
        return {dp[1].first - l * k, cuts};
    }

    return solve_final(pref, n, l, k);
}

pair<ll, vector<int>> solve(int n, int k, vector<ll>& a) {
    k = min(k, n);
    return solve_lagrange(a, k);
}

int main(int argc, char* argv[]) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    if (argc >= 3) {
        // Benchmark mode: ./dp_nlogW <file.fi.diameter> <K>
        string filename = argv[1];
        int k = atoi(argv[2]);

        vector<ll> a;
        try {
            a = convert_utils::parseLinear(filename);
        } catch (const exception& e) {
            cerr << e.what() << "\n";
            return 1;
        }
        int n = (int)a.size();

        auto t0 = chrono::high_resolution_clock::now();
        auto [val, cuts] = solve(n, k, a);
        auto t1 = chrono::high_resolution_clock::now();
        double ms = chrono::duration<double, milli>(t1 - t0).count();

        cout << "k = " << k << ": " << val << "\n";
        cout << fixed << setprecision(3) << "time_ms: " << ms << "\n";
        convert_utils::check_solution(n, a, cuts, val);
    } else {
        // Stdin mode: echo "n\na[0] ... a[n-1]" | ./dp_nlogW <K>
        int k = (argc >= 2) ? atoi(argv[1]) : 0;
        int n;
        cin >> n;
        vector<ll> a(n);
        rep(i, 0, n) cin >> a[i];
        auto [val, cuts] = solve(n, k, a);
        cout << "k = " << k << ": " << val << "\n";
        convert_utils::check_solution(n, a, cuts, val);
    }
}
