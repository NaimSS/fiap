#include <bits/stdc++.h>
#include "convert_utils.h"
#define rep(i,a,b) for(int i=(a);i<(b);++i)
using namespace std;
typedef long long ll;
typedef pair<ll,ll> pll;

ll sq(ll x) { return x*x; }

// Monotone CHT for min queries with decreasing-slope insertion and
// increasing-x queries. Each line represents f(x) = slope*x + cte.
//
// max_tie controls tie-breaking on segment count (qtd):
//   max_tie = true  -> prefer the line with higher qtd (simulates +\eps*qtd)
//   max_tie = false -> prefer the line with lower  qtd (simulates -\eps*qtd)
struct line {
    ll slope, cte;
    int qtd, pos;   // qtd = segment count; pos = originating index (1-based)
    line(ll slope, ll cte, int qtd, int pos)
        : slope(slope), cte(cte), qtd(qtd), pos(pos) {}
    ll eval(ll x) const { return slope * x + cte; }
};

struct cht {
    int p = 0;
    bool max_tie;
    vector<line> func;

    cht(bool max_tie = true) : max_tie(max_tie) {}

    // Returns true if B is never uniquely optimal between A and C.
    // Uses __int128 to prevent overflow: ctes can be O(k*\lambda + sum weights squared), and
    // the cross-product of slopes/intercepts would overflow 64-bit at large lambda.
    bool bad(line A, line B, line C) {
        __int128 cross = (__int128)(B.cte - A.cte) * (A.slope - C.slope)
                       - (__int128)(C.cte - A.cte) * (A.slope - B.slope);
        if (cross != 0) return cross > 0;
        // Lines A and B cross C at the same x: B has zero-width optimal range.
        // Keep B only if it wins the qtd tie-break over A.
        return max_tie ? B.qtd <= A.qtd : B.qtd >= A.qtd;
    }

    void addline(line r) {
        if (!func.empty() && func.back().slope == r.slope) {
            if (func.back().cte <= r.cte) return;   // r is dominated
            func.pop_back();
        }
        while ((int)func.size() >= 2 &&
               bad(func[func.size()-2], func[func.size()-1], r))
            func.pop_back();
        if (p >= (int)func.size()) p = max(0, (int)func.size() - 1);
        func.push_back(r);
    }

    // Returns {min value at x, qtd of winner, pos of winner}.
    tuple<ll,int,int> queryMin(ll x) {
        assert(!func.empty() && p < (int)func.size());
        while (p < (int)func.size() - 1 &&
               (func[p].eval(x) > func[p+1].eval(x) ||
                (func[p].eval(x) == func[p+1].eval(x) &&
                 (max_tie ? func[p+1].qtd > func[p].qtd
                           : func[p+1].qtd < func[p].qtd)))) p++;
        return {func[p].eval(x), func[p].qtd, func[p].pos};
    }
};

// Penalized DP at fixed lambda (inner step of WQS / Aliens binary search).
//
// Segment cost: C(i,j)^2 where C(i,j) = pref[j-1] - pref[i-1].
// Each segment also incurs the penalty lambda, so we minimise \sum C^2 + \lambda*segs.
//
// Recurrence (right-to-left over 1-indexed positions):
//   dp[n+1] = 0
//   dp[i]   = min_{j>i} { C(i,j)^2 + dp[j] } + \lambda
//           = min_j { (C(1,j) - C(1,i))^2 + dp[j] } + \lambda
//           = min_j { C(1,j)^2 - 2*C(1,i)*C(1,j) + dp[j] } + C(1,i)^2 + \lambda
//
// Viewing each j as a line with slope = C(1,j) and intercept = dp[j]+C(1,j)^2,
// the inner min becomes a CHT query at x = -2*C(1,i).
//
// dp[i] stores {penalized cost, segment count} with max-tie on count.
// nxt_out[i] records the chosen j (start of the next segment after i).
void solve_lambda(vector<ll>& pref, int n, vector<pll>& dp, ll lambda,
                  vector<int>& nxt) {
    auto C = [&](int i, int j) -> ll {
        return i == j ? 0 : pref[j-1] - pref[i-1];
    };

    cht CHT; // Convex hull trick - ties for largest number of segments
    for (int i = n+1; i >= 1; i--) {
        if (i > n) {
            dp[i] = {0, 0};
        } else {
            auto [vval, vqtd, vpos] = CHT.queryMin(-2 * C(1, i));
            dp[i] = {vval + sq(C(1,i)) + lambda, vqtd + 1};
            nxt[i] = vpos;
        }
        CHT.addline({C(1,i), dp[i].first + sq(C(1,i)), (int)dp[i].second, i});
    }
}

// Reconstruct a k-segment solution at the critical lambda.
//
// Runs two CHTs in parallel:
//   max-tie CHT -> dp_cnt_max[i]: maximum reachable segment count from i
//   min-tie CHT -> dp_cnt_min[i]: minimum reachable segment count from i
//
// At the critical lambda the true k satisfies dp_cnt_min[1] <= k <= dp_cnt_max[1].
//
// Greedy reconstruction: with `rem` segments remaining at position `pos`,
// try the max-tie transition j_max. If the remaining count rem-1 is still
// achievable from j_max (i.e. in [dp_cnt_min[j_max], dp_cnt_max[j_max]]),
// take it; otherwise take the min-tie transition j_min.
//
// Returns {true optimal cost for k segments, cut positions (0-indexed)}.
pair<ll, vector<int>> solve_final(vector<ll>& pref, int n, ll lambda, int k) {
    auto C = [&](int i, int j) -> ll {
        return i == j ? 0 : pref[j-1] - pref[i-1];
    };

    vector<ll> dp_val(n+2, 0);
    vector<int> dp_cnt_max(n+2, 0), dp_cnt_min(n+2, 0);
    vector<int> nxt_max(n+2, -1), nxt_min(n+2, -1);

    cht CHT_max(true), CHT_min(false);

    for (int i = n+1; i >= 1; i--) {
        if (i > n) {
            dp_val[i] = dp_cnt_max[i] = dp_cnt_min[i] = 0;
        } else {
            auto [vval_max, vqtd_max, vpos_max] = CHT_max.queryMin(-2 * C(1,i));
            auto [vval_min, vqtd_min, vpos_min] = CHT_min.queryMin(-2 * C(1,i));
            dp_val[i]     = vval_max + sq(C(1,i)) + lambda;
            dp_cnt_max[i] = vqtd_max + 1;
            dp_cnt_min[i] = vqtd_min + 1;
            nxt_max[i]    = vpos_max;
            nxt_min[i]    = vpos_min;
        }
        ll base = dp_val[i] + sq(C(1,i));
        CHT_max.addline({C(1,i), base, dp_cnt_max[i], i});
        CHT_min.addline({C(1,i), base, dp_cnt_min[i], i});
    }

    vector<int> cuts;
    int pos = 1;
    for (int rem = k; rem > 1; rem--) {
        int j_max = nxt_max[pos], j_min = nxt_min[pos];
        int j = (dp_cnt_min[j_max] <= rem-1 && rem-1 <= dp_cnt_max[j_max])
                ? j_max : j_min;
        cuts.push_back(j - 2);
        pos = j;
    }
    return {dp_val[1] - lambda * k, cuts};
}

// Binary search over lambda for lagrange optimization
// Binary search: find the largest lambda where the max-tie DP gives >= k segs.
// At that lambda, f(k) = dp[1] - \lambda*k.
// Set when the bisection over lambda never lands on exactly k cuts, i.e. f has a
// flat segment at the target k and the answer comes from solve_final. Used by the
// scan tool that reports which real instances hit this corner case.
bool g_flat_case = false;

pair<ll, vector<int>> solve_lagrange(vector<ll>& a, int k) {
    int n = (int)a.size();
    g_flat_case = false;

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

    g_flat_case = true;
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
