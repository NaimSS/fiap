#include <bits/stdc++.h>
#include "convert_utils.h"
#define rep(i,a,b) for(int i=(a);i<(b);++i)
using namespace std;
typedef long long ll;

ll sq(ll x){ return x*x; }

const ll inf = 1e18;

// convex hull trick structure for minimum queries
// https://github.com/tdas0/lib/blob/master/library/data_structure/CHT_2.cpp
struct line_cht {
    ll slope, intercept;
    int pos;  // 1-indexed position j this line represents
};

struct cht {
    int p = 0;
    vector<line_cht> func;

    bool bad(line_cht A, line_cht B, line_cht C) {
        return (__int128)(B.intercept - A.intercept) * (A.slope - C.slope)
             > (__int128)(C.intercept - A.intercept) * (A.slope - B.slope);
    }
    ll f(line_cht r, ll x) { return r.slope * x + r.intercept; }

    void addline(line_cht r) {
        if (!func.empty() && func.back().slope == r.slope) {
            if (func.back().intercept <= r.intercept) return;
            func.pop_back();
        }
        while ((int)func.size() >= 2 && bad(func[func.size()-2], func[func.size()-1], r)) {
            func.pop_back();
        }
        if (p >= (int)func.size()) p = max(0, (int)func.size() - 1);
        func.push_back(r);
    }

    // returns {min value, pos of the winning line}
    pair<ll, int> queryMin(ll x) {
        assert(!func.empty());
        assert(p < (int)func.size());
        while (p < (int)func.size() - 1 && f(func[p], x) > f(func[p + 1], x)) p++;
        return {f(func[p], x), func[p].pos};
    }
};

pair<ll, vector<int>> solve(int n, int k, vector<ll>& a) {
    k = min(k, n);
    assert((ll)(n + 2) * (k + 1) <= 100'000'000LL);

    vector<ll> pref(n + 1);
    rep(i, 1, n + 1) pref[i] = pref[i - 1] + a[i - 1];

    auto C = [&](int i, int j) -> ll {
        if (i == j) return 0;
        return pref[j - 1] - pref[i - 1];
    };

    // dp[i][kk] = min cost to partition a[i..n] into kk segments
    vector<vector<ll>> dp(n + 2, vector<ll>(k + 1, inf));
    // nxt[i][kk] = 1-indexed start of next segment in optimal split of dp[i][kk]
    vector<vector<int>> nxt(n + 2, vector<int>(k + 1, -1));
    for (int kk = 0; kk <= k; kk++) dp[n + 1][kk] = 0;

    vector<cht> CHT(k + 1);

    for (int i = n + 1; i >= 1; i--) {
        if (i > n) {
            // seed CHT for each level with the n+1 base case
            for (int kk = 1; kk <= k; kk++) {
                line_cht l = { C(1, i), dp[i][kk] + sq(C(1, i)), i };
                CHT[kk].addline(l);
            }
            continue;
        }
        dp[i][0] = inf;
        if (k >= 1) dp[i][1] = sq(C(i, n + 1));
        for (int kk = 2; kk <= k; kk++) {
            auto [val, j] = CHT[kk - 1].queryMin(-2 * C(1, i));
            dp[i][kk] = val + sq(C(1, i));
            nxt[i][kk] = j;  // j = 1-indexed start of next segment
        }
        for (int kk = 1; kk <= k; kk++) {
            line_cht l = { C(1, i), dp[i][kk] + sq(C(1, i)), i };
            CHT[kk].addline(l);
        }
    }
    
    // Reconstruct 0-indexed cut positions.
    // nxt[pos][rem] = j means next segment starts at 1-indexed j,
    // so the 0-indexed cut (last element of this segment) is j-2.
    vector<int> cuts;
    int pos = 1;
    for (int rem = k; rem > 1; rem--) {
        int j = nxt[pos][rem];
        cuts.push_back(j - 2);  // 0-indexed cut position
        pos = j;
    }
    return {dp[1][k], cuts};
}

int main(int argc, char* argv[]) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    if (argc >= 3) {
        // Benchmark mode: ./dp_nk <file.fi.diameter> <K>
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

        // cerr << "n = " << n << ", k = " << k << "\n";

        auto t0 = chrono::high_resolution_clock::now();
        auto [val, cuts] = solve(n, k, a);
        auto t1 = chrono::high_resolution_clock::now();
        double ms = chrono::duration<double, milli>(t1 - t0).count();

        cout << "k = " << k << ": " << val << "\n";
        cout << fixed << setprecision(3) << "time_ms: " << ms << "\n";
        convert_utils::check_solution(n, a, cuts, val);
    } else {
        // Stdin mode: echo "n\na[0] ... a[n-1]" | ./dp_nk <K>
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
