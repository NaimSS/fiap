#include <bits/stdc++.h>
#include "convert_utils.h"
#define rep(i,a,b) for(int i=(a);i<(b);++i)
using namespace std;
typedef long long ll;

ll sq(ll x){ return x*x; }

const ll inf = 1e18;

pair<ll, vector<int>> solve(int n, int k, vector<ll>& a) {
    k = min(k, n);
    assert((ll)(n + 2) * (k + 1) <= 100'000'000LL);

    vector<ll> pref(n + 1);
    rep(i, 1, n + 1) pref[i] = pref[i - 1] + a[i - 1];

    auto C = [&](int i, int j) { return pref[j] - pref[i - 1]; };

    // dp[i][kk] = min cost to partition a[i..n] into kk segments
    vector<vector<ll>> dp(n + 2, vector<ll>(k + 1, inf));
    // nxt_table[i][kk] = 1-indexed start of the next segment in the optimal split
    vector<vector<int>> nxt_table(n + 2, vector<int>(k + 1, -1));
    for (int kk = 0; kk <= k; kk++) dp[n + 1][kk] = 0;

    for (int i = n + 1; i >= 1; i--) {
        if (i > n) continue;
        dp[i][0] = inf;
        dp[i][1] = sq(C(i, n));
        nxt_table[i][1] = n + 1;  // sentinel: segment covers [i..n], no next
        for (int kk = 2; kk <= k; kk++) {
            rep(nxt, i, n) {
                ll cost = dp[nxt + 1][kk - 1] + sq(C(i, nxt));
                if (cost < dp[i][kk]) {
                    dp[i][kk] = cost;
                    nxt_table[i][kk] = nxt + 1;  // 1-indexed start of next segment
                }
            }
        }
    }

    // Reconstruct 0-indexed cut positions.
    // nxt_table[pos][rem] = j means segment ends just before 1-indexed j,
    // so the 0-indexed cut (last element of this segment) is j-2.
    vector<int> cuts;
    int pos = 1;
    for (int rem = k; rem > 1; rem--) {
        int j = nxt_table[pos][rem];
        cuts.push_back(j - 2);  // 0-indexed cut position
        pos = j;
    }
    return {dp[1][k], cuts};
}

int main(int argc, char* argv[]) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    if (argc >= 3) {
        // Benchmark mode: ./dp_n2k <file.fi.diameter> <K>
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
        if (getenv("PRINT_CUTS")) {
            cout << "cuts:";
            for (int c : cuts) cout << ' ' << c;
            cout << "\n";
        }
    } else {
        // Stdin mode: echo "n\na[0] ... a[n-1]" | ./dp_n2k <K>
        int k = (argc >= 2) ? atoi(argv[1]) : 0;
        int n;
        cin >> n;
        vector<ll> a(n);
        rep(i, 0, n) cin >> a[i];
        auto [val, cuts] = solve(n, k, a);
        cout << "k = " << k << ": " << val << "\n";
        convert_utils::check_solution(n, a, cuts, val);
        if (getenv("PRINT_CUTS")) {
            cout << "cuts:";
            for (int c : cuts) cout << ' ' << c;
            cout << "\n";
        }
    }
}
