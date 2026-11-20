#include <bits/stdc++.h>
#include "convert_utils.h"
#define rep(i,a,b) for(int i=(a);i<(b);++i)
using namespace std;
typedef long long ll;
#include "swamk.h"

ll sq(ll x){ return x*x; }

const ll inf = 1e18;

pair<ll, vector<int>> solve(int n, int k, vector<ll>& a) {
    k = min(k, n);
    // avoid memory overflow
    assert((ll)(n + 2) * (k + 1) <= 100'000'000LL);

    vector<ll> pref(n + 1); // 1-index
    rep(i, 1, n + 1) pref[i] = pref[i - 1] + a[i - 1];

    auto C = [&](int i, int j) -> ll {
        if (i >= j){
            // symmetric reflection to keep monge property
            return pref[i-1] - pref[j-1];
        }
        return pref[j - 1] - pref[i - 1];
    };
    auto w = [&](int i,int j) -> ll {
        return sq(C(i,j));
    };

    // dp[i][kk] = min cost to partition a[i..n] into kk segments
    // 1-index
    vector<vector<ll>> dp(n + 2, vector<ll>(k + 1, inf));
    // nxt[i][kk] = 1-indexed start of next segment in optimal split of dp[i][kk]
    vector<vector<int>> nxt(n + 2, vector<int>(k + 1, -1));
    
    // Fill first col
    for(int i=1;i<=n+1;i++){
        dp[i][1] = w(i, n+1);
        dp[i][0] = inf;
    }
    dp[n+1][0] = 0;
    dp[n+1][1] = inf;

    for(int kk=2;kk <= k;kk++){
        // we must now use SMAWK to get the row-minima of
        // M[i][j] = w(i,j) + dp(j,kk-1) if i < j else \infty
        // select returns where M[i][j1] > M[i][j2]
        auto select = [&](int i,int j1,int j2){
            // make 1-index
            i++,j1++,j2++;
            return w(i,j1) + dp[j1][kk-1] > w(i,j2) + dp[j2][kk-1];
        };
        auto row_minima = smawk(n+1, n+1, select);
        // now we know that dp[i][kk] will use j = row_minima[i] for the transition
        for(int i=1;i<=n+1;i++){
            if((n+1)-i < kk){
                dp[i][kk] = inf;
                continue;
            }
            int j = row_minima[i-1] + 1;
            nxt[i][kk] = j;
            dp[i][kk] = w(i, j) + dp[j][kk-1];
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
