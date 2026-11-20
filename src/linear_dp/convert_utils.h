#pragma once
// convert_utils.h
// Parses a .fi.diameter file and returns the linearized weight array
// for the DP solvers.
//
// The array alternates: [subtreeWeight(node), edgeWeight, subtreeWeight(node), ...]
// where subtreeWeight collapses all off-diameter branch edge weights into the
// corresponding diameter node. Length = 2*(num_diameter_nodes) - 1.
//
// Edge weights are stored as floor(km * 1e5) integers.

#include <bits/stdc++.h>
using namespace std;

namespace convert_utils {

static bool is_comment_or_empty(const string& s) {
    for (char c : s) {
        if (!isspace((unsigned char)c)) return c == '#' || c == '%';
    }
    return true;
}

static vector<string> split_ws(const string& s) {
    vector<string> t;
    string cur;
    for (char c : s) {
        if (isspace((unsigned char)c)) {
            if (!cur.empty()) { t.push_back(cur); cur.clear(); }
        } else cur.push_back(c);
    }
    if (!cur.empty()) t.push_back(cur);
    return t;
}

// Parse km string, keep 5 fractional digits: value = floor(km * 10^5).
static long long parse_len(const string& s) {
    const string& a = s;
    if (a.find('e') != string::npos || a.find('E') != string::npos) {
        long double x = stold(a);
        if (x < 0) x = 0;
        return (long long)floor(x * 100000.0L + 1e-18L);
    }
    string intp, frac;
    auto pos = a.find('.');
    if (pos == string::npos) { intp = a; frac = ""; }
    else { intp = a.substr(0, pos); frac = a.substr(pos + 1); }
    if (intp.empty() || intp == "+" || intp == "-") intp = "0";
    while ((int)frac.size() < 6) frac.push_back('0');
    string frac5 = frac.substr(0, 5);
    long long ip = 0;
    try { ip = stoll(intp); } catch (...) { ip = 0; }
    long long fp5 = 0;
    for (char c : frac5)
        fp5 = fp5 * 10 + (isdigit((unsigned char)c) ? (c - '0') : 0);
    if (ip < 0) ip = 0;
    return ip * 100000LL + fp5;
}

static bool is_all_digits(const string& s) {
    if (s.empty()) return false;
    for (char c : s) if (!isdigit((unsigned char)c)) return false;
    return true;
}

inline vector<long long> cleanZerosResult(const vector<long long>& result) {
    vector<long long> cleanResult;
    cleanResult.reserve(result.size());
    for (int i = 0; i < (int)result.size(); i++) {
        if (result[i] != 0) cleanResult.push_back(result[i]);
    }
    return cleanResult;
}

// Returns the linearized weight array from a .fi.diameter file.
// Throws std::runtime_error on parse failure.
inline vector<long long> parseLinear(const string& filename) {
    ifstream fin(filename);
    if (!fin) throw runtime_error("Cannot open file: " + filename);

    // ---- Node mapping ----
    unordered_map<int,int> id_to_idx;
    vector<int> orig_ids;
    auto get_or_add = [&](int id) -> int {
        auto it = id_to_idx.find(id);
        if (it != id_to_idx.end()) return it->second;
        int idx = (int)orig_ids.size();
        id_to_idx[id] = idx;
        orig_ids.push_back(id);
        return idx;
    };

    struct Edge { int u, v; long long w; };
    vector<Edge> eEdges;       // indexed by 0-based e-line order
    vector<int> diamArcIdx;    // 0-based e-line indices from the diameter line
    bool readingDiam = false;

    string line;
    while (getline(fin, line)) {
        if (is_comment_or_empty(line)) continue;
        auto tok = split_ws(line);
        if (tok.empty()) continue;

        if (tok[0] == "p" || tok[0] == "f") {
            // skip
        } else if (tok[0] == "v") {
            if ((int)tok.size() >= 2) get_or_add(stoi(tok[1]));
        } else if (tok[0] == "e") {
            if ((int)tok.size() >= 4) {
                int a = stoi(tok[1]), b = stoi(tok[2]);
                long long w = parse_len(tok[3]);
                get_or_add(a); get_or_add(b);
                eEdges.push_back({a, b, w});
            }
        } else if (tok[0] == "diameter" || tok[0] == "DIAMETER") {
            readingDiam = true;
            for (int i = 1; i < (int)tok.size(); i++)
                diamArcIdx.push_back(stoi(tok[i]));
        } else if (readingDiam) {
            bool ok = true;
            for (auto& t : tok) if (!is_all_digits(t)) { ok = false; break; }
            if (ok) for (auto& t : tok) diamArcIdx.push_back(stoi(t));
            else readingDiam = false;
        }
    }

    if (diamArcIdx.empty())
        throw runtime_error("No diameter section found in: " + filename);

    int n = (int)orig_ids.size();

    // ---- Full adjacency list (e-line index as arc ID) ----
    struct AdjE { int to, eidx; long long w; };
    vector<vector<AdjE>> adj(n);
    for (int i = 0; i < (int)eEdges.size(); i++) {
        int u = id_to_idx[eEdges[i].u], v = id_to_idx[eEdges[i].v];
        long long w = eEdges[i].w;
        adj[u].push_back({v, i, w});
        adj[v].push_back({u, i, w});
    }

    // ---- Diameter subgraph (diamArcIdx are 0-based e-line indices) ----
    set<int> diamSet(diamArcIdx.begin(), diamArcIdx.end());
    vector<vector<pair<int,int>>> diamAdj(n); // {neighbor, eidx}
    for (int ei : diamArcIdx) {
        if (ei < 0 || ei >= (int)eEdges.size())
            throw runtime_error("Diameter arc index out of range: " + to_string(ei));
        int u = id_to_idx[eEdges[ei].u], v = id_to_idx[eEdges[ei].v];
        diamAdj[u].push_back({v, ei});
        diamAdj[v].push_back({u, ei});
    }

    // Find endpoints (degree 1 in the diameter path)
    int endA = -1, endB = -1;
    for (int i = 0; i < n; i++) {
        if ((int)diamAdj[i].size() == 1) {
            if (endA == -1) endA = i;
            else endB = i;
        }
    }
    if (endA == -1 || endB == -1)
        throw runtime_error("Could not find diameter endpoints in: " + filename);

    // Walk diameter path from endA to endB
    vector<int> pathNodes;
    vector<long long> pathEdgeW;
    {
        int cur = endA, prev = -1;
        while (true) {
            pathNodes.push_back(cur);
            if (cur == endB) break;
            bool moved = false;
            for (auto& [nb, ei] : diamAdj[cur]) {
                if (nb != prev) {
                    pathEdgeW.push_back(eEdges[ei].w);
                    prev = cur;
                    cur = nb;
                    moved = true;
                    break;
                }
            }
            if (!moved)
                throw runtime_error("Broken diameter path in: " + filename);
        }
    }

    // ---- Collapse: for each diameter node, sum off-diameter subtree edge weights ----
    auto subtreeWeight = [&](int root) -> long long {
        long long total = 0;
        stack<int> stk;
        stk.push(root);
        vector<bool> vis(n, false);
        vis[root] = true;
        for (int nd : pathNodes)
            if (nd != root) vis[nd] = true;
        while (!stk.empty()) {
            int u = stk.top(); stk.pop();
            for (auto& e : adj[u]) {
                if (vis[e.to] || diamSet.count(e.eidx)) continue;
                vis[e.to] = true;
                total += e.w;
                stk.push(e.to);
            }
        }
        return total;
    };

    // ---- Build result: [nodeW, edgeW, nodeW, edgeW, ..., nodeW] ----
    vector<long long> result;
    result.reserve(2 * pathNodes.size() - 1);
    for (int i = 0; i < (int)pathNodes.size(); i++) {
        result.push_back(subtreeWeight(pathNodes[i]));
        if (i < (int)pathEdgeW.size())
            result.push_back(pathEdgeW[i]);
    }


    // Important: We must remove the zeros from the results
    // This isnt just an optimization: If we dont do this,
    // the DP could have the optimal solution with less than K cuts
    // and this could fail the reconstruction step
    result = cleanZerosResult(result);

    return result;
}

// Verifies that a partition described by 0-indexed cut positions produces the
// expected total cost (sum of squared segment sums).
//
// cuts[i] = p means the i-th segment ends at a[p] (0-indexed); the next segment
// starts at a[p+1].  For k segments there must be exactly k-1 cuts.
// Aborts via assert on any inconsistency.
inline void check_solution(int n, const vector<long long>& a,
                           const vector<int>& cuts, long long expected) {
    for (int i = 0; i < (int)cuts.size(); i++) {
        assert(cuts[i] >= 0 && cuts[i] <= n - 2 && "cut index out of range");
        if (i > 0) assert(cuts[i] > cuts[i - 1] && "cuts not strictly increasing");
    }

    long long actual = 0;
    int start = 0;
    for (int c : cuts) {
        long long seg = 0;
        for (int j = start; j <= c; j++) seg += a[j];
        actual += seg * seg;
        start = c + 1;
    }
    long long seg = 0;
    for (int j = start; j < n; j++) seg += a[j];
    actual += seg * seg;

    assert(actual == expected && "partition cost does not match solver output");
}

} // namespace convert_utils
