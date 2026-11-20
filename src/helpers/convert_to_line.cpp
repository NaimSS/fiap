// .fi.diameter -> array for the linear dps. every branch off the trunk is
// folded into its trunk node, so the array alternates
// node (weight of its branches), edge, node, edge, ..., node
// ./convert_to_line < input.fi.diameter
// output: size, then the weights. EMIT_META=1 adds a line with the node id /
// e-line index of each entry.

#include <bits/stdc++.h>
using namespace std;

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

// floor(km * 1e5), same as converter.cpp
static long long parse_len(const string& s) {
    string a = s;
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

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    // .fi id -> 0..n-1
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

    vector<int> feederNodes;
    struct RawEdge { int u, v; long long w; };
    vector<RawEdge> eEdges;       // in file order
    vector<int> diamArcIdx;
    bool readingDiam = false;     // the diameter list can span several lines

    string line;
    while (getline(cin, line)) {
        if (is_comment_or_empty(line)) continue;
        auto tok = split_ws(line);
        if (tok.empty()) continue;

        if (tok[0] == "p") {
        } else if (tok[0] == "f") {
            for (int i = 1; i < (int)tok.size(); i++)
                feederNodes.push_back(stoi(tok[i]));
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
            // continuation line, only numbers
            bool ok = true;
            for (auto& t : tok) if (!is_all_digits(t)) { ok = false; break; }
            if (ok) for (auto& t : tok) diamArcIdx.push_back(stoi(t));
            else readingDiam = false;
        }
    }

    if (diamArcIdx.empty()) {
        cerr << "Error: no diameter section found.\n";
        return 1;
    }

    // arc i = i-th e line, the indices used by fi_to_diameter
    int n = (int)orig_ids.size();

    struct Arc { int u, v; long long w; };
    vector<Arc> allArcs;
    allArcs.reserve(eEdges.size());

    for (auto& e : eEdges)
        allArcs.push_back({id_to_idx[e.u], id_to_idx[e.v], e.w});

    struct AdjE { int to, arcIdx; long long w; };
    vector<vector<AdjE>> adj(n);
    for (int i = 0; i < (int)allArcs.size(); i++) {
        int u = allArcs[i].u, v = allArcs[i].v;
        long long w = allArcs[i].w;
        adj[u].push_back({v, i, w});
        adj[v].push_back({u, i, w});
    }

    // graph with only the trunk arcs
    set<int> diamArcSet(diamArcIdx.begin(), diamArcIdx.end());

    vector<vector<pair<int,int>>> diamAdj(n); // {neighbor, arcIdx}
    for (int ai : diamArcIdx) {
        if (ai < 0 || ai >= (int)allArcs.size()) {
            cerr << "Error: diameter arc index " << ai << " out of range.\n";
            return 1;
        }
        int u = allArcs[ai].u, v = allArcs[ai].v;
        diamAdj[u].push_back({v, ai});
        diamAdj[v].push_back({u, ai});
    }

    // ends of the trunk have degree 1
    int endA = -1, endB = -1;
    for (int i = 0; i < n; i++) {
        if ((int)diamAdj[i].size() == 1) {
            if (endA == -1) endA = i;
            else endB = i;
        }
    }
    if (endA == -1 || endB == -1) {
        cerr << "Error: could not find diameter endpoints.\n";
        return 1;
    }

    // walk the trunk from endA to endB
    vector<int> pathNodes;
    vector<long long> pathEdgeW;
    vector<int> pathArcEidx;   // e-line index of each trunk arc
    {
        int cur = endA, prev = -1;
        while (true) {
            pathNodes.push_back(cur);
            if (cur == endB) break;
            bool moved = false;
            for (auto& [nb, ai] : diamAdj[cur]) {
                if (nb != prev) {
                    pathEdgeW.push_back(allArcs[ai].w);
                    pathArcEidx.push_back(ai);
                    prev = cur;
                    cur = nb;
                    moved = true;
                    break;
                }
            }
            if (!moved) {
                cerr << "Error: broken diameter path.\n";
                return 1;
            }
        }
    }

    // weight of the branches hanging off each trunk node
    vector<bool> isDiamNode(n, false);
    for (int nd : pathNodes) isDiamNode[nd] = true;

    auto subtreeWeight = [&](int root) -> long long {
        long long total = 0;
        stack<int> stk;
        stk.push(root);
        vector<bool> vis(n, false);
        vis[root] = true;
        // don't walk into the other trunk nodes
        for (int nd : pathNodes)
            if (nd != root) vis[nd] = true;
        while (!stk.empty()) {
            int u = stk.top(); stk.pop();
            for (auto& e : adj[u]) {
                if (vis[e.to]) continue;
                if (diamArcSet.count(e.arcIdx)) continue;
                vis[e.to] = true;
                total += e.w;
                stk.push(e.to);
            }
        }
        return total;
    };

    // node, edge, node, ..., node
    vector<long long> result;
    result.reserve(2 * pathNodes.size() - 1);
    for (int i = 0; i < (int)pathNodes.size(); i++) {
        result.push_back(subtreeWeight(pathNodes[i]));
        if (i < (int)pathEdgeW.size())
            result.push_back(pathEdgeW[i]);
    }

    cout << (int)result.size() << "\n";
    for (int i = 0; i < (int)result.size(); i++) {
        if (i) cout << ' ';
        cout << result[i];
    }
    cout << "\n";

    // EMIT_META: .fi node id at even positions, e-line index at odd ones
    if (getenv("EMIT_META")) {
        bool first = true;
        for (int i = 0; i < (int)pathNodes.size(); i++) {
            if (!first) cout << ' ';
            first = false;
            cout << orig_ids[pathNodes[i]];
            if (i < (int)pathArcEidx.size())
                cout << ' ' << pathArcEidx[i];
        }
        cout << "\n";
    }

    return 0;
}
