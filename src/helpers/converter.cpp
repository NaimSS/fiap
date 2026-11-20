// .fi -> tree for the general dp: one node per network node and one per edge
// (with the edge length), rooted at the feeder, parents before children.
// ./converter [K] < input.fi > tree
// output: "n K", n weights, n-1 parents. K is only a default for the solvers,
// n if not given. EMIT_META=1 adds a line with the .fi id of each node (-1 for
// edge nodes).
#include <bits/stdc++.h>
using namespace std;

static bool is_comment_or_empty(const string& s) {
    for (char c : s) {
        if (!isspace((unsigned char)c)) return c == '#';
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

// floor(km * 1e5), read from the string so there's no float error
static long long parse_len_km_to_val100000(const string& s) {
    string a = s;
    if (a.find('e') != string::npos || a.find('E') != string::npos) {
        long double x = stold(a);
        if (x < 0) x = 0;
        return (long long)floor(x * 100000.0L + 1e-18L);
    }

    string intp, frac;
    auto pos = a.find('.');
    if (pos == string::npos) {
        intp = a; frac = "";
    } else {
        intp = a.substr(0, pos);
        frac = a.substr(pos + 1);
    }

    if (intp.empty() || intp == "+") intp = "0";
    if (intp == "-") intp = "0";

    while ((int)frac.size() < 5) frac.push_back('0');
    string frac5 = frac.substr(0, 5);

    long long ip = 0;
    try { ip = stoll(intp); } catch (...) { ip = 0; }
    long long fp5 = 0;
    for (char c : frac5) fp5 = fp5 * 10 + (isdigit((unsigned char)c) ? (c - '0') : 0);

    if (ip < 0) ip = 0;
    return ip * 100000LL + fp5;
}

int main(int argc, char* argv[]) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    long long headerK = -1;  // -1: n
    if (argc > 1) headerK = atoll(argv[1]);

    int declared_nodes = -1;
    int declared_arcs  = -1;
    int feeder_root_id = 1;
    vector<int> feeders;

    unordered_map<int,int> id_to_orig; // .fi id -> 0..n-1
    vector<int> orig_ids;
    auto get_or_add_orig = [&](int id) -> int {
        auto it = id_to_orig.find(id);
        if (it != id_to_orig.end()) return it->second;
        int idx = (int)orig_ids.size();
        id_to_orig[id] = idx;
        orig_ids.push_back(id);
        return idx;
    };

    // edges without repeats, key = (min << 32) | max
    struct E { int u, v; long long len; };
    unordered_map<long long, E> edges_map;

    string line;
    while (std::getline(cin, line)) {
        if (is_comment_or_empty(line)) continue;
        auto tok = split_ws(line);
        if (tok.empty()) continue;

        if (tok[0] == "p") {
            // p fi <nodes> <arcs> <failure rate> <crew velocity>
            if ((int)tok.size() >= 5) {
                declared_nodes = stoi(tok[2]);
                declared_arcs  = stoi(tok[3]);
            }
        } else if (tok[0] == "f") {
            // f <feeder_1> [feeder_2 ...]
            for (int i = 1; i < (int)tok.size(); i++) {
                int fid = stoi(tok[i]);
                feeders.push_back(fid);
                get_or_add_orig(fid);
            }
            if (!feeders.empty()) feeder_root_id = feeders[0];
        } else if (tok[0] == "v") {
            // v <num_node> <x> <y>
            if ((int)tok.size() >= 2) {
                int id = stoi(tok[1]);
                get_or_add_orig(id);
            }
        } else if (tok[0] == "e") {
            // e <src> <tgt> <length in km>
            if ((int)tok.size() >= 4) {
                int a = stoi(tok[1]);
                int b = stoi(tok[2]);
                string ls = tok[3];
                long long len = parse_len_km_to_val100000(ls);

                get_or_add_orig(a);
                get_or_add_orig(b);

                int mn = min(a,b), mx = max(a,b);
                long long key = ( (long long)(uint32_t)mn << 32 ) | (uint32_t)mx;

                auto it = edges_map.find(key);
                if (it == edges_map.end()) {
                    edges_map[key] = E{a,b,len};
                } else if (it->second.len != len) {
                    cerr << "Warning: duplicate edge " << mn << "-" << mx
                         << " with different lengths; keeping the first ("
                         << it->second.len << " vs " << len << ", x1e-5 km)\n";
                }
            }
        }
    }

    // some files have fewer v lines than the p line says: add the missing ids,
    // starting at the smallest id seen (files can be 0- or 1-indexed)
    if (declared_nodes > 0 && (int)orig_ids.size() < declared_nodes) {
        int base = orig_ids.empty() ? 0 : *min_element(orig_ids.begin(), orig_ids.end());
        for (int id = base; id < base + declared_nodes; id++) get_or_add_orig(id);
    }

    vector<E> edges;
    edges.reserve(edges_map.size());
    for (auto &kv : edges_map) edges.push_back(kv.second);

    const int n_orig = (int)orig_ids.size();
    const int m = (int)edges.size();
    const int N = n_orig + m;

    vector<vector<int>> adj(N);
    vector<long long> val(N, 0);

    // edge j becomes node n_orig + j
    for (int j = 0; j < m; j++) {
        int a_id = edges[j].u, b_id = edges[j].v;
        int a = id_to_orig[a_id];
        int b = id_to_orig[b_id];
        int vtx = n_orig + j;
        val[vtx] = edges[j].len;

        adj[a].push_back(vtx);
        adj[vtx].push_back(a);

        adj[b].push_back(vtx);
        adj[vtx].push_back(b);
    }

    // several feeders: weight 0 root above all of them
    int totalN = N;
    int root;
    if ((int)feeders.size() > 1) {
        int superRoot = N;
        adj.push_back({});
        val.push_back(0);
        totalN = N + 1;
        for (int f : feeders) {
            auto it = id_to_orig.find(f);
            if (it == id_to_orig.end()) continue;
            int fi = it->second;
            adj[superRoot].push_back(fi);
            adj[fi].push_back(superRoot);
        }
        root = superRoot;
    } else {
        root = 0;
        auto itroot = id_to_orig.find(feeder_root_id);
        if (itroot != id_to_orig.end()) root = itroot->second;
    }

    // bfs order, so parents come before children
    vector<int> newIndex(totalN, -1);
    vector<int> parent_old(totalN, -1);
    vector<int> order;
    order.reserve(totalN);

    deque<int> q;
    newIndex[root] = 0;
    parent_old[root] = -1;
    q.push_back(root);
    order.push_back(root);

    while (!q.empty()) {
        int u = q.front(); q.pop_front();
        for (int w : adj[u]) {
            if (newIndex[w] != -1) continue;
            newIndex[w] = (int)order.size();
            parent_old[w] = u;
            order.push_back(w);
            q.push_back(w);
        }
    }

    // other components (shouldn't happen) go under the root
    int components = 1;
    for (int i = 0; i < totalN; i++) {
        if (newIndex[i] != -1) continue;
        components++;
        newIndex[i] = (int)order.size();
        parent_old[i] = root;
        order.push_back(i);
        q.push_back(i);
        while (!q.empty()) {
            int u = q.front(); q.pop_front();
            for (int w : adj[u]) {
                if (newIndex[w] != -1) continue;
                newIndex[w] = (int)order.size();
                parent_old[w] = u;
                order.push_back(w);
                q.push_back(w);
            }
        }
    }

    // forest iff edges = n - components (n - 1 alone misses a cycle when
    // some node is disconnected)
    {
        long long edgeCount = 0;
        for (auto& a : adj) edgeCount += (long long)a.size();
        edgeCount /= 2;
        if (edgeCount != (long long)totalN - components) {
            cerr << "Error: graph has a cycle (nodes=" << totalN
                 << ", edges=" << edgeCount << ", components=" << components << ").";
            if ((int)feeders.size() > 1)
                cerr << " Multiple feeders sharing a component create a cycle.";
            cerr << "\n";
            return 1;
        }
        if (components > 1)
            cerr << "Warning: network has " << components << " components; the "
                 << components - 1 << " not containing the root were attached under it\n";
    }

    cout << totalN << " " << (headerK >= 0 ? headerK : (long long)totalN) << "\n";

    for (int i = 0; i < totalN; i++) {
        long long v = val[ order[i] ];
        if (i) cout << ' ';
        cout << v;
    }
    cout << "\n";

    for (int i = 1; i < totalN; i++) {
        int node = order[i];
        int pold = parent_old[node];
        int pnew = (pold == -1 ? 0 : newIndex[pold]);
        int pout = pnew + 1;
        if (i > 1) cout << ' ';
        cout << pout;
    }
    cout << "\n";

    if (getenv("EMIT_META")) {
        for (int i = 0; i < totalN; i++) {
            int node = order[i];
            if (i) cout << ' ';
            cout << (node < n_orig ? orig_ids[node] : -1);
        }
        cout << "\n";
    }

    return 0;
}
