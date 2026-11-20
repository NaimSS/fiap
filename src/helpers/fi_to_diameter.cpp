// main trunk = longest path from the feeder. writes <input.fi>.diameter, a copy
// of the .fi plus "diameter i0 i1 ..." with the 0-based e-line indices.
// ./fi_to_diameter <input.fi> [source_node_id]
// several feeders: adds a root with 0-length edges to each, also written to the
// output as v / e lines.

#include <bits/stdc++.h>
using namespace std;

static bool is_comment_or_empty(const string& s) {
    for (char c : s)
        if (!isspace((unsigned char)c)) return c == '#' || c == '%' || c == 'c';
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

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <input.fi> [source_node_id]\n";
        return 1;
    }

    string inFile = argv[1];
    int sourceOverride = (argc >= 3) ? stoi(argv[2]) : -1;

    ifstream fin(inFile);
    if (!fin) { cerr << "Error: cannot open " << inFile << "\n"; return 1; }

    vector<string> rawLines;
    string line;
    while (getline(fin, line)) rawLines.push_back(line);
    fin.close();

    // parse
    vector<int> feeders;
    struct Edge { int u, v; double w; };
    vector<Edge> edges;
    set<int> nodeSet;
    map<int, pair<double,double>> coord;

    for (auto& l : rawLines) {
        auto tok = split_ws(l);
        if (tok.empty() || tok[0][0] == 'c' || tok[0][0] == '#' || tok[0][0] == '%') continue;
        if (tok[0] == "f") {
            for (size_t i = 1; i < tok.size(); i++) feeders.push_back(stoi(tok[i]));
        } else if (tok[0] == "v" && tok.size() >= 2) {
            int id = stoi(tok[1]);
            nodeSet.insert(id);
            if (tok.size() >= 4) coord[id] = {stod(tok[2]), stod(tok[3])};
        } else if (tok[0] == "e" && tok.size() >= 4) {
            int u = stoi(tok[1]), v = stoi(tok[2]);
            double w = stod(tok[3]);
            edges.push_back({u, v, w});
            nodeSet.insert(u);
            nodeSet.insert(v);
        }
    }
    int realEdgeCount = (int)edges.size();

    // several feeders: new root (biggest id + 1) with a 0-length edge to each
    int superRoot = -1;
    if ((int)feeders.size() > 1) {
        superRoot = (nodeSet.empty() ? 0 : *nodeSet.rbegin()) + 1;
        for (int f : feeders)
            edges.push_back({superRoot, f, 0.0});
        nodeSet.insert(superRoot);
    }

    // start from the super-root, else the given node, else the feeder
    int source;
    if (superRoot != -1) {
        source = superRoot;
    } else {
        source = (sourceOverride != -1) ? sourceOverride
               : (!feeders.empty())    ? feeders[0]
               : (nodeSet.empty() ? 0 : *nodeSet.begin());
    }
    if (!nodeSet.count(source)) {
        cerr << "Error: source node " << source << " not found in the graph.\n";
        return 1;
    }

    // adj[u] = (neighbor, e-line index, length)
    unordered_map<int, vector<tuple<int,int,double>>> adj;
    for (int i = 0; i < (int)edges.size(); i++) {
        adj[edges[i].u].emplace_back(edges[i].v, i, edges[i].w);
        adj[edges[i].v].emplace_back(edges[i].u, i, edges[i].w);
    }

    // farthest node from the source; bfs is enough since it's a tree
    unordered_map<int, double> dist;
    unordered_map<int, int>    parentNode;
    unordered_map<int, int>    parentEdge; // e-line index

    for (int n : nodeSet) { dist[n] = -1.0; parentNode[n] = -1; parentEdge[n] = -1; }
    dist[source] = 0.0;

    queue<int> q;
    q.push(source);

    int farthest = source;
    double maxDist = 0.0;

    while (!q.empty()) {
        int u = q.front(); q.pop();
        for (auto& [v, eidx, w] : adj[u]) {
            if (dist[v] < 0.0) {
                dist[v] = dist[u] + w;
                parentNode[v] = u;
                parentEdge[v] = eidx;
                q.push(v);
                if (dist[v] > maxDist) { maxDist = dist[v]; farthest = v; }
            }
        }
    }

    // a radial network has no cycle (two feeders in one component make one)
    if ((int)edges.size() > (int)nodeSet.size() - 1) {
        cerr << "Error: graph has a cycle (nodes=" << nodeSet.size()
             << ", edges=" << edges.size() << ").";
        if (superRoot != -1)
            cerr << " Multiple feeders sharing a component create a cycle.";
        cerr << "\n";
        return 1;
    }
    {
        int reached = 0;
        for (int n : nodeSet) if (dist[n] >= 0.0) reached++;
        if (reached != (int)nodeSet.size())
            cerr << "Warning: " << ((int)nodeSet.size() - reached)
                 << " node(s) unreachable from the source are omitted from the "
                    "diameter (a component without a feeder).\n";
    }

    // walk back from the farthest node to the source
    vector<int> pathEdges;
    for (int cur = farthest; parentEdge[cur] != -1; cur = parentNode[cur])
        pathEdges.push_back(parentEdge[cur]);
    reverse(pathEdges.begin(), pathEdges.end());

    // output: the .fi as is, plus the super-root lines and the diameter line
    string outFile = inFile + ".diameter";
    ofstream fout(outFile);
    if (!fout) { cerr << "Error: cannot write " << outFile << "\n"; return 1; }

    int numArtE = (int)edges.size() - realEdgeCount;   // edges to the super-root
    for (auto& l : rawLines) {
        auto tok = split_ws(l);
        if (superRoot != -1 && !tok.empty() && tok[0] == "p" && tok.size() >= 4) {
            // p fi <nodes> <arcs> ...: count the super-root and its edges
            ostringstream os;
            os << tok[0];
            for (size_t i = 1; i < tok.size(); i++) {
                if      (i == 2) os << ' ' << (stoi(tok[2]) + 1);
                else if (i == 3) os << ' ' << (stoi(tok[3]) + numArtE);
                else             os << ' ' << tok[i];
            }
            fout << os.str() << "\n";
        } else {
            fout << l << "\n";
        }
    }

    if (superRoot != -1) {
        // super-root drawn at the mean of the feeder coordinates
        double cx = 0.0, cy = 0.0; int cnt = 0;
        for (int f : feeders) {
            auto it = coord.find(f);
            if (it != coord.end()) { cx += it->second.first; cy += it->second.second; cnt++; }
        }
        if (cnt > 0) { cx /= cnt; cy /= cnt; }
        fout << "v " << superRoot << " " << fixed << setprecision(6)
             << cx << " " << cy << "\n";
        for (int i = realEdgeCount; i < (int)edges.size(); i++)
            fout << "e " << edges[i].u << " " << edges[i].v << " 0\n";
    }

    fout << "diameter";
    for (int idx : pathEdges) fout << " " << idx;
    fout << "\n";
    fout.close();

    cerr << "source=" << source << " farthest=" << farthest
         << " dist=" << maxDist << " km  edges=" << pathEdges.size();
    if (superRoot != -1)
        cerr << "  (super-root " << superRoot << " added for "
             << feeders.size() << " feeders)";
    cerr << "\n";
    cerr << "written: " << outFile << "\n";
    return 0;
}
