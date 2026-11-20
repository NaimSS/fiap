#include <bits/stdc++.h>
using namespace std;
#pragma once

struct Val {
  long long first; // bestValue (MINIMIZE)
  int second;      // numberOfCuts (on ties in first, MAXIMIZE)
};

static inline bool betterMin(const Val& a, const Val& b) {
  if (a.first != b.first) return a.first < b.first;
  return a.second > b.second;
}
static inline Val pickMin(const Val& a, const Val& b) {
  return betterMin(a, b) ? a : b;
}

static inline long long cdiv(long long a, long long b) { return a/b+((a^b)>0&&a%b); } // divide a by b rounded up
static inline long long fdiv(long long a, long long b) { return a/b-((a^b)<0&&a%b); } // divide a by b rounded down

struct Seg {
  long long l; // start x (>=0), inclusive
  long long m; // slope for first: first(x)=m*x+b
  long long b; // intercept
  int k;       // cuts, constant on this segment
   
  // Optional - we could store the id of the taken option for
  // easier recovery of the solution
  int idTaken = -1; // id of the taken cut
};

struct DPFunc {
  // segments sorted by l, covering [0, +inf)
  vector<Seg> s;

  // lazy deltas: actual_m = g.m + lazy_m, actual_b = g.b + lazy_b, actual_k = g.k + lazy_k
  long long lazy_m = 0;
  long long lazy_b = 0;
  int       lazy_k = 0;

  static constexpr long long INF = (1LL << 60);

  DPFunc() {
    // baseline: (0,0) for all x>=0
    s.push_back(Seg{0, 0, 0, 0, -1});
  }
  DPFunc(int id){
    s.push_back(Seg{0, 0, 0, 0, id});
  }

  // apply pending lazy deltas to all segments and reset
  void push() {
    if (lazy_m == 0 && lazy_b == 0 && lazy_k == 0) return;
    for (auto& g : s) { g.m += lazy_m; g.b += lazy_b; g.k += lazy_k; }
    lazy_m = 0; lazy_b = 0; lazy_k = 0;
  }

  // evaluate at integer x>=0
  Val eval(long long x) const {
    if (x < 0) return Val{INF, INT_MIN / 4}; // shouldn't happen (x>=0 by contract)
    int idx = int(upper_bound(s.begin(), s.end(), x,
                              [](long long x, const Seg& seg) { return x < seg.l; }) -
                  s.begin()) -
              1;
    idx = max(idx, 0);
    const Seg& g = s[idx];
    long long f = (g.m + lazy_m) * x + (g.b + lazy_b);
    return Val{f, g.k + lazy_k};
  }

  // coalesce adjacent segments with identical formula
  void normalize() {
    push();
    // Do I really need to sort here?
    // Invariant should be maintained by the other operations.
    // sort(s.begin(), s.end(), [](const Seg& a, const Seg& b) { return a.l < b.l; });

    // drop duplicates by keeping last
    vector<Seg> t;
    for (auto seg : s) {
      if (!t.empty() && t.back().l == seg.l)
        t.back() = seg;
      else
        t.push_back(seg);
    }
    s.swap(t);

    // ensure starts at 0 with baseline if missing/empty
    if (s.empty() || s[0].l != 0) {
      s.insert(s.begin(), Seg{0, 0, 0, 0, -1});
    }

    vector<Seg> out;
    out.reserve(s.size());
    for (auto seg : s) {
      if (seg.l < 0) continue;
      if (out.empty()) {
        out.push_back(seg);
        continue;
      }
      auto& p = out.back();
      if (p.m == seg.m && p.b == seg.b && p.k == seg.k) {
        // same formula, redundant breakpoint
        continue;
      }
      out.push_back(seg);
    }
    s.swap(out);
  }

  // Add on all: f(x) <- f(x) + value
  void addAll(long long value) { lazy_b += value; }

  // Add slope: f(x) <- f(x) + x*w
  void addSlope(long long w) { lazy_m += w; }

  // second (cuts): f(x).second += delta
  void addSecondConst(int delta) { lazy_k += delta; }

  void reset(int id = -1) {
    s.clear();
    s.push_back(Seg{0, 0, 0, 0, id});
    lazy_m = 0; lazy_b = 0; lazy_k = 0;
  }

  // Shift: f(x) <- f(x+w), w >= 0
  void shift(long long w) {
    if (w <= 0) return;
    // lazy_m * w must be folded into lazy_b before the stored-b update below
    lazy_b += lazy_m * w;

    vector<Seg> out;
    out.reserve(s.size());
    for (int i = 0; i < (int)s.size(); i++) {
      long long l = s[i].l;
      long long r = (i + 1 < (int)s.size() ? s[i + 1].l : INF);
      long long nl = l - w;
      long long nr = r - w;
      if (nr <= 0) continue; // entirely shifted into x<0

      // old: m*(x+w)+b => m*x + (b + m*w)
      Seg ng;
      ng.l = max(0LL, nl);
      ng.m = s[i].m;
      ng.b = s[i].b + s[i].m * w;
      ng.k = s[i].k;
      ng.idTaken = s[i].idTaken;
      out.push_back(ng);
    }
    if (out.empty()) out.push_back(Seg{0, 0, 0, 0, -1}); // baseline fallback
    s.swap(out);
    normalize();
  }

  // Merge: f(x) <- min(f(x), other(x) + extra_b), with other's k offset by extra_k
  // Tie-break: if first equal, pick LARGEST second.
  void mergeMinWithOffsets(const DPFunc& other, long long extra_b, int extra_k) {
    push(); // materialize this->s in-place

    const long long om = other.lazy_m;
    const long long ob = other.lazy_b + extra_b;
    const int       ok = other.lazy_k + extra_k;

    const int An = (int)s.size();
    const int Bn = (int)other.s.size();

    vector<Seg> out;
    out.reserve(An + Bn);

    auto emit = [&](long long l, long long r, long long m, long long b, int k, int id) {
      if (l >= r) return;
      if (!out.empty() && out.back().m == m && out.back().b == b && out.back().k == k) return;
      out.push_back(Seg{l, m, b, k, id});
    };

    int i = 0, j = 0;
    while (i < An && j < Bn) {
      long long La = s[i].l;
      long long Ra = (i + 1 < An ? s[i + 1].l : INF);
      long long Lb = other.s[j].l;
      long long Rb = (j + 1 < Bn ? other.s[j + 1].l : INF);

      long long L = max(La, Lb);
      long long R = min(Ra, Rb);
      if (L < R) {
        long long am = s[i].m,           ab = s[i].b,           ak = s[i].k;
        long long bm = other.s[j].m + om, bb = other.s[j].b + ob;
        int       bk = other.s[j].k + ok;

        int idA = s[i].idTaken;
        int idB = other.s[j].idTaken;

        long long dm = am - bm;
        long long db = ab - bb; // d(x)=A-B

        if (dm == 0) {
          if (db < 0) {
            // A < B everywhere on [L,R)
            emit(L, R, am, ab, ak, idA);
          } else if (db > 0) {
            // B < A everywhere on [L,R)
            emit(L, R, bm, bb, bk, idB);
          } else {
            // equal first -> pick larger second
            if (ak >= bk) emit(L, R, am, ab, ak, idA);
            else          emit(L, R, bm, bb, bk, idB);
          }
        } else if (dm > 0) {
          // d increasing: A below B on the left, B wins on the right
          long long split = cdiv(-db, dm); // first x with d(x) >= 0
          if (db % dm == 0) {
            long long xeq = (-db) / dm;
            if (ak >= bk) split = xeq + 1; // integer tie at xeq -> A's side
          }
          long long mid = max(L, min(R, split));
          emit(L, mid, am, ab, ak, idA);
          emit(mid, R, bm, bb, bk, idB);
        } else {
          // d decreasing: B below A on the left, A wins on the right
          long long split = fdiv(-db, dm) + 1;
          if (db % dm == 0) {
            long long xeq = (-db) / dm;
            if (ak >= bk) split = xeq; // integer tie at xeq -> A from xeq onward
          }
          long long mid = max(L, min(R, split));
          emit(L, mid, bm, bb, bk, idB);
          emit(mid, R, am, ab, ak, idA);
        }
      }

      if (Ra < Rb) i++;
      else if (Rb < Ra) j++;
      else { i++; j++; }
    }

    if (out.empty()) out.push_back(Seg{0, 0, 0, 0, -1});
    s.swap(out);
    lazy_m = 0; lazy_b = 0; lazy_k = 0; // merged segments are actual (no lazy on this)
  }

  void mergeMin(const DPFunc& other) {
    mergeMinWithOffsets(other, 0, 0);
  }

  // // Optional helper mirroring your old "keep zero baseline" behavior:
  // // result(x) = min(result(x), 0) with tie going to larger k.
  // void mergeMinAndKeepZero(const DPFunc& other) {
  //   mergeMin(other);
  //   DPFunc base;     // (0,0) everywhere
  //   mergeMin(base);  // clamp with baseline (0,0)
  // }

  void debug() const {
    #ifdef LOCAL
      cout << "DPFunc with " << s.size() << " segments:\n";
      for (const auto& g : s) {
        cout << "  for x >= " << g.l << ": first(x) = " << (g.m + lazy_m) << "*x + " << (g.b + lazy_b)
            << ", second(x) = " << (g.k + lazy_k) << "\n";
      }
    #endif
  }

  Val at0() const { 
    return eval(0); 
  }
};
