#pragma once
#include <bits/stdc++.h>
using namespace std;

// piecewise linear f(x) for the tree dp, (value, cuts) on each segment.
// f(x) = S(x + xoff) + lazy_m*x + lazy_b, so noCut is O(1).
// f is always non-decreasing in x (-DPARANOID_CHECKS checks it).

// small vector, most functions have very few segments.
// uses operator new (not malloc) so mem_cap.h counts it.

template <class T, int CAP>
struct SVec {
  T inl[CAP];
  T* d;
  int n, cap;
  SVec() : d(inl), n(0), cap(CAP) {}
  SVec(const SVec& o) : d(inl), n(0), cap(CAP) { assignFrom(o); }
  SVec& operator=(const SVec& o) { assignFrom(o); return *this; }
  ~SVec() { if (d != inl) ::operator delete(d); }
  void grow(int need) {
    if (need <= cap) return;
    int nc = cap * 2; while (nc < need) nc *= 2;
    T* nd = (T*)::operator new(sizeof(T) * nc);
    memcpy(nd, d, sizeof(T) * n);
    if (d != inl) ::operator delete(d);
    d = nd; cap = nc;
  }
  void assignFrom(const SVec& o) { n = 0; grow(o.n); memcpy(d, o.d, sizeof(T) * o.n); n = o.n; }
  void clear() { n = 0; }
  bool empty() const { return n == 0; }
  int size() const { return n; }
  void push_back(const T& v) { if (n == cap) grow(n + 1); d[n++] = v; }
  T& operator[](int i) { return d[i]; }
  const T& operator[](int i) const { return d[i]; }
  T& back() { return d[n - 1]; }
  const T& back() const { return d[n - 1]; }
  T* begin() { return d; }
  T* end() { return d + n; }
  const T* begin() const { return d; }
  const T* end() const { return d + n; }
  void drop_front(int cnt) {
    if (cnt <= 0) return;
    memmove(d, d + cnt, sizeof(T) * (n - cnt));
    n -= cnt;
  }
  void release() { if (d != inl) ::operator delete(d); d = inl; cap = CAP; n = 0; }
};

struct Val {
  long long first; // cost
  int second;      // cuts
};

// on equal cost: more cuts if true, fewer if false
extern bool g_preferMoreCuts;
static inline bool tieA(int ak, int bk) {
  return g_preferMoreCuts ? (ak >= bk) : (ak <= bk);
}

static inline bool betterMin(const Val& a, const Val& b) {
  if (a.first != b.first) return a.first < b.first;
  return g_preferMoreCuts ? (a.second > b.second) : (a.second < b.second);
}
static inline Val pickMin(const Val& a, const Val& b) {
  return betterMin(a, b) ? a : b;
}

static inline long long cdiv(long long a, long long b) { return a/b+((a^b)>0&&a%b); } // divide a by b rounded up
static inline long long fdiv(long long a, long long b) { return a/b-((a^b)<0&&a%b); } // divide a by b rounded down

struct Seg {
  long long l; // start, in y = x + xoff
  long long m; // m*y + b
  long long b;
  int k;       // cuts (+ lazy_k)
};

#ifdef TREE_STATS
struct TreeStats {
  long long mergeCalls = 0, mergeSegWork = 0, constCalls = 0, constWalk = 0;
  long long maxSegs = 0;
  long long segSamples = 0, segSum = 0;     // for the average |F|
  long long recBuilds = 0, recPositions = 0; // build_full calls in recovery
} extern g_treestats;
#define TREE_SAMPLE_SEGS(sz)                                                     \
  do {                                                                         \
    long long sz_ = (long long)(sz);                                        \
    g_treestats.maxSegs = max(g_treestats.maxSegs, sz_);                        \
    g_treestats.segSamples++;                                                    \
    g_treestats.segSum += sz_;                                                \
  } while (0)
#endif

typedef SVec<Seg, 4> SegVec;

struct DPFunc {
  SegVec s; // segment i covers y in [s[i].l, s[i+1].l)

  long long xoff   = 0;
  long long lazy_m = 0;
  long long lazy_b = 0;
  int       lazy_k = 0;

  static constexpr long long INF = (1LL << 60);

  DPFunc() {
    s.push_back(Seg{0, 0, 0, 0}); // f = 0
  }

  void reset() {
    s.clear();
    s.push_back(Seg{0, 0, 0, 0});
    xoff = lazy_m = lazy_b = 0; lazy_k = 0;
  }

  void release() { s.release(); xoff = lazy_m = lazy_b = 0; lazy_k = 0; }

  // last segment with l <= y
  int idxAt(long long y) const {
    int lo = 0, hi = s.size() - 1;
    while (lo < hi) {
      int mid = (lo + hi + 1) >> 1;
      if (s[mid].l <= y) lo = mid; else hi = mid - 1;
    }
    return lo;
  }

  // segment i as m*x + b in real x
  inline void actualCoef(int i, long long& m, long long& b, int& k) const {
    m = s[i].m + lazy_m;
    b = s[i].b + s[i].m * xoff + lazy_b;
    k = s[i].k + lazy_k;
  }

  inline long long xstart(int i) const {
    long long v = s[i].l - xoff;
    return v < 0 ? 0 : v;
  }

  Val eval(long long x) const {
    if (x < 0) return Val{INF, INT_MIN / 4};
    long long y = x + xoff;
    int i = idxAt(y);
    long long m, b; int k;
    actualCoef(i, m, b, k);
    return Val{m * x + b, k};
  }

  Val at0() const { return eval(0); }

  // f(x) <- f(x+w) + 2wx + w^2
  void noCut(long long w) {
    lazy_b += lazy_m * w + w * w;
    lazy_m += 2 * w;
    xoff   += w;
  }

  void addAll(long long value) { lazy_b += value; }
  void addSlope(long long w) { lazy_m += w; }
  void addSecondConst(int delta) { lazy_k += delta; }

  // drop segments that ended up at x < 0
  void compactDead() {
    int i0 = idxAt(xoff);
    if (i0 > 16 && i0 * 2 > s.size()) s.drop_front(i0);
  }

  // f = min(f, C) on [0, xmax]. f is non-decreasing, so C wins a suffix:
  // binary search where it starts and cut the tail
  void mergeMinConst(long long C, int kc, long long xmax = INF) {
#ifdef TREE_STATS
    g_treestats.constCalls++;
    TREE_SAMPLE_SEGS(s.size());
#endif
    compactDead();
    const int n = s.size();
    const int i0 = idxAt(xoff);

    // first segment starting at a value >= C
    int lo = i0, hi = n;
    while (lo < hi) {
      int mid = (lo + hi) >> 1;
      long long m, b; int k;
      actualCoef(mid, m, b, k);
      long long xs = (mid == i0 ? 0 : xstart(mid));
      if (m * xs + b >= C) hi = mid; else lo = mid + 1;
    }
    int start = (lo > i0 ? lo - 1 : i0); // crossing may be inside lo-1

    // drop segments from firstDrop on and put C from x = xc
    auto emitConstTail = [&](int firstDrop, long long xc) {
      Seg ng;
      ng.l = xc + xoff;
      ng.m = -lazy_m;
      ng.b = C - lazy_b + lazy_m * xoff;
      ng.k = kc - lazy_k;
      s.n = firstDrop;
      if (!s.empty() && s.back().l == ng.l) { s.back() = ng; return; }
      if (!s.empty() && s.back().m == ng.m && s.back().b == ng.b && s.back().k == ng.k) return;
      s.push_back(ng);
    };

    for (int i = start; i < n; i++) {
#ifdef TREE_STATS
      g_treestats.constWalk++;
#endif
      long long am, ab; int ak;
      actualCoef(i, am, ab, ak);
      long long L = (i == i0 ? 0 : xstart(i));
      if (L > xmax) return;
      long long R = (i + 1 < n ? xstart(i + 1) : INF);
      long long dm = am;        // d(x) = f(x) - C
      long long db = ab - C;
#ifdef PARANOID_CHECKS
      assert(am >= 0);
#endif
      if (dm == 0) {
        if (db > 0) { emitConstTail(i, L); return; }
        if (db == 0 && !tieA(ak, kc)) { emitConstTail(i, L); return; }
      } else {
        long long split = cdiv(-db, dm); // first x with d(x) >= 0
        if (db % dm == 0) {
          long long xeq = (-db) / dm;
          if (tieA(ak, kc)) split = xeq + 1;
        }
        long long mid = max(L, min(R, split));
        if (mid > xmax) return;
        if (mid < R) {
          if (mid <= L) emitConstTail(i, L);
          else          emitConstTail(i + 1, mid);
          return;
        }
      }
    }
  }

  // f = min(f, other + extra_b) on [0, xmax], other's cuts + extra_k.
  // result has xoff = 0
  void mergeMinWithOffsets(const DPFunc& other, long long extra_b, int extra_k,
                           long long xmax = INF) {
    const int An = s.size();
    const int Bn = other.s.size();

    // first just check if other ever wins, usually it doesn't at big lambda
    {
      int i = idxAt(xoff);
      int j = other.idxAt(other.xoff);
      bool bWins = false;
      while (i < An && j < Bn && !bWins) {
        long long La = (s[i].l - xoff < 0 ? 0 : s[i].l - xoff);
        long long Ra = (i + 1 < An ? s[i + 1].l - xoff : INF);
        long long Lb = (other.s[j].l - other.xoff < 0 ? 0 : other.s[j].l - other.xoff);
        long long Rb = (j + 1 < Bn ? other.s[j + 1].l - other.xoff : INF);
        long long L = max(La, Lb);
        long long R = min(Ra, Rb);
        if (L > xmax) break;
        if (L < R) {
          long long am, ab; int ak;
          actualCoef(i, am, ab, ak);
          long long bm, bb; int bk;
          other.actualCoef(j, bm, bb, bk);
          bb += extra_b; bk += extra_k;
          long long dm = am - bm;       // d = A - B
          long long db = ab - bb;
          long long xt = (dm > 0 ? min(R - 1, xmax) : L);
          long long dmax = dm * xt + db;
          if (dmax > 0) bWins = true;
          else if (dmax == 0 && !tieA(ak, bk)) bWins = true;
        }
        if (Ra < Rb) i++;
        else if (Rb < Ra) j++;
        else { i++; j++; }
      }
      if (!bWins) return;
    }

    int i = idxAt(xoff);
    int j = other.idxAt(other.xoff);

    static thread_local SegVec out;
    out.clear();

    auto emit = [&](long long l, long long r, long long m, long long b, int k) {
      if (l >= r) return;
      if (!out.empty() && out.back().m == m && out.back().b == b && out.back().k == k) return;
      out.push_back(Seg{l, m, b, k});
    };

    while (i < An && j < Bn) {
      long long La = (s[i].l - xoff < 0 ? 0 : s[i].l - xoff);
      long long Ra = (i + 1 < An ? s[i + 1].l - xoff : INF);
      long long Lb = (other.s[j].l - other.xoff < 0 ? 0 : other.s[j].l - other.xoff);
      long long Rb = (j + 1 < Bn ? other.s[j + 1].l - other.xoff : INF);

      long long L = max(La, Lb);
      long long R = min(Ra, Rb);
      if (L > xmax) break;
      if (L < R) {
        long long am, ab; int ak;
        actualCoef(i, am, ab, ak);
        long long bm, bb; int bk;
        other.actualCoef(j, bm, bb, bk);
        bb += extra_b; bk += extra_k;

        long long dm = am - bm;
        long long db = ab - bb;

        if (dm == 0) {
          if (db < 0)      emit(L, R, am, ab, ak);
          else if (db > 0) emit(L, R, bm, bb, bk);
          else {
            if (tieA(ak, bk)) emit(L, R, am, ab, ak);
            else              emit(L, R, bm, bb, bk);
          }
        } else if (dm > 0) {
          // A below on the left, B on the right
          long long split = cdiv(-db, dm);
          if (db % dm == 0) {
            long long xeq = (-db) / dm;
            if (tieA(ak, bk)) split = xeq + 1;
          }
          long long mid = max(L, min(R, split));
          emit(L, mid, am, ab, ak);
          emit(mid, R, bm, bb, bk);
        } else {
          // B below on the left, A on the right
          long long split = fdiv(-db, dm) + 1;
          if (db % dm == 0) {
            long long xeq = (-db) / dm;
            if (tieA(ak, bk)) split = xeq;
          }
          long long mid = max(L, min(R, split));
          emit(L, mid, bm, bb, bk);
          emit(mid, R, am, ab, ak);
        }
      }

      if (Ra < Rb) i++;
      else if (Rb < Ra) j++;
      else { i++; j++; }
    }

    if (out.empty()) out.push_back(Seg{0, 0, 0, 0});
#ifdef TREE_STATS
    g_treestats.mergeCalls++;
    g_treestats.mergeSegWork += (An - idxAt(xoff)) + (Bn - other.idxAt(other.xoff));
    TREE_SAMPLE_SEGS(out.size());
#endif
    s = out;
    xoff = 0; lazy_m = 0; lazy_b = 0; lazy_k = 0;
  }

  void mergeMin(const DPFunc& other) { mergeMinWithOffsets(other, 0, 0); }

  // leaf hanging off the chain, in place:
  //   f(x) <- min(f(x+w) + 2wx + w^2, f(x) + c with kc more cuts)
  // same as G = *this; noCut(w); mergeMinWithOffsets(G, c, kc, xmax)
  void cutLeafMerge(long long w, long long c, int kc, long long xmax = INF) {
    {
      const int n = s.size();
      int i = idxAt(xoff + w), j = idxAt(xoff);
      bool cutWins = false;
      while (i < n && j < n && !cutWins) {
        long long La = (s[i].l - xoff - w < 0 ? 0 : s[i].l - xoff - w);
        long long Ra = (i + 1 < n ? s[i + 1].l - xoff - w : INF);
        long long Lb = (s[j].l - xoff < 0 ? 0 : s[j].l - xoff);
        long long Rb = (j + 1 < n ? s[j + 1].l - xoff : INF);
        long long L = max(La, Lb);
        long long R = min(Ra, Rb);
        if (L > xmax) break;
        if (L < R) {
          // d = keep - cut
          long long dm = (s[i].m - s[j].m) + 2 * w;
          long long db = (s[i].b + s[i].m * (xoff + w) + lazy_m * w + w * w)
                       - (s[j].b + s[j].m * xoff + c);
          long long xt = (dm > 0 ? min(R - 1, xmax) : L);
          long long dmax = dm * xt + db;
          if (dmax > 0) cutWins = true;
          else if (dmax == 0 && !tieA(s[i].k, s[j].k + kc)) cutWins = true;
        }
        if (Ra < Rb) i++;
        else if (Rb < Ra) j++;
        else { i++; j++; }
      }
      if (!cutWins) { noCut(w); return; }
    }

    const int n = s.size();
    int i = idxAt(xoff + w); // keep side
    int j = idxAt(xoff);     // cut side

    static thread_local SegVec out;
    out.clear();

    auto emit = [&](long long l, long long r, long long m, long long b, int k) {
      if (l >= r) return;
      if (!out.empty() && out.back().m == m && out.back().b == b && out.back().k == k) return;
      out.push_back(Seg{l, m, b, k});
    };

    while (i < n && j < n) {
      long long La = (s[i].l - xoff - w < 0 ? 0 : s[i].l - xoff - w);
      long long Ra = (i + 1 < n ? s[i + 1].l - xoff - w : INF);
      long long Lb = (s[j].l - xoff < 0 ? 0 : s[j].l - xoff);
      long long Rb = (j + 1 < n ? s[j + 1].l - xoff : INF);

      long long L = max(La, Lb);
      long long R = min(Ra, Rb);
      if (L > xmax) break;
      if (L < R) {
        long long am = s[i].m + lazy_m + 2 * w;
        long long ab = s[i].b + s[i].m * (xoff + w) + lazy_b + lazy_m * w + w * w;
        int       ak = s[i].k + lazy_k;
        long long bm = s[j].m + lazy_m;
        long long bb = s[j].b + s[j].m * xoff + lazy_b + c;
        int       bk = s[j].k + lazy_k + kc;

        long long dm = am - bm;
        long long db = ab - bb;

        if (dm == 0) {
          if (db < 0)      emit(L, R, am, ab, ak);
          else if (db > 0) emit(L, R, bm, bb, bk);
          else {
            if (tieA(ak, bk)) emit(L, R, am, ab, ak);
            else              emit(L, R, bm, bb, bk);
          }
        } else if (dm > 0) {
          long long split = cdiv(-db, dm);
          if (db % dm == 0) {
            long long xeq = (-db) / dm;
            if (tieA(ak, bk)) split = xeq + 1;
          }
          long long mid = max(L, min(R, split));
          emit(L, mid, am, ab, ak);
          emit(mid, R, bm, bb, bk);
        } else {
          long long split = fdiv(-db, dm) + 1;
          if (db % dm == 0) {
            long long xeq = (-db) / dm;
            if (tieA(ak, bk)) split = xeq;
          }
          long long mid = max(L, min(R, split));
          emit(L, mid, bm, bb, bk);
          emit(mid, R, am, ab, ak);
        }
      }

      if (Ra < Rb) i++;
      else if (Rb < Ra) j++;
      else { i++; j++; }
    }

    if (out.empty()) out.push_back(Seg{0, 0, 0, 0});
#ifdef TREE_STATS
    g_treestats.mergeCalls++;
    g_treestats.mergeSegWork += (n - idxAt(xoff)) + out.size();
    TREE_SAMPLE_SEGS(out.size());
#endif
    s = out;
    xoff = 0; lazy_m = 0; lazy_b = 0; lazy_k = 0;
  }

  void debug() const {
#ifdef LOCAL
    cout << "DPFunc with " << s.size() << " segments (xoff=" << xoff << "):\n";
    for (int i = 0; i < (int)s.size(); i++) {
      long long m, b; int k;
      actualCoef(i, m, b, k);
      cout << "  for x >= " << xstart(i) << ": first(x) = " << m << "*x + " << b
           << ", second(x) = " << k << "\n";
    }
#endif
  }
};
