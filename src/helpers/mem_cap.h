// heap cap inside the process (ulimit doesn't limit memory on macOS):
// operator new/delete count the bytes and exit(1) past the cap.
//   g++ -O2 -include src/helpers/mem_cap.h dp.cpp -o dp
//   MEM_CAP_MB=3072 ./dp   (default 3072, 0 = no cap)
//   MEM_CAP_REPORT=1 ./dp  prints mem_cap_peak_kb at exit
// counts requested bytes, so real RSS is a bit higher.
// only sees operator new: anything allocated with malloc is invisible.
#ifndef MEM_CAP_H
#define MEM_CAP_H

#include <new>
#include <cstdlib>
#include <cstdio>
#include <cstddef>
#include <cstdint>

namespace mem_cap_detail {
static std::size_t g_cap = 0, g_used = 0, g_peak = 0;
static bool g_capSet = false;

inline void account(std::size_t sz) {
    if (!g_capSet) {
        g_capSet = true;
        const char* e = std::getenv("MEM_CAP_MB");
        long mb = e ? std::atol(e) : 3072;
        g_cap = mb > 0 ? ((std::size_t)mb << 20) : 0;
    }
    if (g_cap && g_used + sz > g_cap) {
        std::fprintf(stderr,
            "mem_cap: memory cap exceeded (%zu MB used + %zu MB requested "
            "> %zu MB; raise MEM_CAP_MB to allow more)\n",
            g_used >> 20, sz >> 20, g_cap >> 20);
        std::exit(1);
    }
    g_used += sz;
    if (g_used > g_peak) g_peak = g_used;
}

inline void* trackedAlloc(std::size_t sz) {
    account(sz);
    void* raw = std::malloc(sz + 16);             // size in a 16 byte header
    if (!raw) { std::fprintf(stderr, "mem_cap: out of memory\n"); std::exit(1); }
    *(std::size_t*)raw = sz;
    return (char*)raw + 16;
}

inline void trackedFree(void* p) {
    if (!p) return;
    void* raw = (char*)p - 16;
    g_used -= *(std::size_t*)raw;
    std::free(raw);
}

// alignas(>16) types use operator new(size_t, align_val_t), count those too.
// layout: raw ... [pad][raw pointer][size] | aligned block
inline void* trackedAllocAligned(std::size_t sz, std::size_t al) {
    if (al < 16) al = 16;
    account(sz);
    void* raw = std::malloc(sz + al + 16);
    if (!raw) { std::fprintf(stderr, "mem_cap: out of memory\n"); std::exit(1); }
    std::uintptr_t p = ((std::uintptr_t)raw + 16 + al - 1) & ~(std::uintptr_t)(al - 1);
    ((void**)p)[-2] = raw;
    ((std::size_t*)p)[-1] = sz;
    return (void*)p;
}

inline void trackedFreeAligned(void* p) {
    if (!p) return;
    g_used -= ((std::size_t*)p)[-1];
    std::free(((void**)p)[-2]);
}

// runs at exit, also on the cap's exit(1), but not on a timeout or abort()
inline void reportPeak() {
    const char* r = std::getenv("MEM_CAP_REPORT");
    if (r && std::atoi(r)) std::fprintf(stderr, "mem_cap_peak_kb: %zu\n", g_peak >> 10);
}
struct Reporter { Reporter() { std::atexit(reportPeak); } };
static Reporter g_reporter;
} // namespace mem_cap_detail

// the peak only grows, so reading it between two phases gives the peak of the
// first one. use under #ifdef MEM_CAP_H so builds without the cap compile
inline std::size_t mem_cap_peak_bytes() { return mem_cap_detail::g_peak; }
inline std::size_t mem_cap_used_bytes() { return mem_cap_detail::g_used; }

void* operator new(std::size_t sz) { return mem_cap_detail::trackedAlloc(sz); }
void* operator new[](std::size_t sz) { return mem_cap_detail::trackedAlloc(sz); }
void operator delete(void* p) noexcept { mem_cap_detail::trackedFree(p); }
void operator delete[](void* p) noexcept { mem_cap_detail::trackedFree(p); }
void operator delete(void* p, std::size_t) noexcept { mem_cap_detail::trackedFree(p); }
void operator delete[](void* p, std::size_t) noexcept { mem_cap_detail::trackedFree(p); }

void* operator new(std::size_t sz, std::align_val_t al) {
    return mem_cap_detail::trackedAllocAligned(sz, (std::size_t)al);
}
void* operator new[](std::size_t sz, std::align_val_t al) {
    return mem_cap_detail::trackedAllocAligned(sz, (std::size_t)al);
}
void operator delete(void* p, std::align_val_t) noexcept { mem_cap_detail::trackedFreeAligned(p); }
void operator delete[](void* p, std::align_val_t) noexcept { mem_cap_detail::trackedFreeAligned(p); }
void operator delete(void* p, std::size_t, std::align_val_t) noexcept { mem_cap_detail::trackedFreeAligned(p); }
void operator delete[](void* p, std::size_t, std::align_val_t) noexcept { mem_cap_detail::trackedFreeAligned(p); }

#endif // MEM_CAP_H
