#ifndef _MemoryTracking_hpp_
#define _MemoryTracking_hpp_

#ifdef STK_MEMORY_TRACKING

#include <cstdlib>
#include <cstdio>

size_t total_bytes_allocated = 0;
size_t high_water_mark_bytes = 0;
size_t high_water_mark_in_ptrs = 0;

static const size_t max_num_ptrs = 1000000;
static size_t num_ptrs = 0;
static std::pair<void*,size_t> ptr_sizes[max_num_ptrs];

inline
void add_ptr(void* ptr, size_t sz) {
    if (num_ptrs >= max_num_ptrs) {
        std::printf("too many ptrs allocated!");
        std::abort();
    }
    ptr_sizes[num_ptrs++] = std::make_pair(ptr,sz);
    high_water_mark_in_ptrs = std::max(high_water_mark_in_ptrs, num_ptrs);
}

inline
size_t remove_ptr_and_return_size(void*ptr) {
    int num = num_ptrs;
    for(int i=num-1; i>=0; --i) {
        if (ptr_sizes[i].first == ptr) {
            size_t size = ptr_sizes[i].second;
            ptr_sizes[i] = ptr_sizes[num-1];
            --num_ptrs;
            return size;
        }
    }
    return 0;
}

inline
void* operator new(std::size_t sz) {
    void* ptr = std::malloc(sz);
    add_ptr(ptr,sz);
    total_bytes_allocated += sz;
    high_water_mark_bytes = std::max(high_water_mark_bytes, total_bytes_allocated);
    return ptr;
}

inline
void operator delete(void* ptr) throw()
{
    total_bytes_allocated -= remove_ptr_and_return_size(ptr);
    std::free(ptr);
}

namespace stk {
inline
size_t get_total_bytes_currently_allocated()
{
    return total_bytes_allocated;
}
inline
size_t get_high_water_mark_in_bytes()
{
    return high_water_mark_bytes;
}
inline
size_t get_num_ptrs_currently_allocated()
{
    return num_ptrs;
}
inline
size_t get_high_water_mark_in_ptrs()
{
    return high_water_mark_in_ptrs;
}

}//namespace stk

#endif

#endif

