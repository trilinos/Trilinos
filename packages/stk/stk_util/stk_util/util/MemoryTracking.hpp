#ifndef _MemoryTracking_hpp_
#define _MemoryTracking_hpp_

#ifdef STK_MEMORY_TRACKING

#include <cstdlib>

namespace stk {

size_t get_total_bytes_currently_allocated();

size_t get_high_water_mark_in_bytes();

size_t get_high_water_mark_in_ptrs();

}//namespace stk

#endif

#endif

