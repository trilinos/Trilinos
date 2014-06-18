#ifndef STK_HEAP_USAGE_H
#define STK_HEAP_USAGE_H

#include <stddef.h>

namespace stk {

  // Similar to Platform.cpp's get_heap_info except that it doesn't provide largest_free
  // and it doesn't have log output.
  size_t get_heap_used();
}

#endif
