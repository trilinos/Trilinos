#include <stk_util/util/heap_usage.hpp>
#include <stk_util/util/FeatureTest.hpp>

#if defined(__GNUC__)
#ifndef __APPLE__
#include <malloc.h>
#else
#include <sys/malloc.h>
#endif

#elif defined(__PGI)
#include <malloc.h>

#elif defined(__sun)
#include <fstream>
#include <procfs.h>

#elif defined(__SUNPRO_CC)
#include <sys/resource.h>
#endif

namespace stk
{

// Similar to Platform.cpp's get_heap_info except that it doesn't provide largest_free
// and it doesn't have log output.
size_t get_heap_used()
{
  size_t heap_size = 0;

#if defined(SIERRA_HEAP_INFO)
 
# if defined(SIERRA_PTMALLOC3_ALLOCATOR) || defined(SIERRA_PTMALLOC2_ALLOCATOR)
  heap_size = malloc_used();
  
# elif ( defined(__linux__) || defined(REDS) ) && ! defined(__IBMCPP__)
  static struct mallinfo minfo;
  minfo = mallinfo();
  heap_size = static_cast<unsigned int>(minfo.uordblks) + static_cast<unsigned int>(minfo.hblkhd);

# elif defined(__sun)
  pstatus_t proc_status;

  std::ifstream proc("/proc/self/status", std::ios_base::in|std::ios_base::binary);
  if (proc) {
    proc.read(reinterpret_cast<char *>(&proc_status), sizeof(proc_status));
    heap_size = proc_status.pr_brksize;
  }
# endif
#endif // defined(SIERRA_HEAP_INFO)
  return heap_size;
}
} // namespace stk
