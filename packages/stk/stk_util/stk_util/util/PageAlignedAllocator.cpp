
#include <sys/mman.h>
#include <unistd.h>
#include <cstdlib>

#include <stk_util/util/PageAlignedAllocator.hpp>

namespace stk { namespace detail {

namespace {

  inline bool use_mmap(size_t num_bytes)
  {
    static const size_t half_page_size = (sysconf( _SC_PAGE_SIZE ) >> 1);

    return (num_bytes >= half_page_size);
  }

}

void * page_aligned_allocate(size_t num_bytes)
{

  void * ptr = NULL;

  if ( use_mmap(num_bytes) ) {
    const int flags = MAP_PRIVATE | MAP_ANON;
    const int protection = PROT_READ | PROT_WRITE;
    ptr = mmap(NULL, num_bytes, protection, flags, -1 /*file descriptor*/, 0 /*offset*/);
    ptr = (ptr != MAP_FAILED) ? ptr : NULL;
  }
  else {
    ptr = malloc(num_bytes);
  }

  return ptr;
}

void page_aligned_deallocate( void * ptr, size_t num_bytes)
{
  if ( use_mmap(num_bytes) ) {
    munmap(ptr,num_bytes);
  }
  else {
    free(ptr);
  }
}


}} // namespace stk::detail

