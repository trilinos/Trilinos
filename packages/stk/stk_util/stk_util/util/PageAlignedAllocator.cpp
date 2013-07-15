
#include <sys/mman.h>
#include <unistd.h>
#include <cstdlib>

#include <stk_util/util/PageAlignedAllocator.hpp>

namespace stk { namespace detail {

const size_t page_aligned_allocator_impl::half_page_size = (sysconf( _SC_PAGE_SIZE ) >> 1);
const int page_aligned_allocator_impl::m_mmap_flags = MAP_PRIVATE | MAP_ANON;
const int page_aligned_allocator_impl::m_mmap_protection = PROT_READ | PROT_WRITE;


void * page_aligned_allocator_impl::allocate(size_t num_bytes)
{
  void * ptr = mmap(NULL, num_bytes, m_mmap_protection, m_mmap_flags, -1 /*file descriptor*/, 0 /*offset*/);
  ptr = (ptr != MAP_FAILED) ? ptr : NULL;

  return ptr;
}

void page_aligned_allocator_impl::deallocate( void * ptr, size_t num_bytes)
{
  munmap(ptr,num_bytes);
}


}} // namespace stk::detail

