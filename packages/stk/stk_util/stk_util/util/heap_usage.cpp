// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

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
