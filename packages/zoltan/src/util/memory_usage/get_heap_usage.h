// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

//////////////////////////////////////////////////////////////////////
// Subroutine to measure heap usage on several different platforms. //
//////////////////////////////////////////////////////////////////////

#include <iostream>
#if defined(__APPLE__)
#include <malloc/malloc.h>
#else
#include <malloc.h>
#endif

static size_t get_heap_usage()
{
  size_t heap_size = 0;

#if defined(__GNUC__) && defined(__linux__)
  static struct mallinfo minfo;
  minfo = mallinfo();
  heap_size = minfo.uordblks + minfo.hblkhd;

#elif defined(__APPLE__)
  malloc_statistics_t t = {0,0,0,0};
  malloc_zone_statistics(NULL, &t);
  heap_size = t.size_in_use;  

#elif defined(__sun)
  pstatus_t proc_status;

  std::ifstream proc("/proc/self/status", std::ios_base::in|std::ios_base::binary);
  if (proc) {
    proc.read((char *)&proc_status, sizeof(proc_status));
    heap_size = proc_status.pr_brksize;
  }
#endif
  return heap_size;
}
