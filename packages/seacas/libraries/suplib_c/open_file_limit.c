// Copyright(C) 1999-2021, 2023 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#if defined(WIN32) || defined(__WIN32__) || defined(_WIN32) || defined(_MSC_VER) ||                \
    defined(__MINGW32__) || defined(_WIN64) || defined(__MINGW64__)
#include <stdio.h>
#else
#include <limits.h>
#include <unistd.h>
#endif

int open_file_limit(void)
{
  // Returns maximum number of files that one process can have open
  // at one time. (POSIX)
#if defined(WIN32) || defined(__WIN32__) || defined(_WIN32) || defined(_MSC_VER) ||                \
    defined(__MINGW32__) || defined(_WIN64) || defined(__MINGW64__)
  int fdmax = _getmaxstdio();
#else
  int fdmax = sysconf(_SC_OPEN_MAX);
  if (fdmax == -1) {
    // POSIX indication that there is no limit on open files...
    fdmax = INT_MAX;
  }
#endif

  // File descriptors are assigned in order (0,1,2,3,...) on a per-process
  // basis.

  // Assume that we have stdin, stdout, stderr.
  // (3 total).

  return fdmax - 3;

  // Could iterate from 0..fdmax and check for the first EBADF (bad
  // file descriptor) error return from fcntl, but that takes too long
  // and may cause other problems.
  //
  // Another possibility is to do an open and see which descriptor is
  // returned -- take that as 1 more than the current count of open
  // files instead of assuming 3.
  //
}
