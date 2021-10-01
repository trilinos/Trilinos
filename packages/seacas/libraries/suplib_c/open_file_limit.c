// Copyright(C) 1999-2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include <limits.h>
#include <unistd.h>

int open_file_limit()
{
  // Returns maximum number of files that one process can have open
  // at one time. (POSIX)
#if !defined(_WIN64) && !defined(WIN32) && !defined(_WINDOWS) && !defined(_MSC_VER)
  int fdmax = sysconf(_SC_OPEN_MAX);
  if (fdmax == -1) {
    // POSIX indication that there is no limit on open files...
    fdmax = INT_MAX;
  }
#else
  int fdmax = _getmaxstdio();
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
