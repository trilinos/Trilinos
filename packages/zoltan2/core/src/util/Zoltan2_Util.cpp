// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file Zoltan2_Util.cpp
 *  \brief Useful namespace methods.
 */

#include <Zoltan2_Util.hpp>
#include <Teuchos_OpaqueWrapper.hpp>

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

#ifndef _MSC_VER
#include <unistd.h>
#endif

namespace Zoltan2{

#ifndef _WIN32
/* On a linux node, find the total memory currently allocated
 * to this process.
 * Return the number of kilobytes allocated to this process.
 * Return 0 if it is not possible to determine this.
 */
long getProcessKilobytes()
{
long pageSize;

#ifdef _SC_PAGESIZE
  pageSize = sysconf(_SC_PAGESIZE);
#else
#warning "Page size query is not possible.  No per-process memory stats."
  return 0;
#endif

  pid_t pid = getpid();
  std::ostringstream fname;
  fname << "/proc/" << pid << "/statm";
  std::ifstream memFile;

  try{
    memFile.open(fname.str().c_str());
  }
  catch (...){
    return 0;
  }

  char buf[128];
  memset(buf, 0, 128);
  while (memFile.good()){
    memFile.getline(buf, 128);
    break;
  }

  memFile.close();

  std::istringstream sbuf(buf);
  long totalPages;
  sbuf >> totalPages;

  long pageKBytes = pageSize / 1024;
  totalPages = atol(buf);

  return totalPages * pageKBytes;
}
#else
long getProcessKilobytes()
{
#pragma message ("Zoltan2_Util.cpp: Page size query is not implemented on windows.  No per-process memory stats.")
  return 0;
}
#endif

} // namespace Zoltan2
