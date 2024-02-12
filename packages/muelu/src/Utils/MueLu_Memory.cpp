// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <sstream>
#include <fstream>
#include "MueLu_Memory.hpp"

#include <iostream>  // TODO: remove
#include <unistd.h>
#include <time.h>
#ifdef MUELU_USE_MALLINFO
#include <malloc.h>
#endif

//#define MUELU_USE_MALLINFO

namespace MueLu {

namespace MemUtils {

std::string PrintMemoryUsage() {
#ifdef MUELU_USE_MALLINFO
  struct mallinfo mem_stats = mallinfo();
  double memory             = mem_stats.hblkhd + mem_stats.usmblks + mem_stats.uordblks;

  char memchar[128];
  sprintf(memchar, "%12.1f MB", memory / 1048576.0);
  std::string mem(memchar);

  return mem;
#else
  std::ostringstream mem;
  std::ifstream proc("/proc/self/status");
  std::string s;

  mem << PrintMemoryInfo() << " ";
  while (getline(proc, s), !proc.fail()) {
    if (s.substr(0, 6) == "VmSize") {
      mem << s;
      return mem.str();
    }
  }
  return mem.str();
#endif
}

std::string PrintMemoryInfo() {
#ifdef MUELU_USE_MALLINFO
  struct mallinfo mem_stats = mallinfo();
  double memory             = mem_stats.hblkhd + mem_stats.usmblks + mem_stats.uordblks;

  char memchar[128];
  sprintf(memchar, "%12.1f MB", memory / 1048576.0);
  std::string mem(memchar);

  return mem;
#else
  std::ostringstream mem;
  std::ifstream proc("/proc/meminfo");
  std::string s;
  while (getline(proc, s), !proc.fail()) {
    if (s.substr(0, 7) == "MemFree") {
      mem << s;
      return mem.str();
    }
  }
  return mem.str();
#endif
}

void ReportTimeAndMemory(Teuchos::Time const &timer, Teuchos::Comm<int> const &Comm) {
  double maxTime = 0, minTime = 0, avgTime = 0;
  double localTime = timer.totalElapsedTime();
#ifdef HAVE_MPI
  int ntimers = 1, root = 0;
  MPI_Reduce(&localTime, &maxTime, ntimers, MPI_DOUBLE, MPI_MAX, root, MPI_COMM_WORLD);
  MPI_Reduce(&localTime, &minTime, ntimers, MPI_DOUBLE, MPI_MIN, root, MPI_COMM_WORLD);
  MPI_Reduce(&localTime, &avgTime, ntimers, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
#else
  maxTime = localTime;
  minTime = localTime;
  avgTime = localTime;
#endif
  avgTime /= Comm.getSize();
  // std::cout << "(" << Comm.getRank() << ") " << localTime << std::endl;
  if (Comm.getRank() == 0) {
    std::cout << "&&&" << timer.name()
              << " max=" << maxTime << " min=" << minTime << " avg=" << avgTime << std::endl;
    std::cout << "&&&" << timer.name() << " " << MemUtils::PrintMemoryUsage() << std::endl;
  }
}  // ReportTimeAndMemory

}  // namespace MemUtils

}  // namespace MueLu
