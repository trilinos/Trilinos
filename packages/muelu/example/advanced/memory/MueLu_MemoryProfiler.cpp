// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <iostream>
#include <sstream>
#include <fstream>

#include <MueLu_config.hpp>

#define HAVE_MUELU_PROC_SELF_STATUS

#ifdef HAVE_MUELU_GOOGLE_PERFTOOLS
//#include <google/heap-profiler.h>
#include </home/jngaida/google-perftools-1.8.3/src/google/heap-profiler.h>  // TODO :-)
#endif

#include <MueLu_MemoryProfiler.hpp>

std::string GetMemoryUsage() {
  std::ostringstream mem;

#ifdef HAVE_MUELU_PROC_SELF_STATUS
  // TODO: test if /proc/self/status exist on the system instead of #ifdef

  std::ifstream proc("/proc/self/status");
  std::string s;
  while (getline(proc, s), !proc.fail()) {
    if (s.substr(0, 6) == "VmSize") {
      mem << s;
      return mem.str();
    }
  }

#endif  // HAVE_MUELU_PROC_SELF_STATUS

  return mem.str();
}

void MemoryUsageStart(const std::string& autoLogPrefix) {
#ifdef HAVE_MUELU_GOOGLE_PERFTOOLS
  HeapProfilerStart("auto-profiling");
#endif
}

void MemoryUsageStop() {
#ifdef HAVE_MUELU_GOOGLE_PERFTOOLS
  HeapProfilerStop();
#endif
}

void PrintMemoryUsage(const std::string& description, const std::string& filename) {
#ifdef HAVE_MUELU_PROC_SELF_STATUS
  std::cout << description << ": " << GetMemoryUsage() << std::endl;
#endif

#ifdef HAVE_MUELU_GOOGLE_PERFTOOLS
  if (IsHeapProfilerRunning()) {
    char* profile = GetHeapProfile();

    std::istringstream iss(profile);
    std::string sub;
    iss >> sub;
    iss >> sub;
    iss >> sub;  // skip 3 first substring
    iss >> sub;
    double MB = atof(sub.c_str()) / (1024 * 1024);

    // print
    if (description != "") {
      std::ostringstream sname;
      sname.precision(1);
      sname << description << ": " << std::fixed << MB << " MB";
      std::cout << sname.str() << std::endl;
    }

    // dump to file
    if (filename != "") {
      std::ofstream out(filename.c_str(), std::ios::out | std::ios::binary);
      if (!out) {
        std::cout << "Cannot open output file: " << filename << std::endl;
        return;
      }
      out.write(profile, strlen(profile));
      out.close();

      // dump to file using HeapProfilerDump:
      // HeapProfilerDump(filename.c_str());
    }

    free(profile);
  }
#endif
}
