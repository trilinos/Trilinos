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
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <iostream>
#include <sstream>
#include <fstream>

#include <MueLu_config.hpp>

#define HAVE_MUELU_PROC_SELF_STATUS

#ifdef HAVE_MUELU_GOOGLE_PERFTOOLS
//#include <google/heap-profiler.h>
#include </home/jngaida/google-perftools-1.8.3/src/google/heap-profiler.h> // TODO :-)
#endif

#include <MueLu_MemoryProfiler.hpp>

std::string GetMemoryUsage() {

  std::ostringstream mem;

#ifdef HAVE_MUELU_PROC_SELF_STATUS
  //TODO: test if /proc/self/status exist on the system instead of #ifdef

  std::ifstream proc("/proc/self/status");
  std::string s;
  while(getline(proc, s), !proc.fail()) {
    if(s.substr(0, 6) == "VmSize") {
      mem << s;
      return mem.str();
    }
  }

#endif // HAVE_MUELU_PROC_SELF_STATUS

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
    iss >> sub; iss >> sub; iss >> sub; // skip 3 first substring
    iss >> sub;
    double MB = atof(sub.c_str()) / (1024*1024);

    // print
    if (description != "") {
      std::ostringstream sname; sname.precision(1);
      sname << description << ": " << std::fixed << MB << " MB";
      std::cout << sname.str() << std::endl;
    }

    // dump to file
    if (filename != "") {
      std::ofstream out(filename.c_str(), std::ios::out | std::ios::binary);
      if(!out) { std::cout << "Cannot open output file: " << filename << std::endl; return; }
      out.write(profile, strlen(profile));
      out.close();

      // dump to file using HeapProfilerDump:
      // HeapProfilerDump(filename.c_str());
    }

    free(profile);

  }
#endif

}
