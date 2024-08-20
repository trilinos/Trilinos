// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

////////////////////////////////////////////////////////////////////////////////
//
// Much of what appears below was taken from getRSS.c from
//   Author:   David Robert Nadeau
//   Site:     http://NadeauSoftware.com/
//   License:  Creative Commons Attribution 3.0 Unported License
//             http://creativecommons.org/licenses/by/3.0/deed.en_US
// and was extended to gather memory usage information from all processors in
// Teuchos::Comm<int>& comm.
//
// Note:  On Windows, link with psapi.lib.
//
////////////////////////////////////////////////////////////////////////////////

#include <Panzer_MemUtils.hpp>
#include <Teuchos_CommHelpers.hpp>

#if       defined(_WIN32)
#  include <windows.h>
#  include <psapi.h>
#elif     defined(__unix__)      || \
          defined(__unix)        || \
          defined(unix)          || \
          (defined(__APPLE__) &&    \
           defined(__MACH__))
#  include <unistd.h>
#  include <sys/resource.h>
#  if       defined(__APPLE__) && \
            defined(__MACH__)
#    include <mach/mach.h>
#  elif     (defined(_AIX)             ||    \
             defined(__TOS__AIX__))       || \
            (defined(__sun__)          ||    \
             defined(__sun)            ||    \
             defined(sun)              &&    \
             (defined(__SVR4)       ||       \
              defined(__svr4__)))
#    include <fcntl.h>
#    include <procfs.h>
#  elif     defined(__linux__)     || \
            defined(__linux)       || \
            defined(linux)         || \
            defined(__gnu_linux__)
#    include <stdio.h>
#  endif // defined(__APPLE__) && ...
#else
#  error "Cannot define getPeakRSS( ) or getCurrentRSS( ) for an unknown OS."
#endif // defined(_WIN32)

namespace panzer
{
  //////////////////////////////////////////////////////////////////////////////
  //
  //  printMemoryUsage()
  //
  //////////////////////////////////////////////////////////////////////////////
  void printMemoryUsage(std::ostream& s, const Teuchos::Comm<int>& comm)
  {
    MemUsage mem = getMemoryUsage(comm);
    printMemoryUsage(s, comm, mem);
    return;
  } // end of printMemoryUsage()

  //////////////////////////////////////////////////////////////////////////////
  //
  //  printMemoryUsage()
  //
  //////////////////////////////////////////////////////////////////////////////
  void printMemoryUsage(std::ostream& s, const Teuchos::Comm<int>& comm,
    const MemUsage& mem)
  {
    using std::endl;
    if (0 == comm.getRank())
    {
      s << "Estimated memory usage across all processors:" << endl
        << "        Current       Peak        "         << endl
        << "        ------------  ------------"         << endl << "  Min:  ";
      pretty(s, mem.currMin); pretty(s, mem.peakMin); s << endl << "  Max:  ";
      pretty(s, mem.currMax); pretty(s, mem.peakMax); s << endl << "  Tot:  ";
      pretty(s, mem.currTot); pretty(s, mem.peakTot); s << endl;
    }
    return;
  } // end of printMemoryUsage()

  //////////////////////////////////////////////////////////////////////////////
  //
  //  pretty()
  //
  //////////////////////////////////////////////////////////////////////////////
  void pretty(std::ostream& s, size_t num)
  {
    s << std::fixed;
    s.precision(3);
    s << std::setw(9) << std::setfill(' ');
    if (num < (1 << 10))
      s << num << "  B  ";
    else if (num < (1 << 20))
      s << (static_cast<double>(num) / (1 << 10)) << " KB  ";
    else if (num < (1 << 30))
      s << (static_cast<double>(num) / (1 << 20)) << " MB  ";
    else
      s << (static_cast<double>(num) / (1 << 30)) << " GB  ";
    return;
  } // end of pretty()

  //////////////////////////////////////////////////////////////////////////////
  //
  //  getMemoryUsage()
  //
  //////////////////////////////////////////////////////////////////////////////
  MemUsage getMemoryUsage(const Teuchos::Comm<int>& comm)
  {
    MemUsage current = getCurrentRSS(comm);
    MemUsage peak    = getPeakRSS(comm);
    return current + peak;
  } // end of getMemoryUsage()

  //////////////////////////////////////////////////////////////////////////////
  //
  //  getPeakRSS()
  //
  //////////////////////////////////////////////////////////////////////////////
  MemUsage getPeakRSS(const Teuchos::Comm<int>& comm)
  {
    size_t mem(0);

    // Windows
#   if       defined(_WIN32)
      PROCESS_MEMORY_COUNTERS info;
      GetProcessMemoryInfo(GetCurrentProcess(), &info, sizeof(info));
      mem = (size_t)(info.PeakWorkingSetSize);

    // AIX and Solaris
#   elif     (defined(_AIX)             ||    \
              defined(__TOS__AIX__))       || \
             (defined(__sun__)          ||    \
              defined(__sun)            ||    \
              defined(sun)              &&    \
              (defined(__SVR4)       ||      \
               defined(__svr4__)))
      struct psinfo psinfo;
      int fd = -1;
      if ((fd = open("/proc/self/psinfo", O_RDONLY)) == -1)
        mem = (size_t)(0L);  // Can't open?
      if (read(fd, &psinfo, sizeof(psinfo)) != sizeof(psinfo))
      {
        close(fd);
        mem = (size_t)(0L);  // Can't read?
      }
      close(fd);
      mem = (size_t)(psinfo.pr_rssize * 1024L);

    // BSD, Linux, and OSX
#   elif     defined(__unix__)      || \
             defined(__unix)        || \
             defined(unix)          || \
             (defined(__APPLE__) &&    \
              defined(__MACH__))
      struct rusage rusage;
      getrusage(RUSAGE_SELF, &rusage);
#     if       defined(__APPLE__) && \
               defined(__MACH__)
        mem = (size_t)(rusage.ru_maxrss);
#     else
        mem = (size_t)(rusage.ru_maxrss * 1024L);
#     endif // defined(__APPLE__) && ...

    // Unknown OS
#   else
      mem = (size_t)(0L);  // Unsupported.
#   endif // defined(_WIN32)
    return reduceMemUsage(mem, comm, MEM_USAGE_PEAK);
  } // end of getPeakRSS()

  //////////////////////////////////////////////////////////////////////////////
  //
  //  getCurrentRSS()
  //
  //////////////////////////////////////////////////////////////////////////////
  MemUsage getCurrentRSS(const Teuchos::Comm<int>& comm)
  {
    size_t mem(0);

    // Windows
#   if       defined(_WIN32)
      PROCESS_MEMORY_COUNTERS info;
      GetProcessMemoryInfo(GetCurrentProcess(), &info, sizeof(info));
      mem = (size_t)(info.WorkingSetSize);

    // OSX
#   elif     defined(__APPLE__) && \
             defined(__MACH__)
      struct mach_task_basic_info info;
      mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;
      if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO,
        (task_info_t)(&info), &infoCount) != KERN_SUCCESS)
        mem = (size_t)(0L);  // Can't access?
      mem = (size_t)(info.resident_size);

    // Linux
#   elif     defined(__linux__)     || \
             defined(__linux)       || \
             defined(linux)         || \
             defined(__gnu_linux__)
      long rss = 0L;
      FILE* fp = NULL;
      if ((fp = fopen("/proc/self/statm", "r")) == NULL)
        mem = (size_t)(0L);  // Can't open?
      if (fscanf(fp, "%*s%ld", &rss) != 1)
      {
        fclose(fp);
        mem = (size_t)(0L);  // Can't read?
      }
      fclose(fp);
      mem = (size_t)(rss) * (size_t)(sysconf(_SC_PAGESIZE));

    // AIX, BSD, Solaris, and Unknown OS
#   else
      mem = (size_t)(0L);  // Unsupported.
#   endif // defined(_WIN32)
    return reduceMemUsage(mem, comm, MEM_USAGE_CURRENT);
  } // end of getCurrentRSS()

  //////////////////////////////////////////////////////////////////////////////
  //
  //  reduceMemUsage()
  //
  //////////////////////////////////////////////////////////////////////////////
  MemUsage reduceMemUsage(size_t& mem, const Teuchos::Comm<int>& comm,
    const MemUsageType& type)
  {
	  size_t min(0), max(0), tot(0);
    Teuchos::reduceAll(comm, Teuchos::REDUCE_MIN, 1, &mem, &min);
    Teuchos::reduceAll(comm, Teuchos::REDUCE_MAX, 1, &mem, &max);
    Teuchos::reduceAll(comm, Teuchos::REDUCE_SUM, 1, &mem, &tot);
    MemUsage result;
    switch (type)
    {
      case MEM_USAGE_CURRENT:
        result.currMin = min;
        result.currMax = max;
        result.currTot = tot;
        break;
      case MEM_USAGE_PEAK:
        result.peakMin = min;
        result.peakMax = max;
        result.peakTot = tot;
        break;
    }
    return result;
  } // end of reduceMemUsage()
} // end namespace panzer
