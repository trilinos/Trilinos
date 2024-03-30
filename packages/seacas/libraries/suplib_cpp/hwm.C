#include <hwm.h>
// For memory utilities...
#if defined(WIN32) || defined(__WIN32__) || defined(_WIN32) || defined(_MSC_VER) ||                \
    defined(__MINGW32__) || defined(_WIN64) || defined(__MINGW64__)
#define __SUP_WINDOWS__ 1
#endif

#if defined(__SUP_WINDOWS__)
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#if 0
#include <psapi.h>
#endif
#include <windows.h>

#elif defined(__unix__) || defined(__unix) || defined(unix) ||                                     \
    (defined(__APPLE__) && defined(__MACH__))
#include <sys/resource.h>
#include <unistd.h>

#if defined(__APPLE__) && defined(__MACH__)
#include <mach/mach.h>

#elif (defined(_AIX) || defined(__TOS__AIX__)) ||                                                  \
    (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
#include <fcntl.h>
#include <procfs.h>

#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
#include <stdio.h>
#endif
#endif

#if defined(BGQ_LWK) && defined(__linux__)
#include <spi/include/kernel/location.h>
#include <spi/include/kernel/memory.h>
#endif

size_t get_hwm_memory_info()
{
  // Code from http://nadeausoftware.com/sites/NadeauSoftware.com/files/getRSS.c
  size_t memory_usage = 0;
#if defined(__SUP_WINDOWS__)
#if 0
  /* Windows -------------------------------------------------- */
  PROCESS_MEMORY_COUNTERS info;
  GetProcessMemoryInfo(GetCurrentProcess(), &info, sizeof(info));
  memory_usage = (size_t)info.PeakWorkingSetSize;
#else
  memory_usage = 0;
#endif

#elif (defined(_AIX) || defined(__TOS__AIX__)) ||                                                  \
    (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
  /* AIX and Solaris ------------------------------------------ */
  struct psinfo psinfo;
  int           fd = -1;
  if ((fd = open("/proc/self/psinfo", O_RDONLY)) == -1)
    return (size_t)0L; /* Can't open? */
  if (read(fd, &psinfo, sizeof(psinfo)) != sizeof(psinfo)) {
    close(fd);
    return (size_t)0L; /* Can't read? */
  }
  close(fd);
  memory_usage = (size_t)(psinfo.pr_rssize * 1024L);

#elif (defined(__APPLE__) && defined(__MACH__)) || (defined(__linux__) && !defined(BGQ_LWK))
  /* BSD, Linux, and OSX -------------------------------------- */
  struct rusage rusage;
  getrusage(RUSAGE_SELF, &rusage);
#if defined(__APPLE__) && defined(__MACH__)
  memory_usage = (size_t)rusage.ru_maxrss;
#else
  memory_usage = (size_t)(rusage.ru_maxrss * 1024L);
#endif
#endif
  return memory_usage;
}
