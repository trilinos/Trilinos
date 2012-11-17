/**   ------------------------------------------------------------
 *    Copyright 2005-2009, 2011 Sandia Corporation.
 *    Under the terms of Contract DE-AC04-94AL85000, there is a
 *    non-exclusive license for use of this work by or on behalf
 *    of the U.S. Government.  Export of this program may require
 *    a license from the United States Government.
 *    ------------------------------------------------------------
 */

#include <pwd.h>
#include <unistd.h>

#include <iostream>
#include <ostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cstdlib>
#include <stdexcept>
#include <iomanip>
#include <algorithm>
#include <locale>

#include <stk_util/util/FeatureTest.hpp>
#include <stk_util/diag/Env.hpp>
#include <stk_util/diag/Platform.hpp>
#include <stk_util/parallel/Exception.hpp>
#include <stk_util/parallel/ExceptionReport.hpp>
#include <stk_util/parallel/MPI.hpp>
#ifdef STK_BUILT_IN_SIERRA
#  include <stk_util/parallel/mpih.hpp>
#endif
#include <stk_util/environment/ProductRegistry.hpp>

#include <stk_util/diag/Writer.hpp>
#include <stk_util/diag/SlibDiagWriter.hpp>
#include <stk_util/diag/Timer.hpp>
#include <stk_util/diag/Trace.hpp>

#include <fcntl.h>

#if defined(__GNUC__)
#include <fstream>
#ifndef __APPLE__
#include <malloc.h>
#else
#include <sys/malloc.h>
#endif
#include <cstdlib>
#include <sys/time.h>
#include <sys/resource.h>
#if __GNUC__ == 3 || __GNUC__ == 4
#include <cxxabi.h>
#endif

#elif defined(__PGI)
#include <fstream>
#include <malloc.h>
#include <cstdlib>
#include <sys/time.h>
#include <sys/resource.h>

#elif defined(__sun)
#include <fstream>
#include <procfs.h>
#include <sys/resource.h>
#include <sys/systeminfo.h>
#include <sys/utsname.h>
#include <sys/time.h>

#elif defined(__SUNPRO_CC)
#include <sys/resource.h>
#include <sys/time.h>
#include <sys/utsname.h>
#include <netdb.h>
#endif

#if defined(__PUMAGON__)
extern "C" {
#include <util.h>
#include <sys/param.h>
}

#elif defined(__sgi)
#include <sys/time.h>
#include <sys/resource.h>

#elif defined(REDS)
#include <sys/param.h>
#include <sys/utsname.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <catamount/catmalloc.h>
extern void *__heap_start;		///< Magic address of end of instruction/static data
extern void *__heap_end;		///< Magic address of end of heap

#elif defined(__JVN)
#include <sys/param.h>
#include <sys/utsname.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/wait.h>
#include <unistd.h>

#elif defined(__IBMC__) || defined(__IBMCPP__)
#include <sys/utsname.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <netdb.h>

#else
#include <sys/utsname.h>
#include <sys/time.h>
#include <netdb.h>
#endif

#include <stk_util/util/MallocUsed.h>

//
//  NKC 2/21/08 Some of the calls in this file bomb out purify.  Uncomment the line below if trying
//  to build a purify executable.
//
//#define PURIFY_BUILD

#if defined(REDS)
namespace {
  size_t get_redstorm_base_available_memory();
}
#endif

namespace sierra {

#ifdef SIERRA_USE_PLATFORM_DEMANGLER


#if defined(__GNUC__)

#if (__GNUC__ == 3)
std::string
demangle(
  const char *	symbol)
{
#ifdef PURIFY_BUILD
  return symbol;
#else
  std::string   s;
  int		status;

  // GR: PathScale versions 4.x.x attempt to write the length of the demangled
  //     string to the third argument of __cxa_demangle even if a NULL value is
  //     provided. This seems to be a compiler bug. The following code change
  //     provides a work around to prevent a segfault.
  //char *demangled_symbol = abi::__cxa_demangle(symbol, 0, 0, &status);
  size_t len = 0;
  char *demangled_symbol = abi::__cxa_demangle(symbol, 0, &len, &status);

  if (demangled_symbol) {
    s = std::string(demangled_symbol);
    free(demangled_symbol);
  }

  if (status != 0)
    s = std::string(symbol);

  return s;
#endif
}

#elif (__GNUC__ == 4)
std::string
demangle(
  const char *	symbol)
{
#ifdef PURIFY_BUILD
  return symbol;
#else
  std::string   s;

  int		status;

  // GR: PathScale versions 4.x.x attempt to write the length of the demangled
  //     string to the third argument of __cxa_demangle even if a NULL value is
  //     provided. This seems to be a compiler bug. The following code change
  //     provides a work around to prevent a segfault.
  //char *demangled_symbol = __cxxabiv1::__cxa_demangle(symbol, 0, 0, &status);
  size_t len = 0;
  char *demangled_symbol = __cxxabiv1::__cxa_demangle(symbol, 0, &len, &status);

  if (demangled_symbol) {
    s = std::string(demangled_symbol);
    free(demangled_symbol);
  }

  if (status != 0)
    s = std::string(symbol);

  return s;
#endif
}
#endif // (__GNUC__ == 3)

#else
std::string demangle(const char *symbol) {
  return symbol;
}
#endif // defined(__GNUC__)

#else
const char *demangle(const char *symbol) {
  return symbol;
}
#endif // SIERRA_USE_PLATFORM_DEMANGLER


namespace Env {

#if defined(REDS)
// Manipulate standard output buffering
void
startup_preparallel_platform()
{
  static char redstorm_cout_buf[32678];
  std::cout.rdbuf()->pubsetbuf(redstorm_cout_buf, sizeof(redstorm_cout_buf));

  static char redstorm_stdout_buf[32678];
  std::setvbuf(stdout, redstorm_stdout_buf, _IOFBF, sizeof(redstorm_stdout_buf));

  static char redstorm_cerr_buf[32678];
  std::cerr.rdbuf()->pubsetbuf(redstorm_cerr_buf, sizeof(redstorm_cerr_buf));

  static char redstorm_stderr_buf[32678];
  std::setvbuf(stderr, redstorm_stderr_buf, _IOFBF, sizeof(redstorm_stderr_buf));
  
  get_available_memory();
}

#elif defined(_AIX)
// Cleanup AIX locale initialization problems
void
startup_preparallel_platform()
{
  std::locale loc("POSIX");

  std::locale::global(loc); // std::locale::classic());
  std::cout.imbue(loc); // std::locale::classic());
  std::cin.imbue(loc); // std::locale::classic());

  std::ostringstream strout;
  strout << "Don't ask why the IBM locale works if I do this " << 10000000 << std::endl;
}

#else
void
startup_preparallel_platform()
{}
#endif


double
wall_now()
{
  timeval tp;
  struct timezone tz;
  gettimeofday(&tp, &tz);
  return (tp.tv_sec + (((double)(tp.tv_usec))/1000000.0));
}


double
cpu_now()
{
#if defined(REDS)
  struct rusage my_rusage;

  getrusage(RUSAGE_SELF, &my_rusage);

  return (double) (my_rusage.ru_utime.tv_sec)
    + ((double)(my_rusage.ru_utime.tv_usec))*1.0e-6;

#elif ! defined(__PGI)
  struct rusage my_rusage;

  getrusage(RUSAGE_SELF, &my_rusage);

  return (double) (my_rusage.ru_utime.tv_sec + my_rusage.ru_stime.tv_sec)
    + ((double)(my_rusage.ru_utime.tv_usec + my_rusage.ru_stime.tv_usec))*1.0e-6;
#else
  return 0;
#endif
}


void
get_heap_info(
  size_t &		heap_size,
  size_t &		largest_free)
{
  heap_size = 0;
  largest_free = 0;

#if defined(SIERRA_HEAP_INFO)

# if defined(SIERRA_PTMALLOC3_ALLOCATOR) || defined(SIERRA_PTMALLOC2_ALLOCATOR)
  heap_size = malloc_used();
  
# elif 0 // if defined(REDS) // Redstorm now links in gnu's malloc
  static size_t reds_fragments;
  static unsigned long reds_total_free;
  static unsigned long reds_heap_size;
  static unsigned long reds_largest_free;

  ::heap_info(&reds_fragments, &reds_total_free, &reds_largest_free, &reds_heap_size);

  heap_size = reds_heap_size;
  largest_free = reds_largest_free;

  slibout.m(Slib::LOG_MEMORY) <<"reds_fragments " << reds_fragments
			      << ", reds_total_free " << reds_total_free
			      << ", reds_largest_free " << reds_largest_free
			      << ", reds_heap_size " << reds_heap_size << Diag::dendl;

# elif ( defined(__linux__) || defined(REDS) ) && ! defined(__IBMCPP__)
  static struct mallinfo minfo;
  minfo = mallinfo();
  heap_size = (unsigned int) minfo.uordblks + (unsigned int) minfo.hblkhd;
  largest_free = (unsigned int) minfo.fordblks;

  slibout.m(Slib::LOG_MEMORY) << "size_t size " << sizeof(size_t)*8 << " bits"
                              << ", heap size " << heap_size
                              << ", arena " << (unsigned int) minfo.arena
			      << ", ordblks " << minfo.ordblks
			      << ", smblks " << minfo.smblks
			      << ", hblks " << minfo.hblks
			      << ", hblkhd " << (unsigned int) minfo.hblkhd
			      << ", usmblks " << minfo.usmblks
			      << ", fsmblks " << minfo.fsmblks
			      << ", uordblks " << (unsigned int) minfo.uordblks
			      << ", fordblks " << (unsigned int) minfo.fordblks
			      << ", keepcost " << minfo.keepcost << Diag::dendl;


# elif defined(__sun)
  pstatus_t proc_status;

  std::ifstream proc("/proc/self/status", std::ios_base::in|std::ios_base::binary);
  if (proc) {
    proc.read((char *)&proc_status, sizeof(proc_status));
    heap_size = proc_status.pr_brksize;
    slibout.m(Slib::LOG_MEMORY) <<"pr_brksize " << proc_status.pr_brksize
				<< ", pr_stksize " << proc_status.pr_stksize << Diag::dendl;
  }
# endif
#endif // defined(SIERRA_HEAP_INFO)
}


size_t get_available_memory()
{
#if !defined(REDS)
  // The value returned for _SC_AVPHYS_PAGES is the amount of memory
  // the application can use without hindering any other process
  // (given that no other process increases its memory usage).
#if !defined(__APPLE__) && !defined(__FreeBSD__)
  static size_t pagesize = getpagesize();
  size_t avail = sysconf(_SC_AVPHYS_PAGES);
  return avail * pagesize;
#else
  // _SC_AVPHYS_PAGES does not exist on FreeBSD/Apple
  return 0;
#endif
#else
  // On redstorm, we get an estimate of the available memory at the
  // time that the application starts up and then we get the currently
  // available memory by subtracting the used memory from that initial
  // value. 
  static size_t initial_memory_size = 0;
  if (initial_memory_size == 0) {
    initial_memory_size = get_redstorm_base_available_memory();
  }
  return initial_memory_size - get_heap_usage();
#endif
}

void
get_memory_info(
  size_t &		memory_usage,
  size_t &		faults)
{
  memory_usage = 0;
  faults = 0;

#if defined(SIERRA_MEMORY_INFO)
# if defined(REDS)
  memory_usage = (size_t) __heap_start + get_heap_usage();

# elif defined(__linux__)
  std::ifstream proc("/proc/self/stat", std::ios_base::in|std::ios_base::binary);
  if (proc) {

    std::string s;
    int i;
    for (i = 0; i < 11; ++i)
      proc >> s;

    proc >> faults;
    ++i;

    for (; i < 22; ++i)
      proc >> s;

    proc >> memory_usage;
    ++i;
  }
# elif defined(__sun)
  {

    psinfo_t proc_info;

    std::ifstream proc("/proc/self/psinfo", std::ios_base::in|std::ios_base::binary);
    if (proc) {
      proc.read((char *)&proc_info, sizeof(proc_info));
      memory_usage = proc_info.pr_size*1024;
    }
  }

  {
    prusage_t proc_usage;

    std::ifstream proc("/proc/self/usage", std::ios_base::in|std::ios_base::binary);
    if (proc) {
      proc.read((char *)&proc_usage, sizeof(proc_usage));
      faults = proc_usage.pr_majf;
    }
  }
# endif
#endif // defined(SIERRA_MEMORY_INFO)
}


double
vm_now()
{
  size_t	memory_usage;
  size_t	faults;

  get_memory_info(memory_usage, faults);

  return (double) memory_usage;
}


std::string
hostname()
{
  char buf[255];
  ::gethostname(buf, sizeof(buf));
  return std::string(buf);
}


std::string
domainname()
{
#if defined(__PUMAGON__) || defined(REDS)
  return std::string(".sandia.gov");

#elif defined(__sun)
  std::string domain(".");
  char buf[255];

  ::sysinfo(SI_SRPC_DOMAIN, buf, sizeof(buf));
  if (std::strlen(buf)) {
    domain += buf;
  }
  return domain;

#else //#elif defined(__linux) || defined(__sgi) || defined(_AIX)
  std::string domain(".");
  char buf[255];

  ::getdomainname(buf, sizeof(buf));
  if (::strlen(buf)) {
    domain += buf;
  }
  return domain;

#endif
}


std::string
username()
{
#if defined(REDS) || defined(__CRAYXT_COMPUTE_LINUX_TARGET)
  if (get_param("username").empty())
    return "unknown";
  else
    return get_param("username");
#else
  struct passwd *user_info = ::getpwuid(::geteuid());

  return (user_info ? user_info->pw_name : "unknown");
#endif
}


std::string
hardware()
{
#ifndef __sgi
  struct utsname	uts_name;

  uname(&uts_name);

  return uts_name.machine;
#else
  std::string s;
  return s;
#endif
}


std::string
osname()
{
#ifndef __sgi
  struct utsname	uts_name;

  uname(&uts_name);

  return uts_name.sysname;
#else
  std::string s;
  return s;
#endif
}


std::string
osversion()
{
#ifndef __sgi
  struct utsname	uts_name;

  uname(&uts_name);

  return uts_name.release;
#else
  std::string s;
  return s;
#endif
}


int
pid()
{
  return ::getpid();
}


int
pgrp()
{
#if defined(__PUMAGON__) || defined(REDS)
  return 0;
#else
  return ::getpgrp();
#endif
}


bool
path_access(
  const std::string &	name,
  int			mode)

{
  return !name.empty() && ::access(name.c_str(), mode) == 0;
}


bool
path_exists(
  const std::string &	name)
{
  return path_access(name, F_OK);
}


bool
path_read_access(
  const std::string &	name)
{
  return path_access(name, R_OK);
}


bool
path_write_access(
  const std::string &	name)
{
  return path_access(name, W_OK);
}


namespace {

struct flock *
file_lock(
  short	type,
  short	whence)
{
//  /* %TRACE[SPEC]% */ Tracespec trace__("sierra::Fmwk::<unnamed>::file_lock( short type, short whence)"); /* %TRACE% */
  static struct flock ret;
  ret.l_type = type;
  ret.l_start = 0;
  ret.l_whence = whence;
  ret.l_len = 0;
  ret.l_pid = 0; //getpid();
  return &ret;
}

} // namespace <unnamed>

bool
write_lock(
  int		fd)
{
  int i =::fcntl(fd, F_SETLK, file_lock(F_WRLCK, SEEK_SET));
//   if (i == -1)
//     fmwkout << "Write lock failed " << errno << dendl;

  return i != -1;
}


bool
release_lock(
  int		fd)
{
  int i =::fcntl(fd, F_SETLK, file_lock(F_UNLCK, SEEK_SET));
//   if (i == -1)
//     fmwkout << "Release lock failed " << errno << dendl;

  return i != -1;
}


bool
read_lock(
  int		fd)
{
  return ::fcntl(fd, F_SETLK, file_lock(F_RDLCK, SEEK_SET)) != -1;
}


bool
append_lock(
  int		fd)
{
  return ::fcntl(fd, F_SETLK, file_lock(F_WRLCK, SEEK_END)) != -1;
}

} // namespace Env
} // namespace sierra

#if defined(REDS) 

#if defined(__GNUC__)
namespace {
  size_t get_redstorm_base_available_memory()
  {
    return 0;
  }
}
#else
// Written by Mike Davis
#include <catamount/data.h>
namespace {

  void stk_ptr (unsigned long *sp) {
    asm ("movq	%rsp, (%rdi)");
  }
  
  size_t get_redstorm_base_available_memory()
  {
    char *p1;
    size_t stack_top;
    size_t avail_mem;
    size_t heap_base;
    size_t stack_base;
    size_t stack_size;
    size_t unmapped_top = 4 * 1048576;
    size_t os_foot =    140 * 1048576;
    size_t avail_heap;
    
    /**
     * Get the current stack pointer,
     * and use that as an estimate of the stack top
     */
    stk_ptr (&stack_top);

    /**
     * Round stack_top up to the next multiple of 2 MB
     */
    stack_top &= ~0x1fffff;
    stack_top +=  0x200000;
    
    /**
     * Compute the available memory on the node,
     * as the stack top plus the size of an unmapped region
     * above the stack.
     */
    avail_mem = stack_top + unmapped_top;

    /**
     * Deduct the size of the OS footprint from available memory
     */
    avail_mem -= os_foot;

    /**
     * Divide the available memory among
     * the number of processes running on the node
     */
    avail_mem /= _my_pcb->upcb_vnm_degree;

    /**
     * Determine the address of the base of the heap,
     * estimated to be the address of an allocated test block
     */
    p1 = (char *) malloc (1);
    heap_base = (size_t) p1;
    free (p1);

    /**
     * Compute the available heap space as the difference between
     * the available memory for the process
     * and the base address of the heap
     */
    avail_heap = avail_mem - heap_base;

    /**
     * Determine the base address of the stack;
     * the address of the PCB is in a 2 MB page below the stack
     */
    stack_base = (size_t) _my_pcb;
    stack_base &= ~0x1fffff;
    stack_base +=  0x400000;

    /**
     * Deduct the size of the stack
     * and the size of the unmapped region at the top of the stack
     * from the available heap
     */
    stack_size = stack_top - stack_base;
    avail_heap -= unmapped_top;
    avail_heap -= stack_size;

    return avail_heap;
  }
}
#endif
#endif
