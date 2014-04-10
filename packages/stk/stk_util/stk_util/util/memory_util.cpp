#include <stdio.h>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <unistd.h>
#include <iomanip>
#include <vector>
#include <stk_util/util/memory_util.hpp>

#if defined(__APPLE__)
#include<mach/task.h>
#include<mach/mach_init.h>
#endif
#ifdef __linux__
#define PROCFS
#endif

namespace stk
{

// return current memory usage in bytes
size_t get_memory_usage_now()
{
  size_t memory = 0;

#if defined(__APPLE__)
  struct task_basic_info t_info;
  mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;

  if (KERN_SUCCESS != task_info(mach_task_self(),
      TASK_BASIC_INFO,
      reinterpret_cast<task_info_t>(&t_info),
      &t_info_count))
  {
    return 0;
  }
  memory = t_info.resident_size;

#elif defined(__LIBCATAMOUNT__) && defined(REDSTORM_HEAP_INFO)
  size_t frags;
  unsigned long total_free, largest_free, total_used;
  heap_info( &frags, &total_free, &largest_free, &total_used );
  memory = total_used * 1024;

#elif defined(PROCFS)
  unsigned long rss_pages = 0;
  // Read memory size data from /proc/pid/stat
  // see "man proc" for details.
  std::ifstream proc_stat("/proc/self/stat");
  if (proc_stat)
  {
    std::string buf;
    proc_stat
       >> buf >> buf >> buf >> buf >> buf
       >> buf >> buf >> buf >> buf >> buf
       >> buf >> buf >> buf >> buf >> buf
       >> buf >> buf >> buf >> buf >> buf
       >> buf >> buf >> buf
       >> rss_pages;
    proc_stat.close();
  }
  memory = rss_pages * sysconf( _SC_PAGESIZE);
#endif

  /* Success */
  return memory;
}

// return current resident set size in bytes
void get_memory_usage(size_t & now, size_t & hwm)
{
#if defined (PROCFS)
  std::string vmhwm;
  std::string vmrss;
  std::string line(128,'\0');

  now = 0;
  hwm = 0;

  /* Read memory size data from /proc/pid/status
   * run "man proc" to get info on the contents of /proc/self/status
   */
  std::ifstream proc_status("/proc/self/status");
  if (!proc_status) return;

  while (vmhwm.empty() || vmrss.empty())
  {

    if(!std::getline(proc_status, line)) return;

    /* Find VmHWM */
    else if (line.substr(0, 6) == "VmHWM:")
    {
      vmhwm = line.substr(7);
      std::istringstream iss(vmhwm);
      iss >> hwm;
      hwm *= 1024;
    }
    /* Find VmRSS */
    else if (line.substr(0, 6) == "VmRSS:")
    {
      vmrss = line.substr(7);
      std::istringstream iss(vmrss);
      iss >> now;
      now *= 1024;
    }
  }
  proc_status.close();

  #else
  now = get_memory_usage_now();
  static size_t s_hwm = 0;
  s_hwm = std::max(now, s_hwm);
  hwm = s_hwm;
#endif
}

void get_processor_count(std::vector<int> &procinfo)
{
  procinfo.clear();

#if defined (PROCFS)
  std::string proc_string;
  std::string core_string;
  std::string line(128,'\0');

  /* Read available cpu data from /proc/cpuinfo
   */
  std::ifstream proc_cpuinfo("/proc/cpuinfo");
  if (!proc_cpuinfo)
  {
    procinfo.push_back(1);
    return;
  }

  bool have_processor = false;
  bool have_cores     = false;

  int  proc, cores;

  while (1)
  {
    if(!std::getline(proc_cpuinfo, line))
      return;
    else if (line.substr(0, 11) == "processor	:")
    {
      proc_string = line.substr(12);
      std::istringstream iss(proc_string);
      iss >> proc;
      have_processor = true;
    }
    else if (line.substr(0, 11) == "cpu cores	:")
    {
      core_string = line.substr(12);
      std::istringstream iss(core_string);
      iss >> cores;
      have_cores = true;
    }

    if(have_cores && have_processor)
    {
      if((size_t) proc == procinfo.size())
	procinfo.push_back(cores);

      have_cores     = false;
      have_processor = false;
    }
  }
  proc_cpuinfo.close();
#else
  procinfo.push_back(1);
#endif
}

// return memory available
void get_memory_available(size_t & avail)
{
  avail = 0;

#if defined (PROCFS)
  std::string memtotal;
  std::string line(128,'\0');

  /* Read available memory size data from /proc/meminfo
   */
  std::ifstream proc_meminfo("/proc/meminfo");
  if (!proc_meminfo) return;

  while (memtotal.empty())
  {
    if(!std::getline(proc_meminfo, line)) return;

    /* Find MemTotal */
    else if (line.substr(0, 9) == "MemTotal:")
    {
      memtotal = line.substr(10);
      std::istringstream iss(memtotal);
      iss >> avail;
      avail *= 1024;
    }
  }
  proc_meminfo.close();
#endif
}


void get_memory_high_water_mark_across_processors(MPI_Comm comm, size_t& hwm_max, size_t& hwm_min, size_t& hwm_avg)
{
  size_t now = 0;
  size_t hwm = 0;
  stk::get_memory_usage(now, hwm);

  get_max_min_avg(comm, hwm, hwm_max, hwm_min, hwm_avg);
}

void get_memory_available_across_processors(MPI_Comm comm, size_t& avail_max, size_t& avail_min, size_t& avail_avg)
{
  size_t avail = 0;
  std::vector<int> coreinfo;

  stk::get_memory_available(avail);
  stk::get_processor_count(coreinfo);

  avail /= coreinfo.size();

  get_max_min_avg(comm, avail, avail_max, avail_min, avail_avg);
}

std::string human_bytes(size_t arg_bytes)
{
  double bytes = arg_bytes;
  const double K = 1024;
  const double M = K*1024;
  const double G = M*1024;

  std::ostringstream out;
  if (bytes < K) {
    out << std::setprecision(4) << bytes << " B";
  } else if (bytes < M) {
    bytes /= K;
    out << std::setprecision(4) << bytes << " K";
  } else if (bytes < G) {
    bytes /= M;
    out << std::setprecision(4) << bytes << " M";
  } else {
    bytes /= G;
    out << std::setprecision(4) << bytes << " G";
  }
  return out.str();
}

} // namespace stk
