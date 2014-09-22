// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Governement retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#include <stdio.h>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <unistd.h>
#include <iomanip>
#include <vector>
#include <stk_util/environment/memory_util.hpp>

#if defined(__APPLE__)
#include<mach/task.h>
#include<mach/mach_init.h>
#endif
#ifdef __linux__
#  ifdef BGQ_LWK
#    include <spi/include/kernel/memory.h>
#    include <spi/include/kernel/location.h>
#  else
#    define PROCFS
#  endif
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

#elif defined(BGQ_LWK)
  uint64_t heap;
  Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAP, &heap);
  memory = heap;

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

  int  proc=0, cores=0;

  while (1)
  {
    if(!std::getline(proc_cpuinfo, line)) 
    {
      proc_cpuinfo.close();

      if(procinfo.size() == 0u)
	procinfo.push_back(1);  

      return;
    }

    cores = 1;

    if (line.substr(0, 9) == "processor")
    {
      size_t found = line.find_last_of(":");
      proc_string = line.substr(found+1);
      std::istringstream iss(proc_string);
      iss >> proc;
      have_processor = true;      
    }
    else if (line.substr(0, 9) == "cpu cores")
    {
      size_t found = line.find_last_of(":");
      core_string = line.substr(found+1);
      std::istringstream iss(core_string);
      iss >> cores;
      have_cores = true;      
    }

    if(have_processor)
    {
      if(static_cast<size_t>(proc) == procinfo.size())
	procinfo.push_back(cores);

      have_processor = false;            
    }

    if(have_cores)
    {
      int ind = procinfo.size()-1;

      if(have_processor)
	ind = proc;

      procinfo[ind] = cores;

      have_cores     = false;
    }
  }
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
    if(!std::getline(proc_meminfo, line)) 
    {
      proc_meminfo.close();
      return;
    }
    
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
} // namespace stk
