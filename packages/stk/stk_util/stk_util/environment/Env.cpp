/*
// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
 */

#include <stk_util/stk_config.h>
#include <stk_util/environment/Env.hpp>
#include <time.h>                       // for localtime, strftime, time_t
#include <limits.h>                     // for PATH_MAX
#include <stddef.h>                     // for size_t
#include <unistd.h>                     // for getcwd, sleep
#include <cstdlib>                      // for exit, EXIT_FAILURE, NULL
#include <cstring>                      // for strlen, strcpy
#include <iomanip>                      // for operator<<, setw
#include <iostream>                     // for cerr, cout
#include <map>                          // for map<>::mapped_type
#include <sstream>                      // for basic_ostream, operator<<, etc
#include <stk_util/util/Signal.hpp>     // for HUP_received
#include <stk_util/environment/EnvData.hpp>  // for EnvData, etc
#include <stk_util/environment/ProgramOptions.hpp>
#include <stk_util/environment/RuntimeMessage.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>  // for all_write_string
#include <string>                       // for string, operator<<, etc
#include "boost/program_options/detail/parsers.hpp"
#include "boost/program_options/errors.hpp"  // for program_options
#include "boost/program_options/variables_map.hpp"  // for variables_map, etc
#include "stk_util/environment/ReportHandler.hpp"  // for ThrowRequire

#if defined(__GNUC__)
#include <cstdlib>
#include <sys/time.h>
#include <sys/resource.h>
#endif

namespace boost { namespace program_options { class options_description; } }





using namespace std;

namespace sierra {

std::string
format_time(
  double	t,
  const char *	format)
{
  time_t	time = static_cast<time_t>(t);
  char s[128];

  ::strftime(s, sizeof(s), format, ::localtime(&time));

  return std::string(s);
}

namespace Env {

  //
  //  Set or get the gemini version, if passed value is not unknown, set the version, either way return the version
  //
  GeminiSCIVersion GetGeminiVersion(GeminiSCIVersion ver) {
    static GeminiSCIVersion GeminiSCIVersionValue = GEMINI_SCI_1;  //This is the default gemini version
    if(ver != GEMINI_SCI_UNKNOWN) {
      GeminiSCIVersionValue = ver;
    }
    ThrowRequire(GeminiSCIVersionValue != GEMINI_SCI_UNKNOWN);
    return GeminiSCIVersionValue;
  }
    


const std::string &
product_name()
{
  return stk::EnvData::instance().m_productName;
}


const std::string &
executable_file()
{
  return stk::EnvData::instance().m_executablePath;
}


const std::string &
startup_date()
{
  static std::string startup_date;

  if (startup_date.empty())
    startup_date = format_time(stk::EnvData::instance().m_startTime).c_str();

  return startup_date;
}


double
start_time()
{
  return stk::EnvData::instance().m_startTime;
}


bool
developer_mode()
{
  return !get_param("developer-mode").empty();
}


// Similar to Platform.cpp's get_heap_info, but no 'largest_free' and no log output.
void get_heap_used(size_t &heap_size)
{
  heap_size = 0;

#if defined(SIERRA_HEAP_INFO)

# if defined(SIERRA_PTMALLOC3_ALLOCATOR) || defined(SIERRA_PTMALLOC2_ALLOCATOR)
  heap_size = malloc_used();
  
# elif ( defined(__linux__) || defined(REDS) ) && ! defined(__IBMCPP__)
  static struct mallinfo minfo;
  minfo = mallinfo();
  heap_size = static_cast<unsigned int>(minfo.uordblks) + static_cast<unsigned int>(minfo.hblkhd);

# elif defined(__sun)
  pstatus_t proc_status;

  std::ifstream proc("/proc/self/status", std::ios_base::in|std::ios_base::binary);
  if (proc) {
    proc.read(reinterpret_cast<char *>(&proc_status), sizeof(proc_status));
    heap_size = proc_status.pr_brksize;
  }
# endif
#endif // defined(SIERRA_HEAP_INFO)
}

void setInputFileName(std::string name) {
  stk::EnvData::instance().m_inputFile = name;
}

std::string getInputFileName() {
  return stk::EnvData::instance().m_inputFile;
}

void set_input_file_required(bool value)
{
    stk::EnvData::instance().m_inputFileRequired = value;
}

void set_sm_preprocessing(bool value)
{
    stk::EnvData::instance().m_checkSubCycle = value;
    stk::EnvData::instance().m_checkSmRegion = value;
}

void set_zapotec(bool value)
{
    stk::EnvData::instance().m_isZapotec = value;
}

bool is_zapotec()
{
    return stk::EnvData::instance().m_isZapotec;
}


const std::string &
architecture()
{
  return get_param("architecture");
}


const std::string
working_directory() {
  char cwd[PATH_MAX];
  std::string directory = get_param("directory");
  if (directory[0] != '/' && getcwd(cwd, PATH_MAX) != NULL) {
    directory = cwd;
    directory += '/';
  }
  return directory;
}


std::ostream &
output()
{
  return stk::EnvData::instance().m_output;
}


std::ostream &
outputP0()
{
  return *stk::EnvData::instance().m_outputP0;
}


std::ostream &
outputNull() {
  return stk::EnvData::instance().m_outputNull;
}


const char *
section_separator()
{
  static const char *s_sectionSeparator = "+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----";

  return s_sectionSeparator;
}


const char *
subsection_separator()
{
  static const char *s_subsectionSeparator = "---------------------------------------------------";

  return s_subsectionSeparator;
}


double
wall_now()
{
  timeval tp;
  struct timezone tz;
  gettimeofday(&tp, &tz);
  return (tp.tv_sec + ((static_cast<double>(tp.tv_usec))/1000000.0));
}


double
cpu_now()
{
#if defined(REDS)
  struct rusage my_rusage;

  getrusage(RUSAGE_SELF, &my_rusage);

  return static_cast<double>(my_rusage.ru_utime.tv_sec) +
    static_cast<double>(my_rusage.ru_utime.tv_usec)*1.0e-6;

#elif ! defined(__PGI)
  struct rusage my_rusage;

  getrusage(RUSAGE_SELF, &my_rusage);

  return static_cast<double>(my_rusage.ru_utime.tv_sec + my_rusage.ru_stime.tv_sec) +
    static_cast<double>(my_rusage.ru_utime.tv_usec + my_rusage.ru_stime.tv_usec)*1.0e-6;
#else
  return 0;
#endif
}


std::string
section_title(
  const std::string &	title)
{
  static size_t s_sectionSeparatorLength = std::strlen(section_separator());

  std::ostringstream strout;

  strout << std::left << std::setw(s_sectionSeparatorLength - 20) << title << std::right << std::setw(20) << format_time(Env::wall_now());
  return strout.str();
}


int parallel_size() {
  return stk::EnvData::instance().m_parallelSize;
}

int parallel_rank() {
  return stk::EnvData::instance().m_parallelRank;
}

MPI_Comm
parallel_comm()
{
  return stk::EnvData::instance().m_parallelComm;
}

MPI_Comm
parallel_intercomm()
{
  return stk::EnvData::instance().m_interComm;
}

MPI_Comm
parallel_world_comm()
{
  return stk::EnvData::instance().m_worldComm;
}

int parallel_lag_master() {
  return stk::EnvData::instance().m_execMap[EXEC_TYPE_LAG].m_master;
}

int parallel_fluid_master() {
  return stk::EnvData::instance().m_execMap[EXEC_TYPE_FLUID].m_master;
}

int peer_group() {
  return stk::EnvData::instance().m_execMap[EXEC_TYPE_PEER].m_master;
}

bool
is_comm_valid()
{
  stk::EnvData &env_data = stk::EnvData::instance();
  if (env_data.m_parallelComm == MPI_COMM_NULL) {
    return false;
  } else {
    return true;
  }
}

void
output_flush()
{
  stk::EnvData &env_data = stk::EnvData::instance();

  stk::report_deferred_messages(Env::parallel_comm());

  stk::all_write_string(Env::parallel_comm(), *env_data.m_outputP0, env_data.m_output.str());
  env_data.m_output.str("");
}


void
request_shutdown(bool shutdown)
{
  stk::EnvData::instance().m_shutdownRequested = shutdown;
}


bool
is_shutdown_requested()
{
  int shutdown_requested_in = stk::EnvData::instance().m_shutdownRequested || Env::HUP_received();
  int shutdown_requested = -1;

#if defined(STK_HAS_MPI)
  MPI_Allreduce(&shutdown_requested_in, &shutdown_requested, 1, MPI_INT, MPI_SUM, Env::parallel_comm());
#else
  shutdown_requested = shutdown_requested_in;
#endif

  return shutdown_requested != 0;
}


void abort() {
  stk::EnvData &env_data = stk::EnvData::instance();

  // Cannot be sure of parallel synchronization status; therefore, no communications can
  // occur.  Grab and dump all pending output buffers to 'std::cerr'.
  std::cerr << std::endl
            << "*** SIERRA ABORT on P" << stk::EnvData::instance().m_parallelRank << " ***"
            << std::endl
            << "*** check " << get_param("output-log")
            << " file for more information ***"
            << std::endl ;

  if (!env_data.m_output.str().empty()) {
    std::cerr << "Buffer contents of deferred output stream on processor " << parallel_rank()
              << std::endl ;
    std::cerr << env_data.m_output.str();
  }

  std::cerr.flush();
  std::cout.flush();

  ::sleep(1);					// Give the other processors a chance at
						// catching up, seems to help hanging problems.
#if defined(STK_HAS_MPI)
  MPI_Abort(env_data.m_parallelComm, MPI_ERR_OTHER);	// First try to die
#endif
  std::exit( EXIT_FAILURE );                    // Second try to die
}


const std::string &
get_param(
  const char * const	option)
{
  if (stk::EnvData::instance().m_vm.count(option)) {
    if (stk::EnvData::instance().m_vm[option].as<std::string>().empty())
      return stk::EnvData::instance().m_onString;
    else
      return stk::EnvData::instance().m_vm[option].as<std::string>();
  }
  else
    return stk::EnvData::instance().m_emptyString;
}


void
set_param(
  const char *          option,
  const std::string &   value) {


  namespace opt = boost::program_options;

  opt::variables_map &vm = stk::get_variables_map();
  opt::options_description &od = stk::get_options_description();

  int argc = 1;
  char *s = std::strcpy(new char[std::strlen(option) + 1], option);

  opt::store(opt::parse_command_line(argc, &s, od), vm);
  opt::notify(vm);

  delete [] s;
}

} // namespace Env
} // namespace sierra
