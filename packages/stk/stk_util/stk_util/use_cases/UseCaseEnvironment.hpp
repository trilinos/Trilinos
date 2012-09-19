/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_util_use_cases_UseCaseEnvironement_hpp
#define stk_util_use_cases_UseCaseEnvironement_hpp

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/BroadcastArg.hpp>

#include <stk_util/diag/Writer_fwd.hpp>
#include <stk_util/diag/Timer.hpp>

#include <stk_util/environment/OutputLog.hpp>
#include <stk_util/environment/ProgramOptions.hpp>
#include <stk_util/environment/ReportHandler.hpp>
#include <stk_util/environment/RuntimeWarning.hpp>
#include <stk_util/environment/RuntimeDoomed.hpp>

#include <iosfwd>

namespace use_case {

enum LogMask {
  LOG_ALWAYS		= stk::LOG_ALWAYS,
  LOG_TRACE		= stk::LOG_TRACE,
  LOG_TRACE_STATS	= stk::LOG_TRACE_STATS,
  LOG_TRACE_SUB_CALLS	= stk::LOG_TRACE_SUB_CALLS,
  LOG_MEMBERS		= stk::LOG_MEMBERS,

  LOG_SEARCH            = 0x0000010,
  LOG_TRANSFER          = 0x0000020,
  LOG_TIMER             = 0x0000040
};

/**
 * @brief Class <b>message_type</b> ...
 *
 */
enum message_type {
  MSG_WARNING = stk::MSG_WARNING,
  MSG_FATAL = stk::MSG_DOOMED,
  MSG_INFORMATION,
  MSG_EXCEPTION,
  MSG_PARALLEL_EXCEPTION
};


/**
 * @brief Class <b>type</b> ...
 *
 */
enum message_throttle_type {
  MSG_APPLICATION = stk::MSG_APPLICATION,
  MSG_TIME_STEP
};

enum TimerSetMask {
  TIMER_MESH		= 0x00000001,		///< Enable mesh timers
  TIMER_MESH_IO		= 0x00000002,		///< Enable mesh I/O timers
  TIMER_SEARCH		= 0x00000004,		///< Enable search timers
  TIMER_TRANSFER	= 0x00000008,		///< Enable transfer timers
  TIMER_ALL		= 0xFFFFFFFF,		///< Force timer to be active

  TIMER_FORCE		= 0x00000000		///< Force timer to be active
};

std::ostream &out();                ///< Normal output stream
std::ostream &dout();               ///< Diagnostic output stream
std::ostream &pout();               ///< Per-processor output stream (See RuntimeDeferredx)
std::ostream &tout();               ///< Regression test textual output stream

std::ostream &dwout();              ///< Diagnostic writer stream

stk::diag::Writer &dw();

stk::diag::TimerSet &timerSet();

stk::diag::Timer &timer();

void my_report_handler(const char *message, int type);

/**
 * Print the sucess or failure of a use_case by combining the success
 * flags for all procs in comm. Prints a string that should be grepped for
 * in the test XML files.
 *
 * Returns true if all procs passed.
 */
bool print_status(stk::ParallelMachine comm, bool success);

struct UseCaseEnvironment
{
  // Will initialize a comm
  UseCaseEnvironment(int *argc, char ***argv);

  // Assumes already-initialized comm
  UseCaseEnvironment(int *argc, char ***argv, stk::ParallelMachine comm);

  // shared constructor implementation; do not call directly
  void initialize(int* argc, char ***argv);

  ~UseCaseEnvironment();

  const stk::ParallelMachine    m_comm;
  std::string                   m_workingDirectory;
  bool                          m_need_to_finalize;
};

} // namespace use_case

#endif // stk_util_use_cases_UseCaseEnvironement_hpp
