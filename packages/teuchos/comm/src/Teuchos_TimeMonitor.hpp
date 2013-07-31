// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef TEUCHOS_TIMEMONITOR_HPP
#define TEUCHOS_TIMEMONITOR_HPP


/*! \file Teuchos_TimeMonitor.hpp
 *
 * \brief Scope protection wrapper for Teuchos::Time, with timer reporting functionality.
 *
 * An instance of the Teuchos::TimeMonitor class wraps a nonconst
 * reference to a Teuchos::Time timer object.  TimeMonitor's
 * constructor starts the timer, and its destructor stops the timer.
 * This ensures scope safety of timers, so that no matter how a scope
 * is exited (whether the normal way or when an exception is thrown),
 * a timer started in the scope is stopped when the scope is left.
 *
 * TimeMonitor also has class methods that create or destroy timers
 * (in such a way that it can track the complete set of created timers
 * on each process) and compute global timer statistics.
 */

/** \example TimeMonitor/cxx_main.cpp
 *
 * This is an example of how to use the Teuchos::TimeMonitor class.
 */

#include "Teuchos_PerformanceMonitorBase.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_Time.hpp"

#include "Teuchos_CommandLineProcessor.hpp"

//! Defines a static non-member function that returns a time monitor.
#define TEUCHOS_TIMER(funcName, strName) \
  static Teuchos::Time& funcName() \
  {static Teuchos::RCP<Time> rtn = \
      Teuchos::TimeMonitor::getNewCounter(strName); return *rtn;}


/** \brief Defines a timer for a specific function (with differentiator).
 *
 * Same as TEUCHOS_FUNC_TIME_MONITOR(...) except required when used more than
 * once in the same function (like a block of code).
 */
#define TEUCHOS_FUNC_TIME_MONITOR_DIFF( FUNCNAME, DIFF ) \
  static Teuchos::RCP<Teuchos::Time> DIFF ## blabla_localTimer; \
  if(!DIFF ## blabla_localTimer.get()) { \
    std::ostringstream oss; \
    oss << FUNCNAME; \
    DIFF ## blabla_localTimer = Teuchos::TimeMonitor::getNewCounter(oss.str()); \
  } \
  Teuchos::TimeMonitor DIFF ## blabla_localTimeMonitor(*DIFF ## blabla_localTimer)


/** \brief Defines a timer for a specific function.
 *
 Note that the name of the timer can be formated with stream inserts.
 For example, we can define a time monitor for a function as follows:

 \code
 template<typename Scalar>
 void foo()
 {
   TEUCHOS_FUNC_TIME_MONITOR(
     "foo<" << Teuchos::ScalarTraits<Scalar>::name () << ">()"
     );
   ...
 }
 \endcode

 The timer can then be printed at the end of the program using any of
 various class methods, including summarize():
 \code
 Teuchos::TimeMonitor::summarize ();
 \endcode
*/
#define TEUCHOS_FUNC_TIME_MONITOR( FUNCNAME ) \
  TEUCHOS_FUNC_TIME_MONITOR_DIFF( FUNCNAME, main )


namespace Teuchos {

/// \typedef stat_map_type
/// \brief Global statistics collected from timer data.
///
/// Key: name of the timer.
///
/// Value: each entry in the vector is a timing and call count for
///   that timer, corresponding to a particular statistic (e.g.,
///   minimum, arithmetic mean, or maximum).  What statistic that is
///   depends on an auxillary array "statNames" which has the same
///   ordering as the entries in this vector.  See the documentation
///   of \c TimeMonitor::computeGlobalTimerStatistics().
typedef std::map<std::string, std::vector<std::pair<double, double> > > stat_map_type;

/// \class TimeMonitor
/// \brief A scope-safe timer wrapper class, that can compute global timer statistics.
///
/// An instance of the \c TimeMonitor class wraps a nonconst reference
/// to a \c Time timer object.  \c TimeMonitor's constructor starts
/// the timer, and its destructor stops the timer.  This ensures scope
/// safety of timers, so that no matter how a scope is exited (whether
/// the normal way or when an exception is thrown), a timer started in
/// the scope is stopped when the scope is left.
///
/// \c TimeMonitor also has class methods that create or destroy
/// timers and compute global timer statistics.  If you create a timer
/// using getNewCounter() (or the deprecated getNewTimer()), it will
/// add that timer to the set of timers for which to compute global
/// statistics.  The summarize() and report() methods will print
/// global statistics for these timers, like the minimum, mean, and
/// maximum time over all processes in the communicator, for each
/// timer.  These methods work correctly even if some processes have
/// different timers than other processes.  You may also use \c
/// computeGlobalTimerStatistics() to compute the same global
/// statistics, if you wish to use them in your program or output them
/// in a different format than that of these methods.
///
/// \warning This class must only be used to time functions that are
///   called only within the main program.  It may _not_ be used in
///   pre-program setup or post-program teardown!
class TEUCHOSCOMM_LIB_DLL_EXPORT TimeMonitor :
    public PerformanceMonitorBase<Time> {
public:

  /** \name Constructor/Destructor */
  //@{

  /// \brief Constructor: starts the timer.
  ///
  /// \param timer [in/out] Reference to the timer to be wrapped.
  ///   This constructor starts the timer, and the destructor stops
  ///   the timer.
  ///
  /// \param reset [in] If true, reset the timer before starting it.
  ///   Default behavior is not to reset the timer.
  TimeMonitor (Time& timer, bool reset=false);

  //! Destructor: stops the timer.
  ~TimeMonitor();
  //@}

  /// \brief Return a new timer with the given name (class method).
  ///
  /// Call getNewCounter() or this method if you want to create a new
  /// named timer, and you would like TimeMonitor to track the timer
  /// for later computation of global statistics over processes.
  ///
  /// This method wraps getNewCounter() (inherited from the base
  /// class) for backwards compatibiity.
  static RCP<Time> getNewTimer (const std::string& name) {
    return getNewCounter (name);
  }

  /// \brief Disable the timer with the given name.
  ///
  /// "Disable" means that the timer (Time instance) will ignore all
  /// calls to start(), stop(), and incrementNumCalls().  The effect
  /// will be as if the TimeMonitor had never touched the timer.
  ///
  /// If the timer with the given name does not exist (was never
  /// created using getNewCounter() or getNewTimer()), then this
  /// method throws std::invalid_argument.  Otherwise, it disables the
  /// timer.  This effect lasts until the timer is cleared or until
  /// the timer is enabled, either by calling enableTimer() (see
  /// below) or by calling the Time instance's enable() method.
  ///
  /// Disabling a timer does <i>not</i> exclude it from the list of
  /// timers printed by summarize() or report().
  static void disableTimer (const std::string& name);

  /// \brief Enable the timer with the given name.
  ///
  /// If the timer with the given name does not exist (was never
  /// created using getNewCounter() or getNewTimer()), then this
  /// method throws std::invalid_argument.  Otherwise, it undoes the
  /// effect of disableTimer() on the timer with the given name.  If
  /// the timer with the given name was not disabled, then this method
  /// does nothing.
  static void enableTimer (const std::string& name);

  /// \brief Reset all global timers to zero.
  ///
  /// This method only affects Time objects created by getNewCounter()
  /// or getNewTimer().
  ///
  /// \pre None of the timers must currently be running.
  static void zeroOutTimers();

  /// \brief Compute global timer statistics for all timers on the given communicator.
  ///
  /// The typical use case for Time and TimeMonitor is that all
  /// processes in a communicator create the same set of timers, and
  /// then want to report summary statistics.  This method supports
  /// that typical use case.  For each timer in the set, this method
  /// computes a list of global statistics.  "Global" means "for all
  /// processes in the communicator."  "Statistic" means the result of
  /// a reduction over the timing and call count values.  Thus, each
  /// statistic includes both a timing and a call count.  The current
  /// list of computed statistics includes the minimum and maximum
  /// timing (and the corresponding call count for each) and the
  /// arithmetic mean (timing and call count).  This list may expand
  /// in the future.
  ///
  /// Different processes may have different sets of timers.  This
  /// method gives you two options for reconciling the sets.  If setOp
  /// is Intersection, it computes the intersection (the common
  /// subset) of timers on all processes in the communicator.
  /// Otherwise, if setOp is Union, this method computes the union of
  /// timers on all processes in the communicator.  Intersection is
  /// the default, since it means that all reported timers exist on
  /// all participating processes.  For setOp=Union, timers that do
  /// not exist on some processes will be given a zero timing and call
  /// count, so that statistics make sense.
  ///
  /// \note This method must called as a collective by all processes
  ///   in the communicator.
  ///
  /// All output arguments are returned redundantly on all processes
  /// in the communicator.  That makes this method an all-reduce.
  ///
  /// \section Teuchos_TimeMonitor_computeGlobalTimerStatistics_stats Statistics collected
  ///
  /// The "MinOverProcs" and "MaxOverProcs" timings are cumulative:
  /// the reported timing is for all calls.  Along with the min resp.
  /// max timing comes the call count of the process who had the min
  /// resp. max.  (If more than one process had the min resp. max
  /// timing, then the call count on the process with the smallest
  /// rank is reported.)
  ///
  /// The "MeanOverProcs" equals the sum of the processes' cumulative
  /// timings, divided by the number of processes.  Thus, it is
  /// cumulative over all calls, and is comparable with the
  /// "MinOverProcs" and "MaxOverProcs" timings.  This differs from
  /// the "MeanOverCallCounts" (see below).  This does <i>not</i>
  /// weight the mean by call counts.
  ///
  /// The "MeanOverCallCounts" is an arithmetic mean of all timings.
  /// It is <i>not</i> cumulative.  It reports the mean timing for a
  /// single invocation over all calls on all processes, not weighting
  /// any one process more than the others.  For each timer, this is
  /// the sum of the cumulative timing over all processes, divided by
  /// the sum of the call counts over all processes for that timing.
  /// (We compute it a bit differently to help prevent overflow.)  The
  /// "MeanOverCallCounts" is <i>not</i> comparable with the min, max,
  /// or "MeanOverProcs".
  ///
  /// We report with both versions of the mean timing the mean call
  /// count over processes.  This may be fractional, which is one
  /// reason why we report call counts as \c double rather than \c
  /// int.  It has no particular connection to the mean timing.
  ///
  /// \section Teuchos_TimeMonitor_computeGlobalTimerStatistics_perf Performance
  ///
  /// This operation requires interprocess communication.  Suppose
  /// there are \f$P\f$ processes in the given communicator, and
  /// \f$N\f$ unique timers in the global union of all processes'
  /// timers.  Then, this method requires \f$O(\log P)\f$ messages
  /// (\f$O(1)\f$ "reductions" and exactly 1 "broadcast") and
  /// \f$O(N)\f$ per-processor storage (in the worst case) when
  /// computing either the intersection or the union of timers (the
  /// algorithm is similar in either case).  The whole algorithm takes
  /// at worst \f$O(N (\log N) (\log P))\f$ time along the critical
  /// path (i.e., on the "slowest process" in the communicator).  The
  /// \f$N \log N\f$ term comes from sorting the timers by label at
  /// each stage of the reduction in order to compute their union or
  /// intersection.
  ///
  /// \param statData [out] On output: Global timer statistics, stored
  ///   as a map with key timer name, and with value the ordered list
  ///   of statistics for that timer.  The \c statNames output has the
  ///   same order as the ordered list of statistics for each timer.
  ///   Each entry of the statistics list is a (timing, call count)
  ///   pair, the meaning of which depends on the particular statistic
  ///   (see above).
  ///
  /// \param statNames [out] On output: Each value in the statData map
  ///   is a vector.  That vector v has the same number of entries as
  ///   statNames.  statNames[k] is the name of the statistic (see
  ///   above) stored as v[k].  Always refer to statNames for the
  ///   number and names of statistics.
  ///
  /// \param comm [in] Communicator whose process(es) will participate
  ///   in the gathering of timer statistics.  This is a Ptr and not
  ///   an RCP, because RCP would suggest that TimeMonitor were
  ///   keeping the communicator around after return of this method.
  ///   Ptr suggests instead that TimeMonitor will only reference the
  ///   communicator during this method.  If you have an RCP, you can
  ///   turn it into a Ptr by calling its ptr() method:
  ///   \code
  ///   RCP<const Comm<int> > myComm = ...;
  ///   TimeMonitor::computeGlobalTimerStatistics (statData, statNames, myComm.ptr());
  ///   \endcode
  ///
  /// \param setOp [in] If \c Intersection, compute statistics for the
  ///   intersection of all created timers over all processes in the
  ///   communicator.  If \c Union, compute statistics for the union
  ///   of all created timers over all processes in the communicator.
  ///
  /// \param filter [in] Filter for timer labels.  If filter is not
  ///   empty, this method will only compute statistics for timers
  ///   whose labels begin with this string.
  static void
  computeGlobalTimerStatistics (stat_map_type& statData,
                                std::vector<std::string>& statNames,
                                Ptr<const Comm<int> > comm,
                                const ECounterSetOp setOp=Intersection,
                                const std::string& filter="");

  /// \brief Compute global timer statistics for all timers on all (MPI) processes.
  ///
  /// This is an overload of the above computeGlobalTimerStatistics()
  /// method for when the caller does not want to provide a
  /// communicator explicitly.  This method "does the right thing" in
  /// that case.  Specifically:
  /// - If Trilinos was not built with MPI support, this method
  ///   assumes a serial "communicator" containing one process.
  /// - If Trilinos was built with MPI support and MPI has been
  ///   initialized (via MPI_Init() or one of the wrappers in
  ///   Epetra or Teuchos), this method uses MPI_COMM_WORLD as the
  ///   communicator.  This is the most common case.
  /// - If Trilinos was built with MPI support and MPI has <i>not</i>
  ///   been initialized, this method will use a "serial" communicator
  ///   (that does not actually use MPI).  This may produce output on
  ///   all the MPI processes if you are running with Trilinos as an
  ///   MPI job with more than one process.  Thus, if you intend to
  ///   use this method in parallel, you should first initialize MPI.
  ///   (We cannot initialize MPI for you, because we have no way to
  ///   know whether you intend to run an MPI-enabled build serially.)
  ///
  /// \warning If you call this method when MPI is running, you
  ///   <i>must</i> call it on all processes in \c MPI_COMM_WORLD.
  ///   Otherwise, the method will never finish, since it will be
  ///   waiting forever for the non-participating processes.  If you
  ///   want to use computeGlobalTimerStatistics() on a
  ///   subcommunicator, please use the overloaded version above that
  ///   takes a communicator as an input argument.
  static void
  computeGlobalTimerStatistics (stat_map_type& statData,
                                std::vector<std::string>& statNames,
                                const ECounterSetOp setOp=Intersection,
                                const std::string& filter="");

  /// \brief Print summary statistics for all timers on the given communicator.
  ///
  /// If writeGlobalStatus=true, this method computes the same
  /// statistics as computeGlobalTimerStatistics(), using the same
  /// collective algorithm.  (<tt>writeGlobalStatus=false</tt> means
  /// that only the process with rank 0 in the communicator reports
  /// its timers' data.)  It then reports the results to the given
  /// output stream on the process with rank 0 in the given
  /// communicator.  Output follows a human-readable tabular form.
  ///
  /// \param comm [in] Communicator whose process(es) will participate
  ///   in the gathering of timer statistics.  This is a Ptr and not
  ///   an RCP, because RCP would suggest that TimeMonitor were
  ///   keeping the communicator around after return of this method.
  ///   Ptr suggests instead that TimeMonitor will only reference the
  ///   communicator during this method.  If you have an RCP, you can
  ///   turn it into a Ptr by calling its ptr() method:
  ///   \code
  ///   RCP<const Comm<int> > myComm = ...;
  ///   TimeMonitor::summarize (myComm.ptr());
  ///   \endcode
  ///
  /// \param out [out] Output stream to which to write.  This will
  ///   only be used on the process with rank 0 in the communicator.
  ///
  /// \param alwaysWriteLocal [in] If true, the process with Rank 0 in
  ///   the communicator will write its local timings to the given
  ///   output stream.  Defaults to false, since the global statistics
  ///   are more meaningful.  If the local set of timers differs from
  ///   the global set of timers (either the union or the
  ///   intersection, depending on \c setOp), Proc 0 will create
  ///   corresponding local timer data (<i>not</i> corresponding
  ///   timers) with zero elapsed times and call counts, just to pad
  ///   the table of output.
  ///
  /// \param writeGlobalStats [in] If true (the default), compute and
  ///   display the statistics that \c computeGlobalTimerStatistics()
  ///   computes.  If there is only one MPI process or if this is a
  ///   non-MPI build of Trilinos, only compute and show the "global"
  ///   timings, without the "statistics" that would be all the same
  ///   anyway.
  ///
  /// \param writeZeroTimers [in] If false, do not display results for
  ///   timers that have never been called (numCalls() == 0).  If
  ///   true, display results for all timers, regardless of their call
  ///   count.  Note that \c setOp and \c writeGlobalStats might
  ///   reintroduce timers with zero call counts.
  ///
  /// \param setOp [in] If \c Intersection, compute and display the
  ///   intersection of all created timers over all processes in the
  ///   communicator.  If \c Union, compute and display the union of
  ///   all created timers over all processes in the communicator.
  ///
  /// \param filter [in] Filter for timer labels.  If filter is not
  ///   empty, this method will only print timers whose labels begin
  ///   with this string.
  ///
  /// \note If \c writeGlobalStats is true, this method <i>must</i> be
  ///   called as a collective by all processes in the communicator.
  ///   This method will <i>only</i> perform communication if
  ///   <tt>writeGlobalStats</tt> is true.
  static void
  summarize (Ptr<const Comm<int> > comm,
             std::ostream &out=std::cout,
             const bool alwaysWriteLocal=false,
             const bool writeGlobalStats=true,
             const bool writeZeroTimers=true,
             const ECounterSetOp setOp=Intersection,
             const std::string& filter="");

  /// \brief Print summary statistics for all timers on all (MPI) processes.
  ///
  /// This is an overload of the above summarize() method for when the
  /// caller does not want to provide a communicator explicitly.  This
  /// method "does the right thing" in that case.  For an explanation
  /// of what that means, see the documentation of the overload of
  /// computeGlobalTimerStatistics() that does not require a
  /// communicator argument.
  ///
  /// \warning If you call this method when MPI is running, you
  ///   <i>must</i> call it on all processes in \c MPI_COMM_WORLD.
  ///   Otherwise, the method will never finish, since it will be
  ///   waiting forever for the non-participating processes.  If you
  ///   want to use \c summarize() on a subcommunicator, please use
  ///   the overloaded version above that takes a communicator as an
  ///   input argument.
  static void
  summarize (std::ostream& out=std::cout,
             const bool alwaysWriteLocal=false,
             const bool writeGlobalStats=true,
             const bool writeZeroTimers=true,
             const ECounterSetOp setOp=Intersection,
             const std::string& filter="");

  /// \brief Report timer statistics to the given output stream.
  ///
  /// This is like summarize(), but gives you more control over the
  /// output format.  To get the default parameters, either call
  /// getValidReportParameters(), or call this method with params
  /// nonnull but empty (it will fill in default parameters).
  ///
  /// \param comm [in] Communicator whose process(es) will participate
  ///   in the gathering of timer statistics.  This is a Ptr and not
  ///   an RCP, because RCP would suggest that TimeMonitor were
  ///   keeping the communicator around after return of this method.
  ///   Ptr suggests instead that TimeMonitor will only reference the
  ///   communicator during this method.  If you have an RCP, you can
  ///   turn it into a Ptr by calling its ptr() method:
  ///   \code
  ///   RCP<const Comm<int> > myComm = ...;
  ///   TimeMonitor::report (myComm.ptr (), ...);
  ///   \endcode
  ///
  /// \param out [out] Output stream to which to write.  This will
  ///   only be used on the process with rank 0 in the communicator.
  ///
  /// \param filter [in] Filter for timer labels.  If filter is not
  ///   empty, this method will only print timers whose labels begin
  ///   with this string.
  ///
  /// \param params [in/out] Parameters to control output format and
  ///   which statistics to generate.  If null, we use default
  ///   parameters if this method was not yet called with params
  ///   nonnull, otherwise we use the previous set of parameters.  If
  ///   nonnull, we read the given parameters, filling in defaults,
  ///   and use the resulting parameters for all subsequent calls to
  ///   report() (until new parameters are set).
  ///
  /// \section Teuchos_TimeMonitor_report_SupportedParams Supported parameters
  ///
  /// Here is the current set of supported parameters:
  /// - "Report format": "Table" (default), "YAML"
  /// - "YAML style": "spacious" (default), "compact"
  /// - "How to merge timer sets": "Intersection" (default), "Union"
  /// - "alwaysWriteLocal": true, false (default)
  /// - "writeGlobalStats": true (default), false
  /// - "writeZeroTimers": true (default), false
  ///
  /// This method currently supports two different output formats.
  /// "Table" format is the same tabular format which summarize()
  /// uses.  It displays times and call counts in a table that is easy
  /// for humans to read, but hard to parse.  "YAML" format uses a
  /// standard, structured, human-readable output format called YAML.
  /// <a href="http://yaml.org">YAML</a> stands for YAML Ain't Markup
  /// Language.
  ///
  /// "YAML style" refers to two variants of YAML output that report()
  /// can generate.  The "compact" mode attempts to put as much data
  /// on each line as possible.  It may be more readable when there
  /// are a small number of timers.  The "spacious" mode prefers one
  /// line per datum whenever possible.  Both modes have the same
  /// schema, that is, their output has the same hierarchical
  /// structure and thus the same parse tree.
  ///
  /// (In technical terms: compact mode uses YAML's so-called "flow
  /// style" for sequences and mappings whenever possible, except at
  /// the outermost level where it would hinder readability.  Spacious
  /// mode does not use "flow style" for lists or mappings.  For an
  /// explanation of YAML's flow style, see <a
  /// href="http://www.yaml.org/spec/1.2/spec.html#style/flow/">Chapter
  /// 7 of the YAML 1.2 spec</a>.)
  ///
  /// "How to merge timer sets" refers to the set operation by which
  /// processors should combine their sets of timers in order to
  /// compute global timer statistics.  This corresponds to the
  /// <tt>setOp</tt> argument of summarize().
  ///
  /// The remaining Boolean parameters are the same as the eponymous
  /// arguments of summarize(), to whose documentation one should
  /// refer.  There are some wrinkles: in particular, YAML output
  /// ignores the "alwaysWriteLocal" parameter and assumes
  /// "writeGlobalStats" is true.
  static void
  report (Ptr<const Comm<int> > comm,
          std::ostream& out,
          const std::string& filter,
          const RCP<ParameterList>& params=null);

  /// \brief Report timer statistics to the given output stream.
  ///
  /// This is like the 4-argument version of report(), but with a
  /// default filter.
  static void
  report (Ptr<const Comm<int> > comm,
          std::ostream& out,
          const RCP<ParameterList>& params=null);

  /// \brief Report timer statistics to the given output stream.
  ///
  /// This is like the 4-argument version of report(), but with a
  /// default communicator.
  static void
  report (std::ostream& out,
          const std::string& filter,
          const RCP<ParameterList>& params=null);

  /// \brief Report timer statistics to the given output stream.
  ///
  /// This is like the 4-argument version of report(), but with a
  /// default communicator and a default filter.
  static void
  report (std::ostream& out,
          const RCP<ParameterList>& params=null);

  //! Default parameters (with validators) for report().
  static RCP<const ParameterList> getValidReportParameters ();

 private:
  /// \brief Valid output formats for report().
  ///
  /// \warning This is an implementation detail of TimeMonitor.  It is
  ///   subject to change at any time without notice.
  enum ETimeMonitorReportFormat {
    REPORT_FORMAT_YAML,
    REPORT_FORMAT_TABLE
  };

  /// \brief Valid YAML output formats for report().
  ///
  /// \warning This is an implementation detail of TimeMonitor.  It is
  ///   subject to change at any time without notice.
  enum ETimeMonitorYamlFormat {
    YAML_FORMAT_COMPACT,
    YAML_FORMAT_SPACIOUS
  };

  /// \brief Like summarize(), but with YAML-format output.
  ///
  /// \param comm [in] Communicator over which to compute timer
  ///   statistics.
  /// \param out [out] Output stream to which to write (on Proc 0 of
  ///   the given communicator only).
  /// \param yamlStyle [in] Whether to print YAML output in "compact"
  ///   or "spacious" style.
  /// \param filter [in] Filter for timer labels.  If filter is not
  ///   empty, this method will only print timers whose labels begin
  ///   with this string.
  ///
  /// \warning This is an experimental interface.  It may change or
  ///   disappear without warning.
  static void
  summarizeToYaml (Ptr<const Comm<int> > comm,
                   std::ostream& out,
                   const ETimeMonitorYamlFormat yamlStyle,
                   const std::string& filter="");

  /// \brief Like summarize(), but with YAML-format output and default communicator.
  ///
  /// \warning This is an experimental interface.  It may change or
  ///   disappear without warning.
  static void
  summarizeToYaml (std::ostream& out,
                   const ETimeMonitorYamlFormat yamlStyle,
                   const std::string& filter="");

  /// \brief Add the "Report format" parameter to plist.
  ///
  /// \note Call this in getValidReportParameters() to set a default
  ///   value and validator for this parameter.
  static void setReportFormatParameter (ParameterList& plist);

  /// \brief Add the "YAML style" parameter to plist.
  ///
  /// \note Call this in getValidReportParameters() to set a default
  ///   value and validator for this parameter.
  static void setYamlFormatParameter (ParameterList& plist);

  /// \brief Add the "How to merge timer sets" parameter to plist.
  ///
  /// \note Call this in getValidReportParameters() to set a default
  ///   value and validator for this parameter.
  static void setSetOpParameter (ParameterList& plist);

  /// \brief Set parameters for report().  Call only from report().
  ///
  /// If this method completes successfully, it sets setParams_ to
  /// true as a flag.
  ///
  /// \param params [in/out] Parameters for report().  This may be
  ///   null, in which case we use defaults or the last set of
  ///   parameters.
  ///
  /// \warning This method is not thread safe, in the sense that it
  ///   does not set the class data atomically.  Behavior when calling
  ///   this method from multiple threads is undefined.  Calling this
  ///   routine with different parameter lists from different threads
  ///   will certainly not accomplish what you want to accomplish.
  static void setReportParameters (const RCP<ParameterList>& params);

  //! Parameters for the report() class method.
  //@{

  //! Current output format for report().  Set via setReportParameters().
  static ETimeMonitorReportFormat reportFormat_;

  /// Current output style for report(), when using YAML output.
  /// Set via setReportParameters().
  static ETimeMonitorYamlFormat yamlStyle_;

  //! Whether report() should use the intersection or union of timers over processes.
  static ECounterSetOp setOp_;

  //! Whether report() should always report Proc 0's local timer results.
  static bool alwaysWriteLocal_;

  /// Whether report() should always compute global timer statistics.
  /// This requires communication equivalent to O(1) all-reduces.
  static bool writeGlobalStats_;

  //! Whether report() should report timers with zero call counts.
  static bool writeZeroTimers_;
  //@}

  /// \brief Whether setReportParameters() completed successfully.
  ///
  /// \note Keeping this helps us avoid keeping the whole
  ///   ParameterList around.
  static bool setParams_;
};


} // namespace Teuchos


namespace Teuchos {

/// \class TimeMonitorSurrogateImpl
/// \brief Implementation of TimeMonitorSurrogate that invokes TimeMonitor.
/// \warning Users should not use this class or rely on it in any way.
///   It is an implementation detail.
///
/// Please refer to the documentation of
/// TimeMonitorSurrogateImplInserter and TimeMonitorSurrogate for an
/// explanation of the purpose of this class.
class TimeMonitorSurrogateImpl : public CommandLineProcessor::TimeMonitorSurrogate
{
  virtual void summarize (std::ostream& out) {
    TimeMonitor::summarize (out);
  }
};

/// \class TimeMonitorSurrogateImplInserter
/// \brief Injects run-time dependency of a class on TimeMonitor.
/// \warning Users should not use this class or rely on it in any way.
///   It is an implementation detail.
///
/// \section Teuchos_TimeMonitorSurrogateImplInserter_Summary Summary
///
/// Classes and functions with the name "TimeMonitorSurrogate" in them
/// let CommandLineProcessor optionally call TimeMonitor::summarize(),
/// without needing to know that the TimeMonitor class exists.  This
/// allows Teuchos to put CommandLineProcessor in a separate package
/// from TimeMonitor.  We want to do this because TimeMonitor depends
/// on Comm, and is therefore in the TeuchosComm subpackage (which
/// depends on TeuchosCore), but CommandLineProcessor is in a
/// different subpackage which does not depend on Comm.
///
/// The TimeMonitorSurrogateImplInserter class' constructor ensures
/// that CommandLineProcessor gets informed about TimeMonitor even
/// before the program starts executing main().  This happens
/// automatically, without changes to main(), because we declare an
/// instance of this class in the header file.  If the TeuchosComm
/// subpackage was built and its libraries were linked in,
/// CommandLineProcessor will know about TimeMonitor.
///
/// \section Teuchos_TimeMonitorSurrogateImplInserter_Note Note to Teuchos developers
///
/// This is an instance of the
/// <a href="http://en.wikipedia.org/wiki/Dependency_injection">Dependency injection</a>
/// design pattern.  CommandLineProcessor is not supposed to know
/// about TimeMonitor, because CommandLineProcessor's subpackage does
/// not depend on TimeMonitor's subpackage.  Thus,
/// CommandLineProcessor interacts with TimeMonitor through the
/// TimeMonitorSurrogate interface.  TimeMonitorSurrogateImplInserter
/// "injects" the dependency at run time, if the TeuchosComm
/// subpackage was enabled and the application linked with its
/// libraries.
///
/// Teuchos developers could imitate the pattern of this class in
/// order to use TimeMonitor's class methods (such as summarize())
/// from any other class that does not depend on the TeuchosComm
/// subpackage.
class TimeMonitorSurrogateImplInserter {
public:
  //! Constructor: inject dependency on TimeMonitor into CommandLineProcessor.
  TimeMonitorSurrogateImplInserter () {
    if (is_null (CommandLineProcessor::getTimeMonitorSurrogate ())) {
      CommandLineProcessor::setTimeMonitorSurrogate (Teuchos::rcp (new TimeMonitorSurrogateImpl));
    }
  }
};

} // end namespace Teuchos


namespace {

// Inject the implementation in every translation unit.
Teuchos::TimeMonitorSurrogateImplInserter timeMonitorSurrogateImplInserter;

} // namespace (anonymous)

#endif // TEUCHOS_TIMEMONITOR_H
