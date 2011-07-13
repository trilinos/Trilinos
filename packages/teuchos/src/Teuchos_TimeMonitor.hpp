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
 * \brief Timer class that starts when constructed and stops when the
 * destructor is called
 */

/** \example TimeMonitor/cxx_main.cpp
 *
 * This is an example of how to use the Teuchos::TimeMonitor class.
 */

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_PerformanceMonitorBase.hpp"
#include "Teuchos_Time.hpp"

/** \brief Defines a static non-member function that returns a time monitor.
 */ 
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
 * Note that the name of the timer can be formated with stream inserts.
 * For example, we can define a time monitor for a function as follows:
 
 \code

 template<typename Scalar>
 void foo()
 {
 TEUCHOS_FUNC_TIME_MONITOR(
 "foo<"<<Teuchos::ScalarTraits<Scalar>::name()<<">()"
 );
 ...
 }

 \endcode

 * The timer can then be printed at the end of the program using

 \code

 Teuchos::TimeMonitor::summarize(std::cout);

 \endcode
 
*/
#define TEUCHOS_FUNC_TIME_MONITOR( FUNCNAME ) \
  TEUCHOS_FUNC_TIME_MONITOR_DIFF( FUNCNAME, main )


namespace Teuchos {


/** \brief A scope-safe timer wrapper class.
 *
 * TimeMonitor objects start the timer when constructed, and stop the
 * timer when the destructor is called.  Termination upon destruction
 * lets this timer behave correctly even if scope is exited because of
 * an exception.  TimeMonitor also keeps track of the set of all
 * timers, and has a method (\c summarize()) for printing out global
 * statistics (min, mean, and max over all MPI processes, in an MPI
 * build).
 *
 * \warning This class must only be used to time functions that are
 *   called only within the main program.  It may <i>not</i> be used in
 *   pre-program setup or post-program teardown!
 *
 * \note Teuchos::TimeMonitor uses the \c Teuchos::Time class internally.
 */
class TEUCHOS_LIB_DLL_EXPORT TimeMonitor : public PerformanceMonitorBase<Time>
{
public:

  /** \name Constructor/Destructor */
  //@{
 
  /// \brief Constructor starts the timer.
  ///
  /// \param timer [in/out] Reference to the timer to be wrapped.
  ///
  /// \param reset [in] If true, reset the timer before starting it.
  ///   Default behavior is not to reset the timer.
  TimeMonitor (Time& timer, bool reset=false);
  
  //! Destructor causes timer to stop.
  ~TimeMonitor();
  //@}

  /** \name Static functions */
  //@{
 
  /** \brief Wrapping of getNewCounter() for backwards compatibiity with old
   * code.
   */
  static Teuchos::RCP<Time> getNewTimer (const std::string& name) {
    return getNewCounter(name);
  }

  /** \brief Reset the global timers to zero.
   *
   * <b>Preconditions:</b><ul>
   * <li>None of the timers must currently be running!
   * </ul>
   */
  static void zeroOutTimers();
 
  /// \brief Print summary statistics for all timers. 
  ///
  /// The typical use case for timers is that all MPI processes create
  /// the same set of timers, and then want to report summary
  /// statistics.  This method's default behavior
  /// (writeGlobalStats=true) is to report the mininum, arithmetic
  /// mean, and maximum for each timer.  Duplicate timers get merged
  /// additively.
  ///
  /// Note that different MPI processes may have different sets of
  /// timers.  If writeGlobalStats is true, we have to reconcile the
  /// different sets of timers somehow.  This method gives you two
  /// options: if globalUnionOfTimers is true, it computes the
  /// intersection (the common subset) of timers on all MPI processes,
  /// otherwise it computes the union of timers on all MPI processes.
  ///
  /// Suppose there are \f$P\f$ MPI processes, \f$N\f$ unique timers
  /// in the global union, and \f$n\f$ unique timers in the global
  /// intersection.  This method requires \f$O(\log P)\f$ messages (a
  /// single "reduction") and \f$O(N)\f$ per-processor storage (in the
  /// worst case) when computing either the intersection or the union
  /// of timers (the algorithm is similar in either case).  The whole
  /// algorithm takes at worst \f$O(N (\log N) (\log P))\f$ time along
  /// the critical path.
  ///
  /// \param out [out] Output stream to which to write.  This will
  ///   only be used on MPI Rank 0.
  ///
  /// \param alwaysWriteLocal [in] If true, MPI Proc 0 will write its
  ///   local timings to the given output stream.  Defaults to false,
  ///   since the global statistics are more meaningful.  We also
  ///   exclude this case if the local set of timers differs from the
  ///   global set of timers (either the union or the intersection,
  ///   depending on \c globalUnionOfTimers), since that would mess up
  ///   the table's formatting.
  ///
  /// \param writeGlobalStats [in] If true (the default), compute and
  ///   display the min, average (arithmetic mean), and max of all
  ///   timings over all processors (in MPI_COMM_WORLD).  If there is
  ///   only one MPI process or if this is a non-MPI build of
  ///   Trilinos, we only show the "global" timings, without the
  ///   "statistics" that would be all the same anyway.
  ///
  /// \param writeZeroTimers [in] If false, do not display results for
  ///   timers that have never been called (numCalls() == 0).  If
  ///   true, display results for all timers.
  ///
  /// \param globalUnionOfTimers [in] If true, compute and display the
  ///   union of all created timers over all processors.  If false
  ///   (the default), compute and display the intersection of all
  ///   created timers over all processors.
  ///
  /// \note If writeGlobalStats is true, this method <i>must</i> be
  ///   called by all processors.
  static void 
  summarize (std::ostream &out=std::cout, 
	     const bool alwaysWriteLocal=false,
	     const bool writeGlobalStats=true,
	     const bool writeZeroTimers=true,
	     const bool globalUnionOfTimers=false);
  //@}

 private:

  typedef std::pair<std::string, std::pair<double, int> > timer_datum_t;

  //! Collect and sort local timer data by timer names.
  static void collectLocalTimerData (Array<timer_datum_t>& localData);

  //! Locally filter out timer data with zero call counts.
  static void filterZeroData (Array<timer_datum_t>& timerData);

  /// \brief Merge timer data over all processors.
  ///
  /// \param comm [in] Communicator over which to merge.
  /// \param localTimerData [in] Each processor's timer data.
  /// \param globalTimerData [out] On output, on MPI Proc 0: the
  ///   results of merging the timer data.
  /// \param intersect [in] If true, globalTimerData on output
  ///   contains the intersection of all timers.  If false,
  ///   globalTimerData on output contains the union of all timers.
  static void
  mergeTimers (const Comm<int>& comm, 
	       const Array<timer_datum_t>& localTimerData,
	       Array<timer_datum_t>& globalTimerData,
	       const bool intersect);

  //! Recursive helper function for \c mergeTimers().
  static void 
  mergeTimersHelper (const Comm<int>& comm, 
		     const int myRank,
		     const int left,
		     const int right, // inclusive range [left, right]
		     const Array<timer_datum_t>& localTimerData,
		     Array<timer_datum_t>& globalTimerData,
		     const bool intersect);

  //! Helper function for \c mergeTimersHelper().
  static void
  mergeTimersPair (const Comm<int>& comm, 
		   const int myRank,
		   const int left,
		   const int mid,
		   const Array<timer_datum_t>& localTimerData,
		   Array<timer_datum_t>& globalTimerData,
		   const bool intersect);




};


} // namespace Teuchos


#endif // TEUCHOS_TIMEMONITOR_H
