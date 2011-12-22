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
 * \brief Scope protection wrapper for a Teuchos::Time object.
 *
 * The \c TimeMonitor class wraps a nonconst reference to a \c
 * Teuchos::Time timer object.  TimeMonitor's constructor starts the
 * timer, and its destructor stops the timer.  This ensures scope
 * safety of timers, so that no matter how a scope is exited (whether
 * the normal way or when an exception is thrown), a timer started in
 * the scope is stopped when the scope is left.
 */

/** \example TimeMonitor/cxx_main.cpp
 *
 * This is an example of how to use the \c Teuchos::TimeMonitor class.
 */

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_PerformanceMonitorBase.hpp"
#include "Teuchos_Time.hpp"

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


/// \class TimeMonitor
/// \brief A scope-safe timer wrapper class.
///
/// This class wraps a nonconst reference to a \c Teuchos::Time timer
/// object.  The TimeMonitor's constructor starts the timer, and its
/// destructor stops the timer.  This ensures scope safety of timers,
/// so that no matter how a scope is exited (whether the normal way or
/// when an exception is thrown), a timer started in the scope is
/// stopped when the scope is left.
///
/// This class also keeps track of the set of all timers as class (not
/// instance) data, and has a class method (\c summarize()) for
/// printing out global statistics (min, mean, and max over all MPI
/// processes, in an MPI build).  The \c summarize() method works
/// correctly even if some MPI processes have different timers than
/// other processes.
///
/// \warning This class must only be used to time functions that are
///   called only within the main program.  It may <i>not</i> be used in
///   pre-program setup or post-program teardown!
///
/// \note Teuchos::TimeMonitor uses the \c Teuchos::Time class internally.
///
class TEUCHOS_LIB_DLL_EXPORT TimeMonitor : public PerformanceMonitorBase<Time>
{
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

  /** \name Static functions */
  //@{

  /// \brief Return a new timer with the given name.
  ///
  /// This method wraps \c getNewCounter() (inherited from the base
  /// class) for backwards compatibiity.
  static RCP<Time> getNewTimer (const std::string& name) {
    return getNewCounter (name);
  }

  /// \brief Reset all global timers to zero.
  ///
  /// This method only affects \c Time objects created by \c
  /// getNewCounter() or \c getNewTimer().
  ///
  /// <b>Preconditions:</b><ul>
  /// <li>None of the timers must currently be running.
  /// </ul>
  static void zeroOutTimers();

  /// \brief Print summary statistics for all timers on the given communicator.
  ///
  /// The typical use case for timers is that all MPI processes on a
  /// communicator create the same set of timers, and then want to
  /// report summary statistics.  This method's default behavior
  /// (writeGlobalStats=true) is to report the mininum, arithmetic
  /// mean, and maximum for each timer.  Duplicate timers get merged
  /// additively.
  ///
  /// Note that different MPI processes may have different sets of
  /// timers.  If writeGlobalStats is true, we have to reconcile the
  /// different sets of timers somehow.  This method gives you two
  /// options: if setOp is Intersection, it computes the intersection
  /// (the common subset) of timers on all MPI processes in the
  /// communicator.  Otherwise, if setOp is Union, this method
  /// computes the union of timers on all MPI processes in the
  /// communicator.  Intersection is the default, since it means that
  /// all reported timers exist on all participating processes.
  ///
  /// Suppose there are \f$P\f$ MPI processes in the communicator and
  /// \f$N\f$ unique timers in the global union.  This method requires
  /// \f$O(\log P)\f$ messages (\f$O(1)\f$ "reductions" and exactly 1
  /// "broadcast") and \f$O(N)\f$ per-processor storage (in the worst
  /// case) when computing either the intersection or the union of
  /// timers (the algorithm is similar in either case).  The whole
  /// algorithm takes at worst \f$O(N (\log N) (\log P))\f$ time along
  /// the critical path (i.e., on the "slowest MPI process").
  ///
  /// \param comm [in] Communicator whose process(es) will participate
  ///   in the gathering of timer statistics.
  ///
  /// \param out [out] Output stream to which to write.  This will
  ///   only be used on the process with Rank 0 in the communicator.
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
  ///   display the min, average (arithmetic mean), and max of all
  ///   timings over all processes in the communicator.  If there is
  ///   only one MPI process or if this is a non-MPI build of
  ///   Trilinos, we only show the "global" timings, without the
  ///   "statistics" that would be all the same anyway.
  ///
  /// \param writeZeroTimers [in] If false, do not display results for
  ///   timers that have never been called (numCalls() == 0).  If
  ///   true, display results for all timers.
  ///
  /// \param setOp [in] If Intersection, compute and display the
  ///   intersection of all created timers over all processes in the
  ///   communicator.  If Union, compute and display the union of all
  ///   created timers over all processes in the communicator.
  ///
  /// \note If writeGlobalStats is true, this method <i>must</i> be
  ///   called by all processes in the communicator.  This method will
  ///   <i>only</i> perform communication if writeGlobalStats is true.
  static void 
  summarize (const RCP<const Comm<int> >& comm,
             std::ostream &out=std::cout, 
	     const bool alwaysWriteLocal=false,
	     const bool writeGlobalStats=true,
	     const bool writeZeroTimers=true,
	     const ECounterSetOp setOp=Intersection);

  /// \brief Print summary statistics for all timers on all (MPI) processes.
  ///
  /// This is an overload of the above \c summarize() method for when
  /// the caller does not want to provide a communicator explicitly.
  /// This method "does the right thing" in that case.  Specifically:
  /// - If Trilinos was not built with MPI support, this method
  ///   assumes a serial "communicator" containing one process.  
  /// - If Trilinos was built with MPI support and MPI has been
  ///   initialized (via \c MPI_Init() or one of the wrappers in
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
  ///   <i>must</i> call it on all processes in MPI_COMM_WORLD.
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
	     const ECounterSetOp setOp=Intersection);
  //@}

};


} // namespace Teuchos


#endif // TEUCHOS_TIMEMONITOR_H
