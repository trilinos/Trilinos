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

#ifndef TEUCHOS_PERFORMANCEMONITORBASE_H
#define TEUCHOS_PERFORMANCEMONITORBASE_H

/*! \file Teuchos_PerformanceMonitorBase.hpp
  \brief Provides common capabilities for collecting and reporting
  performance data across processors
*/

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_TableFormat.hpp"

namespace Teuchos
{
  /// \brief Set operation type for \c mergeCounterNames() to perform.
  ///
  /// The \c mergeCounterNames() function merges sets of counter names
  /// over all MPI processes in a communicator.  Different MPI
  /// processes may have created different sets of counters.  This
  /// enum allows the caller to specify how mergeCounterNames() picks
  /// the global set of timers.
  enum ECounterSetOp { Intersection, Union };

  /// \brief Merge counter names over all processors.
  ///
  /// Different MPI processes may have created different sets of
  /// counters.  Use this function to reconcile the sets among
  /// processes, either by computing their intersection or their
  /// union.  This is done using a reduction to MPI Rank 0 (relative
  /// to the given communicator) and a broadcast to all processes
  /// participating in the communicator.  We use a
  /// reduce-and-broadcast rather than just a reduction, so that all
  /// participating processes can use the resulting list of global
  /// names as lookup keys for computing global statistics.
  ///
  /// \param comm [in] Communicator over which to merge.
  ///
  /// \param localNames [in] The calling MPI process' list of (local)
  ///   counter names.
  ///
  /// \param globalNames [out] On output, on each MPI process: the
  ///   results of merging the counter names.
  ///
  /// \param setOp [in] If Intersection, globalNames on output
  ///   contains the intersection of all sets of counter names.  If
  ///   Union, globalNames on output contains the union of all sets of
  ///   counter names.
  void
  mergeCounterNames (const Comm<int>& comm,
                     const Array<std::string>& localNames,
                     Array<std::string>& globalNames,
                     const ECounterSetOp setOp);

  /// \class PerformanceMonitorBase
  ///
  /// \brief Common capabilities for collecting and reporting
  ///   performance data across processors.
  ///
  /// PerformanceMonitorBase is templated on a counter type T (which
  /// might be a timer or a flop counter). The common capability of the
  /// counter type is a counter for the number of calls.  Derived counter
  /// types may supply additional features.
  ///
  /// PerformanceMonitorBase's constructor increments its counter's
  /// call count.  Subclasses of PerformanceMonitorBase may do more
  /// upon construction or destruction; for example, TimeMonitor
  /// starts its timer on construction and stops it on destruction.
  ///
  /// This class keeps a static list of all counters created using the
  /// getNewCounter() method during the course of a run. Counts from
  /// this list can then be printed out at the end of the run.
  /// Subclasses of PerformanceMonitorBase, such as TimeMonitor, may
  /// use this list to do things like compute global timer statistics
  /// over all the MPI processes.
  ///
  /// PerformanceMonitorBase requires that the counter type T provide at
  /// least the following interface:
  ///
  /// \code
  /// // Constructor taking an std::string argument (the counter name).
  /// T (const std::string&);
  ///
  /// // Return the name of the counter.
  /// const std::string& name () const;
  ///
  /// // Add one to the number of calls (the number of times the counter
  /// // was started).
  /// void incrementNumCalls ();
  ///
  /// // Return the number of calls (see incrementNumCalls () above).
  /// int numCalls () const;
  ///
  /// // Indicate whether the counter is already running.
  /// bool isRunning () const;
  /// \endcode
  ///
  template <class T>
  class PerformanceMonitorBase
  {
  public:
    /** \brief Construct with a counter. */
    PerformanceMonitorBase(T& counter_in, bool reset=false)
      : counter_(counter_in), isRecursiveCall_(counter_.isRunning())
    {
      (void) reset;  // get rid of "unused parameter" warning
      counter_.incrementNumCalls ();
    }

    /// \brief Destructor.
    ///
    /// The destructor for the base class does nothing.  We provide a
    /// virtual destructor for memory safety of derived classes.
    virtual ~PerformanceMonitorBase() {}

    /** \brief Create a new counter with the specified name and add it
     *   to a global set of counters of this type.
     *
     * If the counter already exists, just return the existing
     * counter.  If the counter doesn't already exist, create a new
     * counter with that name and return it.
     *
     * New counters should usually be created in this way rather than
     * through a direct constructor call.  This lets
     * PerformanceMonitorBase keep track of them, so that methods like
     * summarize() and report() know about them.  Timers created in
     * other ways are not included in the reports printed by these
     * methods.
     */
    static RCP<T> getNewCounter (const std::string& name);

    //! Table format that will be used to print a summary of timer results.
    static TableFormat& format()
    {
      // WARNING This is not thread safe!  In particular, if multiple
      // threads call this method for the "first time" at the same
      // time, the RCP may be initialized incorrectly.
      static RCP<TableFormat> rtn=rcp(new TableFormat());
      return *rtn;
    }

    /// \brief Return the first counter with the given name, or null if none.
    ///
    /// It is currently possible to create multiple counters with the
    /// same name using \c getNewCounter().  If multiple counters with
    /// the given name exist, this method simply returns the first in
    /// the list.  Do not rely on the ability to create multiple
    /// counters with the same name; this may go away in the future.
    static RCP<T>
    lookupCounter (const std::string& name);

    /// \brief "Forget" about all counters created with getNewCounter().
    ///
    /// This removes all counters from the current set of counters (as
    /// would be returned by counters()).
    static void clearCounters ();

    /// \brief "Forget" about all counters created with getNewCounter().
    ///
    /// This removes all counters from the current list of counters
    /// (as would be returned by counters()).
    ///
    /// \warning This method is DEPRECATED, because the name is
    ///   inaccurate (the template parameter of PerformanceMonitorBase
    ///   may be any kind of performance counter, not just a timer).
    ///   Use clearCounters() instead.
    static TEUCHOS_DEPRECATED void clearTimers ();

    /// \brief "Forget" about any counters with the given name.
    ///
    /// If one or more counters with the given name was created using
    /// getNewCounter(), calling this method with that name will
    /// remove them from the global list of counters.
    static void clearCounter (const std::string& name);

    /// \brief "Forget" about any counters with the given name.
    ///
    /// If one or more counters with the given name was created using
    /// getNewCounter(), calling this method with that name will
    /// remove them from the global list of counters.
    ///
    /// \warning This method is DEPRECATED, because the name is
    ///   inaccurate (the template parameter of PerformanceMonitorBase
    ///   may be any kind of performance counter, not just a timer).
    ///   Use clearCounter() instead.
    static TEUCHOS_DEPRECATED void clearTimer (const std::string& name);

  protected:

    //! Constant access to the instance's counter reference.
    const T& counter() const { return counter_; }

    //! Nonconstant access to the instance's counter reference.
    T& counter() { return counter_; }

    /// \brief Whether we are currently in a recursive call of the counter.
    ///
    /// Subclasses of PerformanceMonitorBase may use this information
    /// to control whether to start or stop the given counter.  This
    /// matters in cases such as timing, where we don't want to start
    /// and stop timers multiple times within a single call stack.
    bool isRecursiveCall() const { return isRecursiveCall_; }

    //! The class' array of all counters that were created with getNewCounter().
    static std::map<std::string, RCP<T> >& counters()
    {
      // Use the "Meyers Trick" to create static data safely.
      //
      // WARNING This is not thread safe!  In particular, if multiple
      // threads call counters() for the "first time" at the same
      // time, the array may be initialized incorrectly.
      static std::map<std::string, RCP<T> > theCounters;
      return theCounters;
    }

  private:

    //! Reference to the counter being wrapped.
    T& counter_;

    //! Whether we are currently in a recursive call of the counter.
    bool isRecursiveCall_;
  };

  template<class T>
  RCP<T>
  PerformanceMonitorBase<T>::getNewCounter (const std::string& name)
  {
    typedef std::map<std::string, RCP<T> > map_type;
    typedef typename map_type::iterator iter_type;

    map_type& ctrs = counters ();
    iter_type it = ctrs.find (name);
    if (it == ctrs.end ()) {
      RCP<T> newCounter = rcp (new T (name));
      ctrs.insert (it, std::make_pair (name, newCounter));
      return newCounter;
    } else {
      return it->second;
    }
  }

  template<class T>
  RCP<T>
  PerformanceMonitorBase<T>::lookupCounter (const std::string& name)
  {
    typedef std::map<std::string, RCP<T> > map_type;
    typedef typename map_type::iterator iter_type;

    map_type& ctrs = counters ();
    iter_type it = ctrs.find (name);
    if (it == ctrs.end ()) {
      return null;
    } else {
      return it->second;
    }
  }

  template<class T>
  void
  PerformanceMonitorBase<T>::clearCounter (const std::string& name)
  {
    counters ().erase (name);
  }

  template<class T>
  void
  PerformanceMonitorBase<T>::clearTimer (const std::string& name)
  {
    clearCounter (name);
  }

  template<class T>
  void
  PerformanceMonitorBase<T>::clearTimers ()
  {
    clearCounters ();
  }

  template<class T>
  void
  PerformanceMonitorBase<T>::clearCounters ()
  {
    counters ().clear ();
  }

} // namespace Teuchos

#endif // TEUCHOS_PERFORMANCEMONITORBASE_H
