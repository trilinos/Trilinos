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
#include "Teuchos_Comm.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
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
  /// (Otherwise, different processes
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
  /// A PerformanceMonitorBase will increment its call counter upon every
  /// constructor call. Derived types might do more upon construction or
  /// destruction; for example, a timer will start upon construction and
  /// stop upon destruction.
  ///
  /// The class keeps a static list of all counters created using
  /// the \c getNewCounter() method during the course of a run. Counts
  /// from this list can then be printed out at the end of the run. 
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
      (void)reset;  // get rid of "unused parameter" warning
      counter_.incrementNumCalls();
    }

    /// \brief Destructor.
    ///
    /// The destructor for the base class does nothing.  We provide a
    /// virtual destructor for memory safety of derived classes.
    virtual ~PerformanceMonitorBase() {}
    
    /** \brief Create a new counter with the specified name and append it to a
     * global list of counters of this type.
     *
     * New counters should usually be created in this way rather than
     * through a direct constructor call.  This lets
     * PerformanceMonitorBase keep track of them, for example, so that
     * you can access the list of timers by calling \c counters().
     * Timers created in other ways are not included in the list
     * returned by \c counters().
     */
    static RCP<T> getNewCounter(const std::string& name)
    {
      RCP<T> rtn = rcp(new T(name), true);
      counters().append(rtn);
      return rtn;
    }

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
    /// the list.
    ///
    /// If you want to create a counter with a given name if it doesn't
    /// yet exist, do the following:
    /// \code
    /// RCP<T> counter = PerformanceMonitorBase<T>::lookupCounter (name);
    /// if (counter.is_null())
    ///   counter = PerformanceMonitorBase<T>::getNewCounter (name);
    /// \endcode
    /// We do not offer this functionality as a single class method
    /// ("lookupOrCreateCounter()"), because it is easy to misspell a
    /// timer name.  If timers that don't exist yet are created
    /// silently, misspelled timer names can cause hard-to-find bugs.
    static RCP<T> 
    lookupCounter (const std::string& name);

    /// \brief "Forget" about all timers created with \c getNewCounter().
    ///
    /// This removes all timers from the current list of timers (as
    /// would be returned by \c counters()).
    static void clearTimers ();

    //! "Forget" about any timers named "name" created with \c getNewCounter().
    static void clearTimer (const std::string& name);

  protected:
    
    //! Constant access to the counter reference.
    const T& counter() const { return counter_; }
    
    //! Nonconstant access to the counter reference.
    T& counter() { return counter_; }

    /// \brief Whether we are currently in a recursive call of the counter.
    ///
    /// This can matter in cases such as timing, where we don't want to
    /// start and stop timers multiple times within a single call stack.
    bool isRecursiveCall() const { return isRecursiveCall_; }
      
    //! The array of all counters created using \c getNewCounter().
    static Array<RCP<T> >& counters() 
    {
      // Use the "Meyers Trick" to create static data safely.
      //
      // WARNING This is not thread safe!  In particular, if multiple
      // threads call counters() for the "first time" at the same
      // time, the array may be initialized incorrectly.
      static Array<RCP<T> > rtn;
      return rtn;
    }

  private:

    //! Reference to the counter being wrapped.
    T& counter_;
    
    //! Whether we are currently in a recursive call of the counter.
    bool isRecursiveCall_;
  };


  template<class T>
  RCP<T> 
  PerformanceMonitorBase<T>::lookupCounter (const std::string& name)
  {
    // If we change how we store the counters, the find() as well as
    // the iterator type will have to change.
    typename Array<RCP<T> >::const_iterator it = 
      std::find (counters().begin(), counters.end(), name);
    if (it == counters().end())
      return null;
    else
      return *it;
  }

  template<class T>
  void
  PerformanceMonitorBase<T>::clearTimer (const std::string& name)
  {
    Array<RCP<T> > newCounters;
    // Only fill newCounters with counters whose name is not name.
    //
    // Alas, standard STL lacks compose1.
    // std::remove_copy_if (counters().begin(), 
    // 			 counters().end(), 
    // 			 std::back_inserter (newCounters),
    // 			 std::compose1 (std::bind2nd (std::equal_to<std::string>, name),
    // 					std::mem_fun (&T::name)));
    typedef typename Array<RCP<T> >::const_iterator iter_t;
    for (iter_t it = counters().begin(); it != counters().end(); ++it)
      {
	if (it->name() != name)
	  newCounters.push_back (*it);
      }
    counters().swap (newCounters);
  }

  template<class T>
  void
  PerformanceMonitorBase<T>::clearTimers ()
  {
    // Just resizing an Array to have length zero may not necessarily
    // free its storage.  The standard idiom is to swap with an empty
    // array.
    Array<RCP<T> > newCounters;
    counters().swap (newCounters);
  }




  
} // namespace Teuchos


#endif
