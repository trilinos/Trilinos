// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_PERFORMANCEMONITORBASE_H
#define TEUCHOS_PERFORMANCEMONITORBASE_H

/// \file Teuchos_PerformanceMonitorBase.hpp
/// \brief Common capabilities for collecting and reporting
///   performance data collectively across MPI processes.

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_TableFormat.hpp"
#include <cstdlib> // atexit

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

  /**
   * merge for unsorted lists.  New entries are at the bottom of the list
   * @param localNames - The calling MPI process' list of (local)
   * counter names.
   * @param globalNames - Global list of names
   * @param setOp If Intersection, globalNames on output
   *   contains the intersection of all sets of counter names.  If
   *   Union, globalNames on output contains the union of all sets of
   *   counter names.
   */
  void
  unsortedMergePair(const Array<std::string>& localNames,
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
    //! Construct with a counter.
    PerformanceMonitorBase(T& counter_in, bool reset=false)
      : counter_(counter_in), isRecursiveCall_(counter_.isRunning())
    {
      (void) reset;  // get rid of "unused parameter" warning
      counter_.incrementNumCalls ();
    }

    //! Default constructor is deleted, since it would be unsafe.
    PerformanceMonitorBase () = delete;

    /// \brief Destructor.
    ///
    /// The destructor for the base class does nothing.  We provide a
    /// virtual destructor for memory safety of derived classes.
    virtual ~PerformanceMonitorBase() = default;

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

  private:
    /// \brief Free the singleton returned by format().
    ///
    /// \warning Only for use as atexit() handler.
    ///
    /// \warning This method is not reentrant.  In particular, if
    ///   multiple threads call this method at the same time, they
    ///   might manage to double-delete format_.  This could only
    ///   happen if the format() method below is called twice by
    ///   different threads.
    static void freeTableFormat () {
      if (format_ != nullptr) {
        delete format_;
        format_ = nullptr;
      }
    }

    /// \brief Free the singleton returned by counters().
    ///
    /// \warning Only for use as atexit() handler.
    ///
    /// \warning This method is not reentrant.  In particular, if
    ///   multiple threads call this method at the same time, they
    ///   might manage to double-delete counters_.  This could only
    ///   happen if the counters() method below is called twice by
    ///   different threads.
    static void freeCounters () {
      if (counters_ != nullptr) {
        delete counters_;
        counters_ = nullptr;
      }
    }

  public:
    /// \brief Table format that will be used to print a summary of
    ///   timer results.
    ///
    /// \warning This method is not reentrant.  In particular, if
    ///   multiple threads call this method at the same time, they
    ///   might manage to double-register the atexit() handler for
    ///   format_.  This could only happen if this method is called
    ///   twice by different threads.
    static TableFormat& format ()
    {
      if (format_ == nullptr) {
        format_ = new TableFormat ();
        // It _is_ possible for atexit() to fail (e.g., because it ran
        // out of memory for storing callbacks).  We could throw an
        // exception here in that case, but I think it's better just
        // to let the minor memory leak happen.
        static_cast<void>( atexit(freeTableFormat) );
      }
      TEUCHOS_TEST_FOR_EXCEPTION(
        format_ == nullptr, std::logic_error, "Teuchos::PerformanceMonitorBase::"
        "format: Should never get here!  format_ is nullptr.");

      return *format_;
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

    /// \brief "Forget" about any counters with the given name.
    ///
    /// If one or more counters with the given name was created using
    /// getNewCounter(), calling this method with that name will
    /// remove them from the global list of counters.
    static void clearCounter (const std::string& name);

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

    /// \brief Array of all counters that were created with
    ///   getNewCounter() on the calling (MPI) process.
    ///
    /// \warning This method is not reentrant.
    static std::map<std::string, RCP<T> >& counters ()
    {
      if (counters_ == nullptr) {
        counters_ = new std::map<std::string, RCP<T> > ();
        // It _is_ possible for atexit() to fail (e.g., because it ran
        // out of memory for storing callbacks).  We could throw an
        // exception here in that case, but I think it's better just
        // to let the minor memory leak happen.
        static_cast<void>( atexit(freeCounters) );
      }
      TEUCHOS_TEST_FOR_EXCEPTION(
        counters_ == nullptr, std::logic_error, "Teuchos::PerformanceMonitorBase::"
        "counters: Should never get here!  counters_ is nullptr.");

      return *counters_;
    }

  private:
    //! Singleton object returned by format().
    static TableFormat* format_;

    //! Singleton object returned by counters().
    static std::map<std::string, RCP<T> >* counters_;

    //! Reference to the counter being wrapped.
    T& counter_;

    //! Whether we are currently in a recursive call of the counter.
    bool isRecursiveCall_;
  };

  template<class T>
  TableFormat*
  PerformanceMonitorBase<T>::format_ = nullptr;

  template<class T>
  std::map<std::string, RCP<T> >*
  PerformanceMonitorBase<T>::counters_ = nullptr;

  template<class T>
  RCP<T>
  PerformanceMonitorBase<T>::getNewCounter (const std::string& name)
  {
    typedef std::map<std::string, RCP<T> > map_type;
    typedef typename map_type::iterator iter_type;

    map_type& ctrs = counters ();
    iter_type it = ctrs.find (name);
    RCP<T> newCounter = null;
    if (it == ctrs.end ()) {
      newCounter = rcp (new T (name));
#ifdef HAVE_TEUCHOS_DEBUG
      const bool wasNotThere = ctrs.insert (std::make_pair (name, newCounter)).second;
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! wasNotThere, std::logic_error,
        "getNewCounter: insert() claims that timer \"" << name << "\" was "
        "already there in the map, even though find() claims that it was not.  "
        "Please report this bug to the Teuchos developers.");
#else
      // Use the returned iterator to optimize insertion.
      ctrs.insert (it, std::make_pair (name, newCounter));
#endif // HAVE_TEUCHOS_DEBUG
    } else {
      newCounter = it->second;
#ifdef HAVE_TEUCHOS_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION(
        it->second.is_null (), std::logic_error,
        "getNewCounter: Timer \"" << name << "\" was already there in the map, "
        "but looking it up by name resulted in a null timer.  "
        "Please report this bug to the Teuchos developers.");
      TEUCHOS_TEST_FOR_EXCEPTION(
        name != it->second->name (), std::logic_error,
        "getNewCounter: Timer \"" << name << "\" was already there in the map, "
        "but looking it up by name resulted in a timer with a different name \""
        << it->second->name () << "\".  Please report this bug to the Teuchos "
        "developers.");
#endif // HAVE_TEUCHOS_DEBUG
    }

#ifdef HAVE_TEUCHOS_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION(
      newCounter.is_null (), std::logic_error,
      "getNewCounter: At end of method, when creating timer \"" << name
      << "\", newCounter is null.  Please report this bug to the Teuchos "
      "developers.");
#endif // HAVE_TEUCHOS_DEBUG
    return newCounter;
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
  PerformanceMonitorBase<T>::clearCounters ()
  {
    counters ().clear ();
  }

} // namespace Teuchos

#endif // TEUCHOS_PERFORMANCEMONITORBASE_H
