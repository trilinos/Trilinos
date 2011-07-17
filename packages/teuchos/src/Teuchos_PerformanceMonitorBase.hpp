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
#include "Teuchos_MPIComm.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_PerformanceMonitorUtils.hpp"
#include "Teuchos_TableFormat.hpp"

namespace Teuchos
{
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

  /** \brief Get the format that will be used to print a summary of
   * results.
   */
  static TableFormat& format()
    {
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
  static RCP<T> 
  lookupCounter (const std::string& name)
  {
    typename Array<RCP<T> >::const_iterator it = std::find (counters().begin(), counters.end(), name);
    if (it == counters().end())
      return null;
    else
      return *it;
  }

protected:
    
  /** \brief Access to the counter. */
  const T& counter() const { return counter_; }
    
  /** \brief Access to the counter. */
  T& counter() { return counter_; }

  /** \brief Indicate whether the current call is recursive.
   *
   * This can matter in cases such as timing where we don't want to start and
   * stop timers multiple times within a single call stack.
   */
  bool isRecursiveCall() const { return isRecursiveCall_; }
      
  //! The array of counters which were created using \c getNewCounter().
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
    
  T& counter_;
    
  bool isRecursiveCall_;
    
};

  
} // namespace Teuchos


#endif
