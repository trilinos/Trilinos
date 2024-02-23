// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_MUTUALLYEXCLUSIVETIME_HPP
#define MUELU_MUTUALLYEXCLUSIVETIME_HPP

#include <string>
#include <stack>
#include "Teuchos_RCPDecl.hpp"  // for RCP

// types for explicit instantiation
#include "MueLu_FactoryBase_fwd.hpp"
#include "MueLu_Level_fwd.hpp"
#include "MueLu_BaseClass.hpp"

namespace MueLu {

/*! @class MutuallyExclusiveTime

    @brief This class wraps a Teuchos::Time and maintains a mutually exclusive property between wrapped timers.

    This class is useful to exclude from a timer the execution time of a subroutine;
    when a MutuallyExclusiveTime is running, other timers are not running.
    Timers have three states (running, stopped, paused) to enforce the mutually exclusive property.
    When the running timer is stopped, the last active timer is restarted. A stack of timers is used internally to support this functionality.
    The template parameter of this class can be used to define several sets of mutually exclusive timers.

    Note: Only one timer can be active at a time but all timers can be inactive at the same time. Timers cannot be destroyed when they are in 'paused'.

    @todo inheritence from PerformanceMonitorBase<Time> ?
*/
template <class TagName>
class MutuallyExclusiveTime : public BaseClass {
 public:
  //! @name Constructor/Destructor
  //@{
  //! Constructor
  MutuallyExclusiveTime(const std::string& name, bool startFlag = false);

  //! Destructor
  ~MutuallyExclusiveTime();
  //@}

  //! @brief Starts the timer. If a MutuallyExclusiveTime timer is running, it will be stopped.
  //! @pre Timer is not already paused.
  //! @post Timer is running. Other MutuallyExclusiveTime timers are paused or stopped.
  void start(bool reset = false);

  //! @name Functions that can only be called on the most recent timer (i.e., the running or last paused timer)
  //@{

  //! @brief Stops the timer.
  //! The previous MutuallyExclusiveTime that has been paused when this timer was started will be resumed.
  //! This method can be called on an already stopped timer or on the currently running timer.
  double stop();

  //! Pause running timer. Used internally by start().
  void pause();

  //! @brief Resume paused timer.
  //! Used internally by stop().  Timer is not reset.
  //! @pre Timer is at the top of the stack.

  void resume();

  //@}

  //! @name Query methods.
  //@{

  bool isRunning();

  bool isPaused();

  //@}

  //! Return a new MutuallyExclusiveTime that is registered with the Teuchos::TimeMonitor (for timer summary).
  // Note: this function is provided by the timer class, not by a monitor (!= Teuchos)
  static RCP<MutuallyExclusiveTime<TagName> > getNewTimer(const std::string& name);

  //! Increment the number of times this timer has been called.
  void incrementNumCalls();

  //! Name of this mutually exclusive timer.
  std::string name_;

  /*! @brief Print std::map of (child,parent) pairs for post-run analysis.

      This method print a std::map.  Each key is a child, and the value is its parent.
      The (child,parent0) pairs can be used by the MueLu script ``mueprof.sh" to print a
      nice hierarchical tree that shows both runtime flow and time for each section.
  */
  static void PrintParentChildPairs();

 private:
  //! This constructor is not public to prevent users from using Teuchos::Time::start()/stop()
  //! instead of MutuallyExclusiveTime::start()/stop(), if they have access to the underlying Teuchos::Time object.
  MutuallyExclusiveTime(RCP<Teuchos::Time> timer);

  // MutuallyExclusiveTime() { }

  //! Using an RCP allows to use Teuchos::TimeMonitor to keep track of the timer.
  RCP<Teuchos::Time> timer_;
  bool isPaused_;

  //! Stack of created timers (active or paused timers).
  // - empty when no active timer
  // - head is the active timer
  // - other timers are timers paused to enforce the mutually exclusive property of the timer set.
  static std::stack<MutuallyExclusiveTime<TagName>*> timerStack_;
  // static std::map<std::string,std::string> myParent_;

  //! Check if 'this' is the head of the stack.
  void TopOfTheStack();

  // TODO: test integrity of the stack:
  //  Head = running or paused
  //  Other timers of the stack = paused
};

extern template class MutuallyExclusiveTime<FactoryBase>;
extern template class MutuallyExclusiveTime<Level>;
extern template class MutuallyExclusiveTime<BaseClass>;

}  // namespace MueLu

#endif  // MUELU_MUTUALLYEXCLUSIVETIME_HPP
