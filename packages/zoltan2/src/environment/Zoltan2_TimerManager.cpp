// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
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
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

/*! \file Zoltan2_TimerManager.cpp
 *  \brief Definitions for TimerManager.
 */

#include <Zoltan2_TimerManager.hpp>

namespace Zoltan2{

TimerManager::TimerManager(const RCP<const Comm<int> > &comm, 
  std::ofstream *os, TimerType tt):
  comm_(comm), myOS_(NULL), fileOS_(os), ttype_(tt), 
  typeSelector_(), timers_(), timerMap_(), stopHint_(-1) 
{
  if (fileOS_ == NULL)
    ttype_ = NO_TIMERS;

  typeSelector_.reset();
  typeSelector_.set(ttype_);

  if (ttype_ == BOTH_TIMERS){
    typeSelector_.set(MACRO_TIMERS);
    typeSelector_.set(MICRO_TIMERS);
  }
}

TimerManager::TimerManager(const RCP<const Comm<int> > &comm, 
  std::ostream *os, TimerType tt):
  comm_(comm), myOS_(os), fileOS_(NULL), ttype_(tt),
  typeSelector_(), timers_(), timerMap_(), stopHint_(-1) 
{
  if (myOS_ == NULL)
    ttype_ = NO_TIMERS;

  typeSelector_.reset();
  typeSelector_.set(ttype_);

  if (ttype_ == BOTH_TIMERS){
    typeSelector_.set(MACRO_TIMERS);
    typeSelector_.set(MICRO_TIMERS);
  }
}

TimerManager::~TimerManager()
{
  if (fileOS_ != NULL){
    fileOS_->close();
    fileOS_=NULL;
  }
}

void TimerManager::stop(TimerType tt, const string &name)
{
  if (!typeSelector_[tt])
    return;

  if (stopHint_>0 && timers_[stopHint_]->name() == name){
    timers_[stopHint_]->stop();
    stopHint_--;
    return;
  }

  std::map<string, int>::iterator curr = timerMap_.find(name);
  if (curr != timerMap_.end()){
    timers_[curr->second]->stop();
  }
  else{  // Stopping a timer that doesn't exist.  Just create it.
    RCP<Teuchos::Time> newTimer = Teuchos::TimeMonitor::getNewTimer(name);
    newTimer->reset();     // reset to zero
    timerMap_[name] = timers_.size();
    timers_.push_back(newTimer);
    cerr << comm_->getRank() << ": warning, stop with no start" << endl;
  }
}

void TimerManager::start(TimerType tt, const string &name)
{
  if (!typeSelector_[tt])
    return;

  std::map<string, int>::iterator curr = timerMap_.find(name);
  int index = -1;
  if (curr == timerMap_.end()){
    RCP<Teuchos::Time> newTimer = Teuchos::TimeMonitor::getNewTimer(name);
    index = timers_.size();
    timerMap_[name] = index;
    timers_.push_back(newTimer);
  }
  else{
    index = curr->second;
  }

  timers_[index]->start();
  timers_[index]->incrementNumCalls();
  stopHint_ = index;
}

void TimerManager::print() const
{
  if (fileOS_)
    Teuchos::TimeMonitor::summarize(comm_.ptr(), *fileOS_);
  else if (myOS_)
    Teuchos::TimeMonitor::summarize(comm_.ptr(), *myOS_);
}

void TimerManager::printAndResetToZero() 
{
  print();
  Teuchos::TimeMonitor::zeroOutTimers();
  if (fileOS_){
    fileOS_->close();
    fileOS_ = NULL;
  }
}

}  // namespace Zoltan2
