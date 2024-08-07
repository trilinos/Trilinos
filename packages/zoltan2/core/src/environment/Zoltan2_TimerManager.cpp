// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

void TimerManager::stop(TimerType tt, const std::string &name)
{
  if (!typeSelector_[tt])
    return;

  if (stopHint_>0 && timers_[stopHint_]->name() == name){
    timers_[stopHint_]->stop();
    stopHint_--;
    return;
  }

  std::map<std::string, int>::iterator curr = timerMap_.find(name);
  if (curr != timerMap_.end()){
    timers_[curr->second]->stop();
  }
  else{  // Stopping a timer that doesn't exist.  Just create it.
    RCP<Teuchos::Time> newTimer = Teuchos::TimeMonitor::getNewTimer(name);
    newTimer->reset();     // reset to zero
    timerMap_[name] = timers_.size();
    timers_.push_back(newTimer);
    std::cerr << comm_->getRank() << ": warning, stop with no start:"
      << name.c_str() << std::endl;
  }
}

void TimerManager::start(TimerType tt, const std::string &name)
{
  if (!typeSelector_[tt])
    return;

  std::map<std::string, int>::iterator curr = timerMap_.find(name);
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
