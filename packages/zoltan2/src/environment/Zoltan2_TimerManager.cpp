// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_TimerManager.cpp
 *  \brief Definitions for TimerManager.
 */

#include <Zoltan2_TimerManager.hpp>

namespace Zoltan2{

TimerManager::TimerManager(const RCP<const Comm<int> > &comm, 
  std::ofstream *os):
  comm_(comm), myOS_(NULL), fileOS_(os), timers(), timerMap(), stopHint(-1) {}

TimerManager::TimerManager(const RCP<const Comm<int> > &comm, 
  std::ostream *os):
  comm_(comm), myOS_(os), fileOS_(NULL), timers(), timerMap(), stopHint(-1) {}

TimerManager::~TimerManager()
{
  if (fileOS_ != NULL){
    fileOS_->close();
    fileOS_=NULL;
  }
}

void TimerManager::stop(const string &name)
{
  if (stopHint>0 && timers[stopHint]->name() == name){
    timers[stopHint]->stop();
    stopHint--;
    return;
  }

  std::map<string, int>::iterator curr = timerMap.find(name);
  if (curr != timerMap.end()){
    timers[curr->second]->stop();
  }
  else{  // Stopping a timer that doesn't exist.  Just create it.
    RCP<Teuchos::Time> newTimer = Teuchos::TimeMonitor::getNewTimer(name);
    newTimer->reset();     // reset to zero
    timerMap[name] = timers.size();
    timers.push_back(newTimer);
    cerr << comm_->getRank() << ": warning, stop with no start" << endl;
  }
}

void TimerManager::start(const string &name)
{
  std::map<string, int>::iterator curr = timerMap.find(name);
  int index = -1;
  if (curr == timerMap.end()){
    RCP<Teuchos::Time> newTimer = Teuchos::TimeMonitor::getNewTimer(name);
    index = timers.size();
    timerMap[name] = index;
    timers.push_back(newTimer);
  }
  else{
    index = curr->second;
  }

  timers[index]->start();
  timers[index]->incrementNumCalls();
  stopHint = index;
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
