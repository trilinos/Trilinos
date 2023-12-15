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

#include <map>
#include <iostream>                      // for basic_ostream, etc
#include <utility>                       // for pair
#include "Teuchos_FancyOStream.hpp"      // for basic_FancyOStream, etc
#include "Teuchos_RCP.hpp"               // for RCP::operator->, etc
#include "Teuchos_TestForException.hpp"  // for TEUCHOS_TEST_FOR_EXCEPTION
#include "Teuchos_Time.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_BaseClass.hpp"
#include "MueLu_VerbosityLevel.hpp"  // for MsgType::Debug, etc
#include "MueLu_MutuallyExclusiveTime.hpp"
#include "MueLu_FactoryBase.hpp"
#include "MueLu_Level.hpp"

namespace MueLu {

std::map<std::string, std::string> myParent_;

template <class TagName>
MutuallyExclusiveTime<TagName>::MutuallyExclusiveTime(const std::string& name, bool startFlag)
  : name_(name)
  , timer_(rcp(new Teuchos::Time(name, false)))
  ,  // second argument is false in any case, because if start==true,
     // timer has to be started by MutuallyExclusiveTime::start() instead of Teuchos::Time::start().
  isPaused_(false) {
  if (startFlag == true) timer_->start();
}

template <class TagName>
MutuallyExclusiveTime<TagName>::~MutuallyExclusiveTime() {
  // This timer can only be destroyed if it is not in the stack
  if (isPaused()) {
    // error message because cannot throw an exception in destructor
    GetOStream(Errors) << "MutuallyExclusiveTime::~MutuallyExclusiveTime(): Error: destructor called on a paused timer." << std::endl;
    // TODO: Even if timing results will be wrong, the timer can be removed from the stack to avoid a segmentation fault.
  }

  stop();  // if isRunning(), remove from the stack, resume previous timer
}

template <class TagName>
void MutuallyExclusiveTime<TagName>::start(bool reset) {
  TEUCHOS_TEST_FOR_EXCEPTION(isPaused(), Exceptions::RuntimeError, "MueLu::MutuallyExclusiveTime::start(): timer is paused. Use resume().");

  if (isRunning()) {
    return;
  }  // If timer is already running, do not pause/push-in-the-stack/start the timer.
     // Otherwise, something bad will happen when this.stop() will be called

  // pause currently running timer
  if (!timerStack_.empty()) {
    GetOStream(Debug) << "pausing parent timer " << timerStack_.top()->name_ << std::endl;
    timerStack_.top()->pause();
    GetOStream(Debug) << "starting child timer " << this->name_ << std::endl;
    myParent_[this->name_] = timerStack_.top()->name_;
  } else {
    GetOStream(Debug) << "starting orphan timer " << this->name_ << std::endl;
    myParent_[this->name_] = "no parent";
  }

  // start this timer
  timer_->start(reset);
  timerStack_.push(this);
}

template <class TagName>
double MutuallyExclusiveTime<TagName>::stop() {
  if (isPaused())
    GetOStream(Errors) << "MueLu::MutuallyExclusiveTime::stop(): timer is paused. Use resume()" << std::endl;

  if (!isRunning()) {
    return timer_->stop();
  }  // stop() can be called on stopped timer

  // Here, timer is running, so it is the head of the stack
  TopOfTheStack();

  timerStack_.pop();
  double r = timer_->stop();

  if (!timerStack_.empty()) {
    GetOStream(Debug) << "resuming timer " << timerStack_.top()->name_ << std::endl;
    timerStack_.top()->resume();
  }

  return r;
}

template <class TagName>
void MutuallyExclusiveTime<TagName>::pause() {
  if (isPaused())  // calling twice pause() is allowed
    return;

  TopOfTheStack();

  timer_->stop();
  isPaused_ = true;
}

template <class TagName>
void MutuallyExclusiveTime<TagName>::resume() {
  TopOfTheStack();

  // no 'shortcut' test necessary:
  // - if timer is stop, it is in pause (cannot be stop and not in pause because this timer is the head of the stack).
  // - if timer is running, nothing is changed by this function.

  timer_->start(false);
  isPaused_ = false;
}

template <class TagName>
bool MutuallyExclusiveTime<TagName>::isRunning() {
  if (timer_->isRunning()) {
    // TEUCHOS_TEST_FOR_EXCEPTION(timerStack_.top() != this, Exceptions::RuntimeError,
    //            "MueLu::MutuallyExclusiveTime::isRunning(): this timer is active so it is supposed to be the head of the stack");
  }
  return timer_->isRunning();
}

template <class TagName>
bool MutuallyExclusiveTime<TagName>::isPaused() {
  TEUCHOS_TEST_FOR_EXCEPTION(isPaused_ && timer_->isRunning(), Exceptions::RuntimeError, "");
  return isPaused_;
}

template <class TagName>
RCP<MutuallyExclusiveTime<TagName> > MutuallyExclusiveTime<TagName>::getNewTimer(const std::string& name) {
  RCP<MutuallyExclusiveTime<TagName> > timer = rcp(new MutuallyExclusiveTime<TagName>(Teuchos::TimeMonitor::getNewTimer(name)));
  timer->name_                               = name;
  return timer;
}

template <class TagName>
void MutuallyExclusiveTime<TagName>::incrementNumCalls() { timer_->incrementNumCalls(); }

template <class TagName>
void MutuallyExclusiveTime<TagName>::PrintParentChildPairs() {
  // key is child, value is parent
  RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  fos->setOutputToRootOnly(0);
  *fos << "Parent Child Map" << std::endl;
  std::map<std::string, std::string>::const_iterator iter;
  for (iter = ::MueLu::myParent_.begin(); iter != ::MueLu::myParent_.end(); ++iter) {
    *fos << "Key: " << iter->first << "  Value: " << iter->second << std::endl;
  }
}

template <class TagName>
MutuallyExclusiveTime<TagName>::MutuallyExclusiveTime(RCP<Teuchos::Time> timer)
  : timer_(timer)
  , isPaused_(false) {}

template <class TagName>
void MutuallyExclusiveTime<TagName>::TopOfTheStack() {
  TEUCHOS_TEST_FOR_EXCEPTION(timerStack_.empty(), Exceptions::RuntimeError, "MueLu::MutuallyExclusiveTime::TopOfTheStack(): timer is not the head of the stack (stack is empty).");
  // TEUCHOS_TEST_FOR_EXCEPTION(timerStack_.top() != this, Exceptions::RuntimeError,  "MueLu::MutuallyExclusiveTime::TopOfTheStack(): timer is not the head of the stack.");
  TEUCHOS_TEST_FOR_EXCEPTION(!(isRunning() || isPaused()), Exceptions::RuntimeError, "MueLu::MutuallyExclusiveTime::TopOfTheStack(): head of the stack timer is neither active nor paused.");
}

template <class TagName>
std::stack<MutuallyExclusiveTime<TagName>*> MutuallyExclusiveTime<TagName>::timerStack_;

// FIXME: move this:
template class MutuallyExclusiveTime<FactoryBase>;
template class MutuallyExclusiveTime<Level>;
template class MutuallyExclusiveTime<BaseClass>;

}  // namespace MueLu
