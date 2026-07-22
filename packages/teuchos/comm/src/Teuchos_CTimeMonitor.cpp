// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_CTimeMonitor.h"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_implicit_cast.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_CompilerCodeTweakMacros.hpp"


namespace {


typedef Teuchos::Array< Teuchos::RCP<Teuchos::Time> >  TimerArray_t;
TimerArray_t timerArray;


} // namespace


int Teuchos_startTimer( char timerName[], int timerID )
{
  using Teuchos::implicit_cast;
  bool success = true;
  try {
    if( timerID < 0 ) {
      // The timer does not exist so create it!
      timerArray.push_back(Teuchos::TimeMonitor::getNewCounter(timerName));
      timerArray.back()->start();
      return timerArray.size()-1;
    }
    // Else, the timer already exists so return it
    TEUCHOS_TEST_FOR_EXCEPTION(
      timerID >=  implicit_cast<int>(timerArray.size()), std::logic_error,
      "Teuchos_startTimer(...): Error, timerID="<<timerID
      <<" is >= timerArray.size()="<<timerArray.size()
      <<" for timerName=\""<<timerName<<"\"!"
      );
    Teuchos::RCP<Teuchos::Time> timer = timerArray[timerID];
    TEUCHOS_TEST_FOR_EXCEPTION(
      timer->isRunning(), std::logic_error,
      "Teuchos_startTimer(...): Error, timerID="<<timerID
      <<", timerName=\""<<timerName<<"\" is already running!"
      );
    timer->start();
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,
    *Teuchos::VerboseObjectBase::getDefaultOStream(), success);
  if (!success) {
    return -1;
  }
  return timerID;
}


void Teuchos_stopTimer( int timerID )
{
  using Teuchos::implicit_cast;
  bool success = true;
  try {
    TEUCHOS_TEST_FOR_EXCEPTION(
      timerID < 0 || timerID >= implicit_cast<int>(timerArray.size()),
      std::logic_error,
      "Teuchos_stopTimer(...): Error, timerID="<<timerID<<" is invalid!"
      );
    Teuchos::RCP<Teuchos::Time> timer = timerArray[timerID];
    timer->stop();
    // Increment the number of times the timer has been used (start to stop).
    timer->incrementNumCalls();
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,
    *Teuchos::VerboseObjectBase::getDefaultOStream(), success);
  if (!success) {} // Avoid warnings
}
