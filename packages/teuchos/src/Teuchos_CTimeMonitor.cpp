// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#include "Teuchos_CTimeMonitor.h"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_implicit_cast.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_VerboseObject.hpp"


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
      timerArray.push_back(Teuchos::TimeMonitor::getNewTimer(timerName));
      timerArray.back()->start();
      return timerArray.size()-1;
    }
    // Else, the timer already exists so return it
    TEST_FOR_EXCEPTION(
      timerID >=  implicit_cast<int>(timerArray.size()), std::logic_error,
      "Teuchos_startTimer(...): Error, timerID="<<timerID
      <<" is >= timerArray.size()="<<timerArray.size()
      <<" for timerName=\""<<timerName<<"\"!"
      );
    Teuchos::RCP<Teuchos::Time> timer = timerArray[timerID];
    TEST_FOR_EXCEPTION(
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
    TEST_FOR_EXCEPTION(
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
}
