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
  (void)success; // To avoid wrong compiler warning on GCC 4.6.1 
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
  (void)success; // GCC 4.6.1 says this variable is unused?
}
