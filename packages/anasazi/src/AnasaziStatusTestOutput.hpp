
// @HEADER
// ***********************************************************************
//
//                 Belos: Block Linear Solvers Package
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
//

#ifndef ANASAZI_STATUS_TEST_OUTPUT_HPP
#define ANASAZI_STATUS_TEST_OUTPUT_HPP

/*!
  \file AnasaziStatusTestOutput.hpp
  \brief Special StatusTest for printing status tests.
*/


#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"
#include "AnasaziEigensolver.hpp"

#include "AnasaziStatusTest.hpp"


  /*! 
    \class Anasazi::StatusTestOutput
    \brief A special Anasazi::StatusTest for printing other status tests. 
    
    Anasazi::StatusTestOutput is a wrapper around another Anasazi::StatusTest that calls 
    Anasazi::StatusTest::print() on the underlying object on calls to Anasazi::StatusTestOutput::checkStatus().
    The frequency and occasion of the printing can be dictated according to some parameters passed to the 
    Anasazi::StatusTestOutput class.
  */

namespace Anasazi {

template <class ScalarType, class MV, class OP>
class StatusTestOutput : public StatusTest<ScalarType,MV,OP> {

 public:
  //! @name Constructors/destructors
  //@{ 

  /*! \brief Constructor
   *
   * The StatusTestOutput requires an OutputManager for printing the underlying StatusTest on
   * calls to checkStatus(), as well as an underlying StatusTest.
   *
   * The other parameters, in addition to the verbosity level of the OutputManager, control when printing is 
   * called. 
   *
   * @param[in] mod A positive number describes how often the output should be printed. Default: 1 (attempt to print on every call to checkStatus())
   * @param[in] printStates A sum of TestStatus returns for which the output should be printed. Default: ::Passed (attempt to print whenever checkStatus() will return ::Passed)
   *
   */
  StatusTestOutput(const Teuchos::RefCountPtr<OutputManager<ScalarType> > &printer, 
                   Teuchos::RefCountPtr<StatusTest<ScalarType,MV,OP> > test,
                   int mod = 1,
                   int printStates = Passed)
    : printer_(printer), test_(test), state_(Undefined), stateTest_(printStates), modTest_(mod), numCalls_(0) {}

  //! Destructor
  virtual ~StatusTestOutput() {};
  //@}

  //! @name Status methods
  //@{ 
  /*! Check and return status of underlying StatusTest.

    This method calls checkStatus() on the StatusTest object passed in the constructor. If appropriate, the
    method will follow this call with a call to print() on the underlying object, using the OutputManager passed via the constructor
    with verbosity level ::StatusTestDetails.

    The internal counter will be incremented during this call, but only after
    performing the tests to decide whether or not to print the underlying
    StatusTest. This way, the very first call to checkStatus() following
    initialization or reset() will enable the underlying StatusTest to be
    printed, regardless of the mod parameter, as the current number of calls
    will be zero.
    
    \return TestStatus indicating whether the underlying test passed or failed.
  */
  TestStatus checkStatus( Eigensolver<ScalarType,MV,OP>* solver ) {
    state_ = test_->checkStatus(solver);

    if (numCalls_ % modTest_ == 0) {
      if ( (state_ & stateTest_) == state_) {
        if ( printer_->isVerbosity(StatusTestDetails) ) {
          test_->print( printer_->stream(StatusTestDetails), 0 );
        }
        else if ( printer_->isVerbosity(Debug) ) {
          test_->print( printer_->stream(Debug), 0 );
        }
      }
    }
    numCalls_++;

    return state_;
  }

  //! Return the result of the most recent checkStatus call, or undefined if it has not been run.
  TestStatus getStatus() const {
    return state_;
  }
  //@}

  //! @name Reset methods
  //@{ 
  /*! \brief Informs the status test that it should reset its internal configuration to the uninitialized state.
   *
   *  This resets the cached state to an ::Undefined state and calls reset() on the underlying test. It also 
   *  resets the counter for the number of calls to checkStatus().
   */
  void reset() { 
    state_ = Undefined;
    test_->reset();
    numCalls_ = 0;
  }

  //! Clears the results of the last status test.
  //! This resets the cached state to an ::Undefined state and calls clearStatus() on the underlying test.
  void clearStatus() {
    state_ = Undefined;
    test_->clearStatus();
  }

  //@}

  //! @name Print methods
  //@{ 
  
  //! Output formatted description of stopping test to output stream.
  ostream& print(ostream& os, int indent = 0) const {
    string ind(indent,' ');
    os << ind << "- StatusTestOutput: ";
    switch (state_) {
    case Passed:
      os << "Passed" << endl;
      break;
    case Failed:
      os << "Failed" << endl;
      break;
    case Undefined:
      os << "Undefined" << endl;
      break;
    }
    os << ind << "(Num calls,Mod test,State test): " << "(" << numCalls_ << ", " << modTest_ << ",";
    if (stateTest_ == 0) {
      os << " none )" << endl;
    }
    else {
      if ( stateTest_ & Passed == Passed ) os << " Passed";
      if ( stateTest_ & Failed == Failed ) os << " Failed";
      if ( stateTest_ & Undefined == Undefined ) os << " Undefined";
      os << " )" << endl;
    }
    // print child, with extra indention
    test_->print(os,indent+3);
    return os;
  }
 
  //@}

  private:
    Teuchos::RefCountPtr<OutputManager<ScalarType> > printer_;
    Teuchos::RefCountPtr<StatusTest<ScalarType,MV,OP> > test_;
    TestStatus state_;
    int stateTest_;
    int modTest_;
    int numCalls_;
};

} // end of Anasazi namespace

#endif /* ANASAZI_STATUS_TEST_OUTPUT_HPP */
