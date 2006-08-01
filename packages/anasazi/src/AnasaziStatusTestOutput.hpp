
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
    : _printer(printer), _test(test), _state(Undefined), _stateTest(printStates), _modTest(mod), _numCalls(0) {}

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
    _state = _test->checkStatus(solver);

    if (_numCalls % _modTest == 0) {
      if (_state & _stateTest == _state) {
        if ( _printer->isVerbosity(StatusTestDetails) ) {
          _test->print( _printer->stream(StatusTestDetails), 0 );
        }
        else if ( _printer->isVerbosity(Debug) ) {
          _test->print( _printer->stream(Debug), 0 );
        }
      }
    }
    _numCalls++;

    return _state;
  }

  //! Return the result of the most recent checkStatus call, or undefined if it has not been run.
  TestStatus getStatus() const {
    return _state;
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
    _state = Undefined;
    _test->reset();
    _numCalls = 0;
  }

  //! Clears the results of the last status test.
  //! This resets the cached state to an ::Undefined state and calls clearStatus() on the underlying test.
  void clearStatus() {
    _state = Undefined;
    _test->clearStatus();
  }

  //@}

  //! @name Print methods
  //@{ 
  
  //! Output formatted description of stopping test to output stream.
  ostream& print(ostream& os, int indent = 0) const {
    string ind(indent,' ');
    os << ind << "- StatusTestOutput: ";
    switch (_state) {
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
    os << ind << "(Num calls,Mod test,State test): " << "(" << _numCalls << ", " << _modTest << ",";
    if (_stateTest == 0) {
      os << " none )" << endl;
    }
    else {
      if ( _stateTest & Passed == Passed ) os << " Passed";
      if ( _stateTest & Failed == Failed ) os << " Failed";
      if ( _stateTest & Undefined == Undefined ) os << " Undefined";
      os << " )" << endl;
    }
    // print child, with extra indention
    _test->print(os,indent+3);
    return os;
  }
 
  //@}

  private:
    Teuchos::RefCountPtr<OutputManager<ScalarType> > _printer;
    Teuchos::RefCountPtr<StatusTest<ScalarType,MV,OP> > _test;
    TestStatus _state;
    int _stateTest;
    int _modTest;
    int _numCalls;
};

} // end of Anasazi namespace

#endif /* ANASAZI_STATUS_TEST_OUTPUT_HPP */
