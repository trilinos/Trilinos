
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

#ifndef ANASAZI_STATUS_TEST_MAXITER_HPP
#define ANASAZI_STATUS_TEST_MAXITER_HPP

/*!
  \file AnasaziStatusTestMaxIters.hpp
  \brief Status test for testing the number of iterations.
*/


#include "AnasaziStatusTest.hpp"


  /*! 
    \class Anasazi::StatusTestMaxIters
    \brief A status test for testing the number of iterations.
    
    Anasazi::StatusTestMaxIters will test true when an eigensolver has reached some number 
    of iterations. Specifically, 
    <div>
                     { true,   if solver->getNumIters() >= maxIter
    status(solver) = {
                     { false,  if solver->getNumIters()  < maxIter
    </div>
    where maxIter is the parameter given to the status tester.
    
    This status test also supports negation, so that it negates the need for a 
    StatusTestMinIters status tester. In this way, all tests on the range of iterations
    can be constructed through the appropriate use of StatusTestMaxIters and StatusTestCombo.
  */

namespace Anasazi {


template <class ScalarType, class MV, class OP>
class StatusTestMaxIters : public StatusTest<ScalarType,MV,OP> {

 public:
   //! @name Constructors/destructors
  //@{ 

  //! Constructor
  StatusTestMaxIters(int maxIter, bool negate = false) : _state(Undefined), _negate(negate) {
    setMaxIters(maxIter);
  };

  //! Destructor
  virtual ~StatusTestMaxIters() {};
  //@}

  //! @name Status methods
  //@{ 

  /*! \brief Check status as defined by test.
    \return TestStatus indicating whether the test passed or failed.
  */
  TestStatus checkStatus( Eigensolver<ScalarType,MV,OP>* solver ) {
    _state = (solver->getNumIters() >= _maxIters) ? Passed : Failed;
    if (_negate) {
      if (_state == Passed) _state = Failed;
      else _state = Passed;
    }
    return _state;
  }

  //! \brief Return the result of the most recent checkStatus call.
  TestStatus getStatus() const {
    return _state;
  }
  //@}

  //! @name Accessor methods
  //@{ 

  /*! \brief Set the maximum number of iterations.
   *  This also resets the test status to Undefined.
   */
  void setMaxIters(int maxIters) {
    _state = Undefined;
    _maxIters = maxIters;
  }

  //! \brief Get the maximum number of iterations.
  int getMaxIters() {return _maxIters;}

  //@}

  //! @name Reset methods
  //@{ 
  //! Informs the status test that it should reset its internal configuration to the uninitialized state.
  /*! The StatusTestMaxIters class has no internal state, so this call is equivalent to calling clearStatus().
    eigenvalue problem. The status test may have information that pertains to a particular problem or solver 
    state. The internal information will be reset back to the uninitialized state. The user specified information 
    that the convergence test uses will remain.
  */
  void reset() {
    _state = Undefined;
  }

  //! \brief Clears the results of the last status test.
  /*! This should be distinguished from the reset() method, as it only clears the cached result from the last 
   * status test, so that a call to getStatus() will return Undefined. This is necessary for the SEQOR and SEQAND
   * tests in the StatusTestCombo class, which may short circuit and not evaluate all of the StatusTests contained
   * in them.
  */
  void clearStatus() {
    _state = Undefined;
  }

  //@}

  //! @name Print methods
  //@{ 
  
  //! Output formatted description of stopping test to output stream.
  ostream& print(ostream& os, int indent = 0) const {
    string ind(indent,' ');
    os << ind << "solver->getNumIters() >= " << _maxIters 
       << (_negate ? " tests " : " negates ");
    switch(_state){
      case Passed:
        os << "Passed" << endl;
        break;
      case Undefined:
        os << "Undefined" << endl;
        break;
      case Failed:
        os << "Failed" << endl;
        break;
    }
    return os;
  }
 
  //@}
  private:
  int _maxIters;
  TestStatus _state;
  bool _negate;

};

} // end of Anasazi namespace

#endif /* ANASAZI_STATUS_TEST_MAXITER_HPP */
