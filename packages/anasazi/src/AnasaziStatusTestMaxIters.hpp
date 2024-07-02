// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
    <pre>
                     { Passed,  if solver->getNumIters() >= maxIter
    status(solver) = {
                     { Failed,  if solver->getNumIters()  < maxIter
    </pre>
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
  StatusTestMaxIters(int maxIter, bool negate = false) : state_(Undefined), negate_(negate) {
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
    state_ = (solver->getNumIters() >= maxIters_) ? Passed : Failed;
    if (negate_) {
      if (state_ == Passed) state_ = Failed;
      else state_ = Passed;
    }
    return state_;
  }

  //! \brief Return the result of the most recent checkStatus call.
  TestStatus getStatus() const {
    return state_;
  }

  //! Get the indices for the vectors that passed the test.
  std::vector<int> whichVecs() const {
    return std::vector<int>(0);
  }

  //! Get the number of vectors that passed the test.
  int howMany() const {
    return 0;
  }

  //@}

  //! @name Accessor methods
  //@{ 

  /*! \brief Set the maximum number of iterations.
   *  \note This also resets the test status to ::Undefined.
   */
  void setMaxIters(int maxIters) {
    state_ = Undefined;
    maxIters_ = maxIters;
  }

  //! \brief Get the maximum number of iterations.
  int getMaxIters() {return maxIters_;}

  /*! \brief Set the negation policy for the status test.
   *  \note This also reset the test status to ::Undefined.
   */
  void setNegate(bool negate) {
    state_ = Undefined;
    negate_ = negate;
  }

  //! \brief Get the negation policy for the status test.
  bool getNegate() const {
    return negate_;
  }

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
    state_ = Undefined;
  }

  //! \brief Clears the results of the last status test.
  /*! This should be distinguished from the reset() method, as it only clears the cached result from the last 
   * status test, so that a call to getStatus() will return ::Undefined. This is necessary for the SEQOR and SEQAND
   * tests in the StatusTestCombo class, which may short circuit and not evaluate all of the StatusTests contained
   * in them.
  */
  void clearStatus() {
    state_ = Undefined;
  }

  //@}

  //! @name Print methods
  //@{ 
  
  //! Output formatted description of stopping test to output stream.
  std::ostream& print(std::ostream& os, int indent = 0) const {
    std::string ind(indent,' ');
    os << ind << "- StatusTestMaxIters: ";
    switch (state_) {
    case Passed:
      os << "Passed" << std::endl;
      break;
    case Failed:
      os << "Failed" << std::endl;
      break;
    case Undefined:
      os << "Undefined" << std::endl;
      break;
    }
    os << ind << "  MaxIters: " << maxIters_ << std::endl;
    return os;
  }
 
  //@}
  private:
  int maxIters_;
  TestStatus state_;
  bool negate_;

};

} // end of Anasazi namespace

#endif /* ANASAZI_STATUS_TEST_MAXITER_HPP */
