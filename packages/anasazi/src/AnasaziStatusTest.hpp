// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//

#ifndef ANASAZI_STATUS_TEST_HPP
#define ANASAZI_STATUS_TEST_HPP

/// \file AnasaziStatusTest.hpp
/// \brief Declaration and definition of Anasazi::StatusTest.

#include "AnasaziTypes.hpp"
#include "AnasaziStatusTestDecl.hpp"
#include "AnasaziEigensolverDecl.hpp"

namespace Anasazi {

  //! @name StatusTest Exceptions
  //@{

  /** \brief Exception thrown to signal error in a status test during Anasazi::StatusTest::checkStatus().
   */
  class StatusTestError : public AnasaziError
  {public: StatusTestError(const std::string& what_arg) : AnasaziError(what_arg) {}};

  //@}

/// \example LOBPCGCustomStatusTest.cpp
/// \brief Use LOBPCG with Tpetra, with custom StatusTest.
///
/// This example shows how to define a custom StatusTest so that
/// Anasazi's solver LOBPCG converges correctly with spectrum folding.
/// Without a custom status test, Anasazi would compute the residual
/// as \f$ R = A^2 X - X \Sigma^2\f$.  The custom status test makes
/// Anasazi use the residual \f$ R = A X - X \Sigma\f$ instead.

template <class ScalarType, class MV, class OP>
class StatusTest {

 public:
   //! @name Constructors/destructors
  //@{

  //! Constructor
  StatusTest() {};

  //! Destructor
  virtual ~StatusTest() {};
  //@}

  //! @name Status methods
  //@{
  /*! Check status as defined by test.

    \return TestStatus indicating whether the test passed or failed.
  */
  virtual TestStatus checkStatus( Eigensolver<ScalarType,MV,OP>* solver ) = 0;

  //! Return the result of the most recent checkStatus call, or undefined if it has not been run.
  virtual TestStatus getStatus() const = 0;

  //! Get the indices for the vectors that passed the test.
  virtual std::vector<int> whichVecs() const = 0;

  //! Get the number of vectors that passed the test.
  virtual int howMany() const = 0;

  //@}

  //! @name Reset methods
  //@{
  //! Informs the status test that it should reset its internal configuration to the uninitialized state.
  /*! This is necessary for the case when the status test is being reused by another solver or for another
    eigenvalue problem. The status test may have information that pertains to a particular problem or solver
    state. The internal information will be reset back to the uninitialized state. The user specified information
    that the convergence test uses will remain.
  */
  virtual void reset() = 0;

  //! Clears the results of the last status test.
  /*! This should be distinguished from the reset() method, as it only clears the cached result from the last
   * status test, so that a call to getStatus() will return ::Undefined. This is necessary for the SEQOR and SEQAND
   * tests in the StatusTestCombo class, which may short circuit and not evaluate all of the StatusTests contained
   * in them.
  */
  virtual void clearStatus() = 0;

  //@}

  //! @name Print methods
  //@{

  //! Output formatted description of stopping test to output stream.
  virtual std::ostream& print(std::ostream& os, int indent = 0) const = 0;

  //@}

};

} // end of Anasazi namespace

#endif /* ANASAZI_STATUS_TEST_HPP */
