// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_STATUSTEST_ABSTRACT_H
#define LOCA_STATUSTEST_ABSTRACT_H

// FIXME Replace by LOCA equivalent
#include "NOX_Common.H"         // for std::ostream
#include "LOCA_Abstract_Iterator.H"
#include <ostream>

namespace LOCA {


  // forward declaration of LOCA::Stepper
  // Modify interfaces using LOCA::Stepper in favor of its LOCA::Abstract::Iterator parent (GAH)
//  class Stepper;

//! %Status checkers.
namespace StatusTest {

/*!

  \brief %Generic status test to check for convergence or failure of
  the LOCA continuation.

*/

  //! Enumerated type for status of the iterator
  enum StatusType {
    Finished = 1,       //!< The iterator is finished
    Failed = 0,         //!< The iterator has failed
    NotFinished = -1,   //!< The iterator is not finished
    Unevaluated = -2    //!< The iterator status has not been evaluated
  };

  //! Type of check that should be done by checkStatus
  enum CheckType
    {
      //! Evaluate every test and subtest
      Complete,
      //! Evaluate all tests and subtests that are necessary to determine the status
      Minimal,
      //! Evaluation is entirely optional
      None
    };

class Abstract {

public:

  //! Constructor.
  Abstract() {};

  //! Destructor.
  virtual ~Abstract() {};

  //! %Test the stopping criterion
  /*!
    The test can (and should, if possible) be skipped if
    checkType is LOCA::StatusType::None.  If the test is skipped, then
    the status should be set to LOCA::StatusTest::Unevaluated.
  */
  // Modify interfaces using LOCA::Stepper in favor of its LOCA::Abstract::Iterator parent (GAH)
  virtual LOCA::StatusTest::StatusType
//  checkStatus(const LOCA::Stepper& stepper,
  checkStatus(const LOCA::Abstract::Iterator& stepper,
            LOCA::StatusTest::CheckType checkType) = 0;

  //! Return the result of the most recent checkStatus call
  virtual LOCA::StatusTest::StatusType
  getStatus() const = 0;

  //! Output formatted description of stopping test to output stream.
  virtual std::ostream&
  print(std::ostream& stream, int indent = 0) const = 0;

};

std::ostream&
operator<<(std::ostream& os, LOCA::StatusTest::StatusType status);

} // namespace StatusTest
} // namespace LOCA


#endif
