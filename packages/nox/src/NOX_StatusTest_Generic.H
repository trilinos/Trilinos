// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_STATUSTEST_GENERIC_H
#define NOX_STATUSTEST_GENERIC_H

#include "NOX_Common.H"        // for std::ostream
#include <ostream>

namespace NOX {

namespace Solver
{
class Generic;            // forward declaration (included below)
}

//! %Status checkers.
namespace StatusTest {

//! Status type.
enum StatusType
  {
    //! Unevaluated
    Unevaluated = -2,
    //! Neither Converged nor Failed
    Unconverged = 0,
    //! Converged
    Converged = 1,
    //! Failed
    Failed = -1
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

/*!

  \brief %Generic status test to check for convergence or failure of
  the nonlinear solver.

*/
class Generic {

public:

  //! Constructor.
  Generic() {};

  //! Destructor.
  virtual ~Generic() {};

  //! %Test the stopping criterion
  /*!
    The test can (and should, if possible) be skipped if
    checkType is NOX::StatusType::None.  If the test is skipped, then
    the status should be set to NOX::StatusTest::Unevaluated.
  */
  virtual NOX::StatusTest::StatusType
  checkStatus(const NOX::Solver::Generic& problem,
          NOX::StatusTest::CheckType checkType) = 0;

  //! Return the result of the most recent checkStatus call
  virtual NOX::StatusTest::StatusType getStatus() const = 0;

  //! Output formatted description of stopping test to output stream.
  virtual std::ostream& print(std::ostream& stream, int indent = 0) const = 0;

};

std::ostream& operator<<(std::ostream& os, NOX::StatusTest::StatusType type);

} // namespace StatusTest
} // namespace NOX


#endif

