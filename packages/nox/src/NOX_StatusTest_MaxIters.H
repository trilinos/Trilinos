// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_STATUS_MAXITERS_H
#define NOX_STATUS_MAXITERS_H

#include "NOX_StatusTest_Generic.H"    // base class
#include "NOX_Utils.H"                  // class data element

namespace NOX {

namespace StatusTest {

//! Failure test based on the maximum number of nonlinear solver iterations.
/*!
  Let \f$ k\f$ denote the current number of iterations (accessed via
  NOX::Solver::getNumIterations) and \f$ k_{\max}\f$ denote the
  tolerance set in the constructor of this status test.
  This test returns NOX::StatusTest::Failed if
  \f$
  k \geq k_{\rm max}.
  \f$
  Otherwise, it returns NOX::StatusTest::Unconverged.

  \note If checkStatus is called with the type set to
  NOX::StatusTest::None, it then the status is set to to
  NOX::Status::Unevaluated and returned. (Also #niters is set to -1.)
*/
class MaxIters : public Generic {

public:

  //! Constructor. Specify the maximum number of nonlinear solver iterations, \f$ k_{\max}\f$ ands optinally an error stream for printing errors.
  MaxIters(int maxIterations, const NOX::Utils* u = NULL);

  //! Destructor.
  virtual ~MaxIters();

  virtual NOX::StatusTest::StatusType
  checkStatus(const NOX::Solver::Generic& problem,
          NOX::StatusTest::CheckType checkType);

  virtual NOX::StatusTest::StatusType getStatus() const;

  virtual std::ostream& print(std::ostream& stream, int indent = 0) const;

  //! Returns the Maximum number of iterations set in the constructor.
  virtual int getMaxIters() const;

  /*!
    \brief Returns the current number of iterations taken by the solver.

    Returns -1 if the status of this test is NOX::StatusTest::Unevaluated.
   */
  virtual int getNumIters() const;

private:

  //! Maximum number of iterations
  int maxiters;

  //! Current number of iterations (if known)
  int niters;

  //! %Status
  NOX::StatusTest::StatusType status;

  //! Ostream used to print errors
  NOX::Utils utils;

};

} // namespace Status
} // namespace NOX

#endif
