// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_STATUS_STAGNATION_H
#define NOX_STATUS_STAGNATION_H

#include "NOX_StatusTest_Generic.H"    // base class

// Forward declarations
namespace NOX {
  namespace Solver {
    class Generic;
  }
}

namespace NOX {

namespace StatusTest {

//! Failure test based on the convergence rate between nonlinear iterations.
/*!

  This status test returns NOX::StatusTest::Failed if we fail to
  reduce the norm of \f$ F\f$ by a specified tolerance for n
  \e consecutive iterations.
  In other words, if

  \f[ \frac{\| F_k \|}{\| F_{k-1} \|} \geq {\rm tolerance}\f]

  for n \e consecutive iterations, the status is set to
  NOX::StatusTest::Failed and returned. Otherwise, the status is set
  to NOX::StatusTest::Uncoverged and returned.  Both n and the
  tolerance are specified in the constructor, by \c n and \c tol,
  respectively.

  Based on experience the following values are recommended:

  <ul>
  <li> For Newton solves:  n = 50, tolerance = 1.0<br>
  <li> For Newton solves with a line search: n = 15, tolerance = 0.99
  </ul>
*/
class Stagnation : public Generic {

public:

  //! Constructor.
  /*!
    \param n - Number of consecutive nonlinear iterations
    \param tol - Tolerance for stagnation test
   */
  Stagnation(int n = 50, double tol = 1.0);

  //! Destructor.
  virtual ~Stagnation();

  /*!
    \brief Tests the stopping criterion.

    \note
    The nature of this test is such that it \e must be executed at
    every nonlinear iteration, so we don't use the checkType argument.
   */
  virtual NOX::StatusTest::StatusType
  checkStatus(const NOX::Solver::Generic& problem,
          NOX::StatusTest::CheckType checkType);

  virtual NOX::StatusTest::StatusType getStatus() const;

  virtual std::ostream& print(std::ostream& stream, int indent = 0) const;

  /*!
    \brief Returns the used specified number of steps that can
    consecutively fail the tolerance test before the test returns a
    failed status.
  */
  virtual int getMaxNumSteps() const;

  //! Returns the current number of steps that have consecutively failed the tolerance test.
  virtual int getCurrentNumSteps() const;

  //! Returns the user specified tolerance.
  virtual double getTolerance() const;

  //! Returns the current convergence rate.
  virtual double getConvRate() const;

private:

  //! User supplied value of n.
  int maxSteps;

  //! Current number of consecutive nonlinear iterations that have failed the specified tolerance.
  int numSteps;

  //! The last nonlinear iteration index.
  /*! This is used to prevent counting a step multiple times if by chance
    the status test is called multiple times between iterations.
  */
  int lastIteration;

  //! User specified tolerance.
  double tolerance;

  //! Currently computed convergence rate.
  double convRate;

  //! %Status
  NOX::StatusTest::StatusType status;


};

} // namespace Status
} // namespace NOX

#endif
