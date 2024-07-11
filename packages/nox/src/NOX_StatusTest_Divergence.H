// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_STATUS_DIVERGENCE_H
#define NOX_STATUS_DIVERGENCE_H

#include "NOX_StatusTest_Generic.H"    // base class

// Forward declarations
namespace NOX {
  namespace Solver {
    class Generic;
  }
}

namespace NOX {

namespace StatusTest {

//! Failure test based on a threshold value of the norm of F.
/*!

  This status test returns NOX::StatusTest::Failed if the norm of F exceeds
  a threshold value for n \e consecutive iterations, i.e.

  \f[ \| F \| > {\rm threshold}\f]

  for n \e consecutive iterations, the status is set to
  NOX::StatusTest::Failed and returned. Otherwise, the status is set
  to NOX::StatusTest::Uncoverged and returned.  Both n and the
  threshold are specified in the constructor, by \c n and \c thresh,
  respectively. While a value for \c thresh must be provided,  a default
  value of \c n = 1 is assumed.

*/
class Divergence : public Generic {

public:

  //! Constructor.
  /*!
    \param thresh - Threshold for divergence test
    \param n - Number of consecutive nonlinear iterations
   */
  Divergence(double thresh, int n = 1);

  //! Destructor.
  virtual ~Divergence();

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
    \brief Returns the user-specified number of steps that can
    consecutively fail the threshold test before the test returns a
    failed status.
  */
  virtual int getMaxNumSteps() const;

  //! Returns the current number of steps that have consecutively failed the threshold test.
  virtual int getCurrentNumSteps() const;

  //! Returns the user-specified threshold.
  virtual double getThreshold() const;

private:

  //! User supplied value of n.
  int maxSteps;

  //! Current number of consecutive nonlinear iterations that have failed the specified threshold.
  int numSteps;

  //! The last nonlinear iteration index.
  /*! This is used to prevent counting a step multiple times if by chance
    the status test is called multiple times between iterations.
  */
  int lastIteration;

  //! User specified threshold.
  double threshold;

  //! %Status
  NOX::StatusTest::StatusType status;


};

} // namespace Status
} // namespace NOX

#endif
