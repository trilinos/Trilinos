// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_STATUSTEST_NSTEP_H
#define NOX_STATUSTEST_NSTEP_H

#include "NOX_StatusTest_Generic.H"    // base class
#include "NOX_Abstract_Vector.H"        // for NormType
#include "NOX_Utils.H"                  // for std::cerr

// Forward declaration
namespace NOX {
namespace Abstract {
class Group;
}
}

namespace NOX {

namespace StatusTest {

/** \brief Takes n Iterations before declaring convergence

   This test is used to control an N-step Newton method.  The test
   declares convergence once n nonlinear iterations have been taken.
   Optionally, an initial phase with a different number of steps can
   be taken.

   This is used where the time steps in a transient method are small
   enough that the step is essentially linear.  The typical use case
   is to start with some initial transient run and take enough Newton
   steps to converge to an accurate solutions.  After a few time
   steps, then switch over to a 1-Newton step per time step.
*/
class NStep : public Generic {

public:

  /** \brief Ctor.

    \param[in] numberOfStepsForConvergence Number of nonlinear
    iterations to take until convergence.

    \param[in] numberOfNonlinearSolvesInRampingPhase Number of
    nonlinear solves in an initial ramping phase.  During the ramping
    phase, convergenc will be declared after using a number of
    nonlinear steps equal to rampingPhaseNumberOfStepsForConvergence.

    \param[in] rampingPhaseNumberOfStepsForConvergence Test will
    declaare convergence in the ramping phase after this number of
    time steps.

  */
  NStep(int numberOfStepsForConvergence,
    int numberOfNonlinearSolvesInRampingPhase = 0,
    int rampingPhaseNumberOfStepsForConvergence = 10,
    const NOX::Utils* u = NULL);

  // derived
  virtual NOX::StatusTest::StatusType
  checkStatus(const NOX::Solver::Generic& problem,
          NOX::StatusTest::CheckType checkType);

  // derived
  virtual NOX::StatusTest::StatusType getStatus() const;

  virtual std::ostream& print(std::ostream& stream, int indent = 0) const;

private:

  NOX::StatusTest::StatusType status_;

  //! Ostream used to print errors
  NOX::Utils utils_;

  int numberOfStepsForConvergence_;
  int numberOfNonlinearSolvesInRampingPhase_;
  int rampingPhaseNumberOfStepsForConvergence_;

  int currentNumberOfSteps_;
  int currentNumberOfNonlinearSolves_;
  bool inRampingPhase_;

};

} // namespace Status
} // namespace NOX

#endif
