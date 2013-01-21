// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#ifndef NOX_STATUSTEST_NSTEP_H
#define NOX_STATUSTEST_NSTEP_H

#include "NOX_StatusTest_Generic.H"	// base class
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
