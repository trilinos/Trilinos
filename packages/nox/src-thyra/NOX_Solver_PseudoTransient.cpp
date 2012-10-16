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

#include "NOX_Solver_PseudoTransient.hpp"	// class definition
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Common.H"
#include "Teuchos_ParameterList.hpp"
#include "NOX_Utils.H"
#include "NOX_GlobalData.H"
#include "NOX_Solver_SolverUtils.H"
#include "NOX_LineSearch_Generic.H"
#include "NOX_LineSearch_Factory.H"
#include "NOX_Direction_Generic.H"
#include "NOX_Direction_Factory.H"

// This relies explicitly on thyra
#include "NOX_Thyra_Group.H"
#include "Thyra_VectorBase.hpp"
#include "Thyra_ModelEvaluatorBase.hpp"

NOX::Solver::PseudoTransient::
PseudoTransient(const Teuchos::RCP<NOX::Abstract::Group>& xGrp, 
		const Teuchos::RCP<NOX::StatusTest::Generic>& t, 
		const Teuchos::RCP<Teuchos::ParameterList>& p) :
  globalDataPtr(Teuchos::rcp(new NOX::GlobalData(p))),
  utilsPtr(globalDataPtr->getUtils()), 
  solnPtr(xGrp),		                       // pointer to xGrp
  oldSolnPtr(xGrp->clone(DeepCopy)),     // create via clone
  dirPtr(xGrp->getX().clone(ShapeCopy)), // create via clone 
  testPtr(t),	
  paramsPtr(p),   
  prePostOperator(utilsPtr, paramsPtr->sublist("Solver Options"))
{
  thyraSolnGroup = Teuchos::rcp_dynamic_cast<NOX::Thyra::Group>(solnPtr,true);
  thyraOldSolnGroup = Teuchos::rcp_dynamic_cast<NOX::Thyra::Group>(oldSolnPtr,true);

  init();
}

// Protected
void NOX::Solver::PseudoTransient::init()
{
  // Initialize 
  stepSize = 0.0;
  nIter = 0;
  status = NOX::StatusTest::Unconverged;
  checkType = parseStatusTestCheckType(paramsPtr->sublist("Solver Options"));

  lineSearchPtr = NOX::LineSearch::
    buildLineSearch(globalDataPtr, paramsPtr->sublist("Line Search"));

  directionPtr = NOX::Direction::
    buildDirection(globalDataPtr, paramsPtr->sublist("Direction"));

  deltaInit = paramsPtr->sublist("Pseudo-Transient").get<double>("deltaInit",1.0e-4);
  delta = deltaInit;
  inv_delta = 1.0 / delta;
  deltaMax = paramsPtr->sublist("Pseudo-Transient").get<double>("deltaMax",1.0e+0);
  deltaMin = paramsPtr->sublist("Pseudo-Transient").get<double>("deltaMin",1.0e-5);
  time = 0.0;

  // Print out parameters
  if (utilsPtr->isPrintType(NOX::Utils::Parameters)) 
  {
    utilsPtr->out() << "\n" << NOX::Utils::fill(72) << "\n";
    utilsPtr->out() << "\n-- Parameters Passed to Nonlinear Solver --\n\n";
    paramsPtr->print(utilsPtr->out(),5);
  }

}

void NOX::Solver::PseudoTransient::
reset(const NOX::Abstract::Vector& initialGuess, 
      const Teuchos::RCP<NOX::StatusTest::Generic>& t)
{
  solnPtr->setX(initialGuess);
  testPtr = t;
  init();
}

void NOX::Solver::PseudoTransient::
reset(const NOX::Abstract::Vector& initialGuess)
{
  solnPtr->setX(initialGuess);
  init();
}

NOX::StatusTest::StatusType NOX::Solver::PseudoTransient::getStatus()
{
  return status;
}

NOX::StatusTest::StatusType NOX::Solver::PseudoTransient::step()
{
  prePostOperator.runPreIterate(*this);

  // On the first step, do some initializations
  if (nIter == 0) {
    // Compute F of initital guess
    NOX::Abstract::Group::ReturnType rtype = solnPtr->computeF();
    if (rtype != NOX::Abstract::Group::Ok) {
      utilsPtr->out() << "NOX::Solver::PseudoTransient::init - "
		      << "Unable to compute F" << endl;
      throw "NOX Error";
    }

    // Test the initial guess
    status = testPtr->checkStatus(*this, checkType);
    if ((status == NOX::StatusTest::Converged) &&
	(utilsPtr->isPrintType(NOX::Utils::Warning))) {
      utilsPtr->out() << "Warning: NOX::Solver::PseudoTransient::init() - "
		      << "The solution passed into the solver (either "
		      << "through constructor or reset method) "
		      << "is already converged!  The solver wil not "
		      << "attempt to solve this system since status is "
		      << "flagged as converged." << endl;
    }

    printUpdate();
  }

  // First check status
  if (status != NOX::StatusTest::Unconverged) {
    prePostOperator.runPostIterate(*this);
    printUpdate();
    return status;
  }

  // Copy pointers into temporary references
  NOX::Abstract::Group& soln = *solnPtr;
  NOX::StatusTest::Generic& test = *testPtr;

  // Pseudo-transient: change the Jacobian evaluation to evaluate a transient version
  {
    deltaOld = delta;
    
    // Update step size
    if (nIter == 0)
      delta = deltaInit;
    else
      delta = deltaOld * oldSolnPtr->getNormF() / solnPtr->getNormF();
    
    inv_delta = 1.0 / delta;
    if (delta > deltaMax)
      inv_delta = 0.0;

    if (delta < deltaMin) {
      delta = deltaMin;
      inv_delta = 1.0 / delta;
    }

    time += delta;
    
    Teuchos::RCP<const ::Thyra::VectorBase<double> > x = 
      thyraSolnGroup->get_current_x();
    
    Teuchos::RCP<const ::Thyra::VectorBase<double> > x_old = 
      thyraOldSolnGroup->get_current_x();

    // Compute x_dot using forward finite difference
    if (is_null(x_dot))
      x_dot = ::Thyra::createMember(x->space());

    if (nIter == 0)
      ::Thyra::put_scalar(0.0,x_dot.ptr());
    else {
      ::Thyra::V_StVpStV(x_dot.ptr(),inv_delta,*x,-inv_delta,*x_old);
    }

    ::Thyra::ModelEvaluatorBase::InArgs<double>& inArgs = thyraSolnGroup->getNonconstInArgs();

    inArgs.set_x_dot(x_dot);
    inArgs.set_alpha(inv_delta);
    inArgs.set_beta(1.0);
    inArgs.set_t(time);
  }

  // Compute the direction for the update vector at the current solution.
  // F is already computed so the only thing to compute is J using our augmented inargs
  bool ok = true;
  ok = directionPtr->compute(*dirPtr, soln, *this);
  if (!ok) 
  {
    utilsPtr->out() << "NOX::Solver::PseudoTransient::iterate - unable to calculate direction" << endl;
    status = NOX::StatusTest::Failed;
    prePostOperator.runPostIterate(*this);
    printUpdate();
    return status;
  }

  // reset the inargs to the correct value for steady state residual evaluations
  {
    ::Thyra::ModelEvaluatorBase::InArgs<double>& inArgs = thyraSolnGroup->getNonconstInArgs();
    inArgs.set_x_dot(Teuchos::null);
    inArgs.set_alpha(0.0);
    inArgs.set_beta(1.0);
    inArgs.set_t(0.0);
  }


  // Update iteration count.
  nIter ++;

  // Copy current soln to the old soln.
  *oldSolnPtr = *solnPtr;

  // Do line search and compute new soln.
  ok = lineSearchPtr->compute(soln, stepSize, *dirPtr, *this);
  if (!ok) 
  {
    if (stepSize == 0.0) 
    {
      utilsPtr->out() << "NOX::Solver::PseudoTransient::iterate - line search failed" << endl;
      status = NOX::StatusTest::Failed;
      prePostOperator.runPostIterate(*this);
      printUpdate();
      return status;
    }
    else if (utilsPtr->isPrintType(NOX::Utils::Warning))
      utilsPtr->out() << "NOX::Solver::PseudoTransient::iterate - using recovery step for line search" << endl;
  }

  // Compute F for new current solution.
  NOX::Abstract::Group::ReturnType rtype = soln.computeF();
  if (rtype != NOX::Abstract::Group::Ok) 
  {
    utilsPtr->out() << "NOX::Solver::PseudoTransient::iterate - unable to compute F" << endl;
    status = NOX::StatusTest::Failed;
    prePostOperator.runPostIterate(*this);
    printUpdate();
    return status;
  }

  // Evaluate the current status.
  status = test.checkStatus(*this, checkType);
 
  prePostOperator.runPostIterate(*this);

  printUpdate();

  return status;
}

NOX::StatusTest::StatusType NOX::Solver::PseudoTransient::solve()
{
  prePostOperator.runPreSolve(*this);

  // Iterate until converged or failed
  while (status == NOX::StatusTest::Unconverged)
    step();

  Teuchos::ParameterList& outputParams = paramsPtr->sublist("Output");
  outputParams.set("Nonlinear Iterations", nIter);
  outputParams.set("2-Norm of Residual", solnPtr->getNormF());

  prePostOperator.runPostSolve(*this);

  return status;
}

const NOX::Abstract::Group& 
NOX::Solver::PseudoTransient::getSolutionGroup() const
{
  return *solnPtr;
}

const NOX::Abstract::Group& 
NOX::Solver::PseudoTransient::getPreviousSolutionGroup() const
{
  return *oldSolnPtr;
}

int NOX::Solver::PseudoTransient::getNumIterations() const
{
  return nIter;
}

const Teuchos::ParameterList& 
NOX::Solver::PseudoTransient::getList() const
{
  return *paramsPtr;
}

// protected
void NOX::Solver::PseudoTransient::printUpdate() 
{
  double normSoln = 0;
  double normStep = 0;

  // Print the status test parameters at each iteration if requested  
  if ((status == NOX::StatusTest::Unconverged) && 
      (utilsPtr->isPrintType(NOX::Utils::OuterIterationStatusTest))) 
  {
    utilsPtr->out() << NOX::Utils::fill(72) << "\n";
    utilsPtr->out() << "-- Status Test Results --\n";    
    testPtr->print(utilsPtr->out());
    utilsPtr->out() << NOX::Utils::fill(72) << "\n";
  }

  // All processes participate in the computation of these norms...
  if (utilsPtr->isPrintType(NOX::Utils::OuterIteration)) 
  {
    normSoln = solnPtr->getNormF();
    normStep = (nIter > 0) ? dirPtr->norm() : 0;
  }

  // ...But only the print process actually prints the result.
  if (utilsPtr->isPrintType(NOX::Utils::OuterIteration)) 
  {
    utilsPtr->out() << "\n" << NOX::Utils::fill(72) << "\n";
    utilsPtr->out() << "-- Nonlinear Solver Step " << nIter << " -- \n";
    utilsPtr->out() << "||F|| = " << utilsPtr->sciformat(normSoln);
    utilsPtr->out() << "  step = " << utilsPtr->sciformat(stepSize);
    utilsPtr->out() << "  dx = " << utilsPtr->sciformat(normStep);
    utilsPtr->out() << "  delta = " << utilsPtr->sciformat(delta);
    utilsPtr->out() << "  inv_delta = " << utilsPtr->sciformat(inv_delta);
    if (status == NOX::StatusTest::Converged)
      utilsPtr->out() << " (Converged!)";
    if (status == NOX::StatusTest::Failed)
      utilsPtr->out() << " (Failed!)";
    utilsPtr->out() << "\n" << NOX::Utils::fill(72) << "\n" << endl;
  }

  // Print the final parameter values of the status test
  if ((status != NOX::StatusTest::Unconverged) && 
      (utilsPtr->isPrintType(NOX::Utils::OuterIteration))) 
  {
    utilsPtr->out() << NOX::Utils::fill(72) << "\n";
    utilsPtr->out() << "-- Final Status Test Results --\n";    
    testPtr->print(utilsPtr->out());
    utilsPtr->out() << NOX::Utils::fill(72) << "\n";
  }
}

double NOX::Solver::PseudoTransient::getStepSize() const
{
  return stepSize;
}
