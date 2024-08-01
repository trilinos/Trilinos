// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "NOX_Direction_NonlinearCG.H"    // class definition
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Common.H"
#include "Teuchos_ParameterList.hpp"
#include "NOX_Solver_Generic.H"
#include "NOX_Utils.H"
#include "NOX_GlobalData.H"

using namespace NOX;
using namespace NOX::Direction;

NonlinearCG::NonlinearCG(const Teuchos::RCP<NOX::GlobalData>& gd,
             Teuchos::ParameterList& params) :
  oldSolnPtr(NULL),
  paramsPtr(0),
  beta(0.0),
  niter(0)
{
  reset(gd, params);
}


bool NonlinearCG::
reset(const Teuchos::RCP<NOX::GlobalData>& gd,
      Teuchos::ParameterList& params)
{
  globalDataPtr = gd;
  utils = gd->getUtils();
  paramsPtr = &params;
  Teuchos::ParameterList& nlcgParams = paramsPtr->sublist("Nonlinear CG");
  restartFrequency = nlcgParams.get("Restart Frequency", 10);
  doPrecondition = false;
  if(  nlcgParams.get("Precondition", "Off") == "On" )
    doPrecondition = true;
  usePRbeta = false;
  if( nlcgParams.get("Orthogonalize", "Fletcher-Reeves") ==  "Polak-Ribiere" )
    usePRbeta = true;

  return true;
}

NonlinearCG::~NonlinearCG()
{
}


bool NonlinearCG::compute(Abstract::Vector& dir, Abstract::Group& soln,
                          const Solver::Generic& solver)
{
  Abstract::Group::ReturnType ok;

  // Initialize vector memory if haven't already
  if(Teuchos::is_null(oldDirPtr))
    oldDirPtr = soln.getX().clone(NOX::ShapeCopy);
  if(Teuchos::is_null(oldDescentDirPtr))
    oldDescentDirPtr = soln.getX().clone(NOX::ShapeCopy);
  // These are conditionally created
  if(Teuchos::is_null(diffVecPtr) && usePRbeta)
    diffVecPtr = soln.getX().clone(NOX::ShapeCopy);
  if(Teuchos::is_null(tmpVecPtr) && doPrecondition)
    tmpVecPtr = soln.getX().clone(NOX::ShapeCopy);

  // Get a reference to the old solution group (const)
  oldSolnPtr = &solver.getPreviousSolutionGroup();
  const Abstract::Group& oldSoln(*oldSolnPtr);

  niter = solver.getNumIterations();

  // Construct Residual and precondition (if desired) as first step in
  // getting new search direction

  ok = soln.computeF();
  if (ok != Abstract::Group::Ok)
  {
    if (utils->isPrintType(Utils::Warning))
      utils->out() << "NOX::Direction::NonlinearCG::compute - Unable to compute F." << std::endl;
    return false;
  }

  dir = soln.getF();

  if(doPrecondition)
  {
    if(!soln.isJacobian())
      ok = soln.computeJacobian();
    if (ok != Abstract::Group::Ok)
    {
      if (utils->isPrintType(Utils::Warning))
        utils->out() << "NOX::Direction::NonlinearCG::compute - Unable to compute Jacobian." << std::endl;
      return false;
    }

    *tmpVecPtr = dir;

    ok = soln.applyRightPreconditioning(false, paramsPtr->sublist("Nonlinear CG").sublist("Linear Solver"), *tmpVecPtr, dir);
    if( ok != Abstract::Group::Ok )
    {
      if (utils->isPrintType(Utils::Warning))
        utils->out() << "NOX::Direction::NonlinearCG::compute - Unable to apply Right Preconditioner." << std::endl;
      return false;
    }
  }

  dir.scale(-1.0);

  // Orthogonalize using previous search direction

  beta = 0.0;

  if( niter!=0 )
  {
    // Two choices (for now) for orthogonalizing descent direction with previous:
    if( usePRbeta )
    {
      // Polak-Ribiere beta
      *diffVecPtr = dir;
      diffVecPtr->update(-1.0, *oldDescentDirPtr, 1.0);

      double denominator = oldDescentDirPtr->innerProduct(oldSoln.getF());

      beta = diffVecPtr->innerProduct(soln.getF()) / denominator;

      // Constrain beta >= 0
      if( beta < 0.0 )
      {
        if (utils->isPrintType(Utils::OuterIteration))
          utils->out() << "BETA < 0, (" << beta << ") --> Resetting to zero" << std::endl;
        beta = 0.0;
      }
    }
    else
    {
      // Fletcher-Reeves beta
      double denominator = oldDescentDirPtr->innerProduct(oldSoln.getF());

      beta = dir.innerProduct(soln.getF()) / denominator;

    }

    //  Allow for restart after specified number of nonlinear iterations
    if( (niter % restartFrequency) == 0 )
    {
      if( utils->isPrintType(Utils::OuterIteration) )
        utils->out() << "Resetting beta --> 0" << std::endl;

      beta = 0 ;  // Restart with Steepest Descent direction
    }
  } // niter != 0

  *oldDescentDirPtr = dir;

  dir.update(beta, *oldDirPtr, 1.0);

  *oldDirPtr = dir;

  return (ok == Abstract::Group::Ok);
}

bool NonlinearCG::compute(Abstract::Vector& dir, Abstract::Group& soln,
                          const Solver::LineSearchBased& solver)
{
  return NOX::Direction::Generic::compute( dir, soln, solver );
}
