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
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//                                                                                 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA                                                                                
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#ifdef WITH_PRERELEASE

#include "NOX_Direction_NonlinearCG.H"	// class definition
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Common.H"
#include "NOX_Parameter_List.H"
#include "NOX_Solver_Generic.H"
#include "NOX_Utils.H"

using namespace NOX;
using namespace NOX::Direction;

NonlinearCG::NonlinearCG(const NOX::Utils& u, Parameter::List& params) :
  utils(u),
  oldSolnPtr(0),			// pointer to old Soln Grp
  tmpVecPtr(0),				// reference to xgrp
  oldDirPtr(0),				// reference to xgrp
  oldDescentDirPtr(0),			// reference to xgrp
  diffVecPtr(0),			// reference to xgrp
  paramsPtr(0)
{
  reset(params);
}


bool NonlinearCG::reset(Parameter::List& params) 
{
  paramsPtr = &params;
  NOX::Parameter::List& nlcgParams = paramsPtr->sublist("Nonlinear CG");
  restartFrequency = nlcgParams.getParameter("Restart Frequency", 10);
  doPrecondition = nlcgParams.isParameterEqual("Precondition", "On");
  usePRbeta = nlcgParams.isParameterEqual("Orthogonalize", "Polak-Ribiere");
  return true;
}

NonlinearCG::~NonlinearCG() 
{
  delete tmpVecPtr; tmpVecPtr = 0;
  delete diffVecPtr; diffVecPtr = 0;
  delete oldDescentDirPtr; oldDescentDirPtr = 0;
  delete oldDirPtr; oldDirPtr = 0;
}


bool NonlinearCG::compute(Abstract::Vector& dir, Abstract::Group& soln,
                          const Solver::Generic& solver)
{
  Abstract::Group::ReturnType ok;

  // Initialize vector memory if haven't already
  if(oldDirPtr==0)
    oldDirPtr = soln.getX().clone(NOX::ShapeCopy);
  if(oldDescentDirPtr==0)
    oldDescentDirPtr = soln.getX().clone(NOX::ShapeCopy);
  // These are conditionally created
  if(diffVecPtr==0 && usePRbeta)
    diffVecPtr = soln.getX().clone(NOX::ShapeCopy);
  if(tmpVecPtr==0 && doPrecondition)
    tmpVecPtr = soln.getX().clone(NOX::ShapeCopy);

  // Create references to vectors for convenience
  Abstract::Vector& oldDir(*oldDirPtr);
  Abstract::Vector& oldDescentDir(*oldDescentDirPtr);
  Abstract::Vector& diffVec(*diffVecPtr);
  Abstract::Vector& tmpVec(*tmpVecPtr);

  // Get a reference to the old solution group (const)
  oldSolnPtr = &solver.getPreviousSolutionGroup();
  const Abstract::Group& oldSoln(*oldSolnPtr);

  niter = solver.getNumIterations();

  // Construct Residual and precondition (if desired) as first step in 
  // getting new search direction

  ok = soln.computeF();
  if (ok != Abstract::Group::Ok) {
    if (utils.isPrintProcessAndType(Utils::Warning))
      cout << "NOX::Direction::NonlinearCG::compute - Unable to compute F." 
           << endl;
    return false;
  }
  dir = soln.getF();  
  if(paramsPtr->sublist("Nonlinear CG").isParameterEqual("Precondition", "On")) {
    if(!soln.isJacobian())
      ok = soln.computeJacobian();
      if (ok != Abstract::Group::Ok) {
        if (utils.isPrintProcessAndType(Utils::Warning))
          cout << "NOX::Direction::NonlinearCG::compute - "
               << "Unable to compute Jacobian." << endl;
        return false;
      }
    tmpVec = dir;
    ok = soln.applyRightPreconditioning(false, paramsPtr->sublist("Nonlinear CG").sublist("Linear Solver"), tmpVec, dir);
    if (ok != Abstract::Group::Ok) {
      if (utils.isPrintProcessAndType(Utils::Warning))
        cout << "NOX::Direction::NonlinearCG::compute - "
             << "Unable to apply Right Preconditioner." << endl;
      return false;
    }
  }

  dir.scale(-1.0);


  // Orthogonalize using previous search direction

  beta = 0.0;

  if(niter!=0){  

// Two choices (for now) for orthogonalizing descent direction with previous:

    if(paramsPtr->sublist("Nonlinear CG").isParameterEqual("Orthogonalize", "Polak-Ribiere"))
    {
//                     Polak-Ribiere beta

      diffVec = dir;
      diffVec.update(-1.0, oldDescentDir, 1.0); 

      double denominator = oldDescentDir.dot(oldSoln.getF());

      beta = diffVec.dot(soln.getF()) / denominator;

    // Constrain beta >= 0
      if(beta < 0.0) {
        if (utils.isPrintProcessAndType(Utils::OuterIteration))
          cout << "BETA < 0, (" << beta << ") --> Resetting to zero" << endl;
        beta = 0.0;
      }
    } 
    else
    {
//                     Fletcher-Reeves beta

      double denominator = oldDescentDir.dot(oldSoln.getF());

      beta = dir.dot(soln.getF()) / denominator;

    } // End of orthogonalization


//  Allow for restart after specified number of nonlinear iterations

    if( (niter % restartFrequency)==0)
    {
       if (utils.isPrintProcessAndType(Utils::OuterIteration))
         cout << "Resetting beta --> 0" << endl;

       beta = 0 ;  // Restart with Steepest Descent direction
    }

  } // niter != 0

  oldDescentDir = dir;

  dir.update(beta, oldDir, 1.0);

  oldDir = dir;

  return (ok == Abstract::Group::Ok);
}

#endif
