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
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
// 
// ************************************************************************
//@HEADER

#include "NOX_Direction_NonlinearCG.H"	// class definition
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Common.H"
#include "NOX_Parameter_List.H"
#include "NOX_Solver_Generic.H"
#include "NOX_Utils.H"

/* Some compilers (in particular the SGI and ASCI Red - TFLOP)
 * fail to find the max and min function.  Therfore we redefine them
 * here.
 */
#ifdef max
#undef max
#endif

#define max(a,b) ((a)>(b)) ? (a) : (b);

#ifdef min
#undef min
#endif

#define min(a,b) ((a)<(b)) ? (a) : (b);

using namespace NOX;
using namespace NOX::Direction;

NonlinearCG::NonlinearCG(Parameter::List& params) :
  oldSolnPtr(NULL),			// pointer to old Soln Grp
  tmpVecPtr(NULL),			// reference to xgrp
  oldDirPtr(NULL),			// reference to xgrp
  oldDescentDirPtr(NULL),		// reference to xgrp
  diffVecPtr(NULL),			// reference to xgrp
  iparams(params),			// copy p
  restartFrequency(params.getParameter("Restart Frequency", 10))
{
  reset(iparams);
}


bool NonlinearCG::reset(Parameter::List& params) 
{
  return true;
}

NonlinearCG::~NonlinearCG() 
{
  delete tmpVecPtr;
  delete diffVecPtr;
  delete oldDescentDirPtr;
  delete oldDirPtr;
}


bool NonlinearCG::compute(Abstract::Vector& dir, Abstract::Group& soln,
                          const Solver::Generic& solver)
{
  bool status = false;

  // Initialize vector memory if haven't already
  if(oldDirPtr==NULL)
    oldDirPtr = soln.getX().clone(NOX::ShapeCopy);
  if(oldDescentDirPtr==NULL)
    oldDescentDirPtr = soln.getX().clone(NOX::ShapeCopy);
  if(diffVecPtr==NULL)
    diffVecPtr = soln.getX().clone(NOX::ShapeCopy);
  if(tmpVecPtr==NULL)
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
  /* NOTE FROM TAMMY: Need to check the return status! */

  soln.computeF();
  dir = soln.getF();  
  if(iparams.isParameterEqual("Precondition", "On")) {
    if(!soln.isJacobian())
      soln.computeJacobian();
    tmpVec = dir;
    status = soln.applyJacobianDiagonalInverse(tmpVec, dir);
  }

  dir.scale(-1.0);


  // Orthogonalize using previous search direction

  beta = 0.0;

  if(niter!=0){  

// Two choices (for now) for orthogonalizing descent direction with previous:

    if(iparams.isParameterEqual("Orthogonalize", "Polak-Ribiere"))
    {
//                     Polak-Ribiere beta

      diffVec = dir;
      diffVec.update(-1.0, oldDescentDir, 1.0); 

      double denominator = oldDescentDir.dot(oldSoln.getF());

      beta = diffVec.dot(soln.getF()) / denominator;

    // Constrain beta >= 0
      if(beta<0.0) {
        if (Utils::doPrint(Utils::OuterIteration))
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
       if (Utils::doPrint(Utils::OuterIteration))
         cout << "Resetting beta --> 0" << endl;

       beta = 0 ;  // Restart with Steepest Descent direction
    }

  } // niter != 0

  oldDescentDir = dir;

  dir.update(beta, oldDir, 1.0);

  oldDir = dir;

  status = true;

  return status;
}

