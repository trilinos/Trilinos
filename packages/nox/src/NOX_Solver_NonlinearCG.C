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

#include "NOX_Solver_NonlinearCG.H"	// class definition

using namespace NOX;
using namespace NOX::Solver;

NonlinearCG::NonlinearCG(Abstract::Group& xgrp, Status::Test& t, const Parameter::List& p) :
  solnptr(&xgrp),		// reference to xgrp
  oldsolnptr(xgrp.clone(DeepCopy)), // create via clone
  oldsoln(*oldsolnptr),		// reference to just-created pointer
  dirptr(xgrp.getX().clone(CopyShape)), // create via clone 
  dir(*dirptr),			// reference to just-created pointer
  olddirptr(xgrp.getX().clone(CopyShape)), // create via clone 
  olddir(*olddirptr),		// reference to just-created pointer
  precdirptr(xgrp.getX().clone(CopyShape)), // create via clone 
  precolddirptr(xgrp.getX().clone(CopyShape)), // create via clone 
  diffVector(xgrp.getX().clone(CopyShape)), // create via clone 
  testptr(&t),			// reference to t
  iparams(p),			// copy p
  oparams(),			// empty list
  linesearch(iparams.sublist("Line Search")), // initialize line search
  step(0.0),			// initialize to zero
  niter(0),			// initialize to zero
  restartFrequency(p.getParameter("Restart Frequency", 100)),
  outputFrequency(p.getParameter("Output Frequency", 1)),
				// initialize local variables to minimize 
				// Parameter::List access
  status(Status::Unconverged)	// initialize convergence status
{
  init();
}

// Protected
void NonlinearCG::init()
{
  // Print out initialization information
  if (Utils::doPrint(Utils::Parameters)) {
    cout << "\n" << Utils::fill(72) << "\n";
    cout << "\n-- Parameters Passed to Nonlinear Solver --\n\n";
    iparams.print(cout,5);
    cout << "\n-- Status Tests Passed to Nonlinear Solver --\n\n";
    testptr->print(cout, 5);
    cout <<"\n" << Utils::fill(72) << "\n";

  }

  // Compute RHS of initital guess
  solnptr->computeRHS();
}

NonlinearCG::~NonlinearCG() 
{
  delete oldsolnptr;
  delete dirptr;
  delete olddirptr;
  delete precdirptr;
  delete precolddirptr;
  delete diffVector;
  delete testptr;
}

bool NonlinearCG::reset(Abstract::Group& xgrp, Status::Test& t, const Parameter::List& p) 
{
  solnptr = &xgrp;
  testptr = &t;
  iparams = p;
  linesearch.reset(iparams.sublist("Line Search"));
  niter = 0;
  status = Status::Unconverged;
  init();
  return true;
}

Status::StatusType NonlinearCG::getStatus()
{
  status = testptr->operator()(*this);
  return status;
}

Status::StatusType NonlinearCG::iterate()
{
  // Copy pointers into temporary references
  Abstract::Group& soln = *solnptr;
  Status::Test& test = *testptr;

  // Construct Residual as first step in getting new search direction
  soln.computeRHS();  

  // Compute NonlinearCG direction for current solution.
  /* NOTE FROM TAMMY: Need to check the return status! */
  //  Two choices available for determining initial descent direction before
  //  orthogonalization: 
  if(iparams.isParameterEqual("Direction", "Richardson"))
  {
    dir = soln.getRHS();  // Richardson direction
    if(niter!=0) 
      olddescentdirptr = &oldsoln.getRHS();
  }
  else
  {
    soln.computeJacobian();
    soln.computeGrad(); 
    dir = soln.getGrad(); // Steepest Descent direction for 
                          // f = 1/2 Trans(R).R
    if(niter!=0) 
      olddescentdirptr = &oldsoln.getGrad();
  }
  dir.scale(-1.0);

  // Diagonally precondition if desired

  *precdirptr = dir;
  if(iparams.isParameterEqual("Diagonal Precondition", "On")) {
    if(!soln.isJacobian())
      soln.computeJacobian();
    soln.applyJacobianDiagonalInverse(dir, *precdirptr);
  }

  // Orthogonalize using previous search direction

  if(niter!=0){  
    *precolddirptr = *olddescentdirptr;
    if(iparams.isParameterEqual("Diagonal Precondition", "On")) 
      soln.applyJacobianDiagonalInverse(*olddescentdirptr, 
                                      *precolddirptr);

// Two choices (for now) for orthogonalizing descent direction with previous:

    if(iparams.isParameterEqual("Orthogonalize", "Polak-Ribiere"))
    {
//                     Polak-Ribiere beta

      *diffVector = *precdirptr;
      diffVector->update(1.0, *precolddirptr, 1.0); 

      double denominator = olddescentdirptr->dot(*precolddirptr);

      beta = dir.dot(*diffVector) / denominator;

    // Constrain beta >= 0
      if(beta<0.0) {
         cout << "BETA < 0, (" << beta << ") --> Resetting to zero" << endl;
         beta = 0.0;
      }
    } 
    else
    {
//                     Fletcher-Reeves beta

      double denominator = olddescentdirptr->dot(*precolddirptr);

      beta = dir.dot(*precdirptr) / denominator;

    } // End of orthogonalization


//  Allow for restart after specified number of nonlinear iterations

    if( (niter % restartFrequency)==0)
    {
       if (Utils::doPrint(Utils::OuterIteration))
         cout << "Resetting beta --> 0" << endl;

       beta = 0 ;  // Restart with Steepest Descent direction
    }

    precdirptr->update(beta, olddir, 1.0);

  } // niter != 0


  // Store direction vector for use in orthogonalization
  dir = *precdirptr;
  olddir = dir; 

  // Copy current soln to the old soln.
  oldsoln = soln;

//  Debugging,  RH
//  Check to see if dir is a descent direction and what dir it points to...
//  if(!oldsoln.isGrad())
//    oldsoln.computeGrad();
//  double testInitialSlope = dir.dot(oldsoln.getGrad());
//  cout << "  Inside NonlinearCG, dginit :" << testInitialSlope << endl;
//  cin.get();
  
  
  // Do line search and compute new soln.
  /* NOTE FROM TAMMY: Need to check the return status! */
  linesearch(soln, step, oldsoln, dir); // niter needs to be added, RH

  // Compute RHS for new current solution.
  soln.computeRHS();

  // Update iteration count.
  niter ++;

  // Evaluate the current status.
  status = test(*this);

  // Return status.
  return status;
}

Status::StatusType NonlinearCG::solve()
{
  status = testptr->operator()(*this);
  printUpdate();

  // Iterate until converged or failed
  while (status == Status::Unconverged) {
    status = iterate();
    if((niter % outputFrequency)==0)
      printUpdate();
  }

  return status;
}

const Abstract::Group& NonlinearCG::getSolutionGroup() const
{
  return *solnptr;
}

const Abstract::Group& NonlinearCG::getPreviousSolutionGroup() const
{
  return oldsoln;
}

int NonlinearCG::getNumIterations() const
{
  return niter;
}

const Parameter::List& NonlinearCG::getOutputParameters() const
{
  oparams.setParameter("Nonlinear Iterations", niter);
  oparams.setParameter("2-Norm of Residual", solnptr->getNormRHS());
  return oparams;
}

// protected
void NonlinearCG::printUpdate() 
{
  double norm_k;
  double norm_update;

  // All processors participate in the computation of these norms...
  if (Utils::doAllPrint(Utils::OuterIteration)) {
    norm_k = solnptr->getNormRHS();
    norm_update = (niter > 0) ? olddirptr->norm() : 0; 
  }

  // ...But only the print processors actually prints the result.
  if (Utils::doPrint(Utils::OuterIteration)) {
    cout << "\n" << Utils::fill(72) << "\n";
    cout << "-- NonlinearCG Step " << niter << " -- \n";
    cout << "f = " << Utils::sci(norm_k);
    cout << "  step = " << Utils::sci(step);
    cout << "  dx = " << Utils::sci(norm_update);
    if (status > 0)
      cout << " (Converged!)";
    if (status < 0)
      cout << " (Failed!)";
    cout << "\n" << Utils::fill(72) << "\n" << endl;
  }

  if ((status != 0) && (Utils::doPrint(Utils::OuterIteration))) {
    cout << Utils::fill(72) << "\n";
    cout << "-- Final Status Test Results --\n";
    testptr->print(cout);
    cout << Utils::fill(72) << "\n";
  }
}




