// $Id$
// $Source$

//@HEADER
// ************************************************************************
//
//                  LOCA Continuation Algorithm Package
//                 Copyright (2005) Sandia Corporation
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
// Questions? Contact Andy Salinger (agsalin@sandia.gov) or Eric Phipps
// (etphipp@sandia.gov), Sandia National Laboratories.
//
// ************************************************************************
//@HEADER

#include "NOX_Parameter_List.H"
#include "LOCA_GlobalData.H"
#include "LOCA_Utils.H"
#include "LOCA_Eigensolver_AnasaziStrategy.H"

#ifdef HAVE_LOCA_ANASAZI
#include "Anasazi_LOCA_MultiVec.H"
#include "Anasazi_LOCA_Matrix.H"
#include "LOCA_AnasaziOperator_Manager.H"
#include "AnasaziBlockKrylovSchur.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#endif

LOCA::Eigensolver::AnasaziStrategy::AnasaziStrategy(
		 const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
		 const Teuchos::RefCountPtr<NOX::Parameter::List>& eigParams,
		 const Teuchos::RefCountPtr<NOX::Parameter::List>& solParams) :
  globalData(global_data)
{
#ifdef HAVE_LOCA_ANASAZI
  eigenParams = eigParams;
  solverParams = solParams;

  // Get values out of parameter list
  blksz = eigenParams->getParameter("Block Size", 1);
  length = eigenParams->getParameter("Arnoldi Size", 30); 
  nev = eigenParams->getParameter("NEV", 4);   
  tol = eigenParams->getParameter("Tol", 1.0e-7);
  step = eigenParams->getParameter("Convergence Check", 1);    
  restart = eigenParams->getParameter("Restarts",1);
  debug = eigenParams->getParameter("Debug Level",1);

  // Set default on which eigenvalues are of interest.  
  // Different defaults for Cayley and others.
  if (eigenParams->getParameter("Operator","Jacobian Inverse") == "Cayley")
    which = eigenParams->getParameter("Sorting Order","CA");   
  else  
    which = eigenParams->getParameter("Sorting Order","LM");

  saveEV = eigenParams->getParameter("Save Eigenvectors", 0);
  cayleyPole = eigenParams->getParameter("Cayley Pole",0.0);
  cayleyZero = eigenParams->getParameter("Cayley Zero",0.0);

  //  Make sure saveEV is an appropriate value
  if (saveEV > nev)
    saveEV = nev;

  // Create parameter list to pass into solver
  LOCA_PL.set( "Block Size", blksz );
  LOCA_PL.set( "Max Blocks", length );
  LOCA_PL.set( "Max Restarts", restart );
  LOCA_PL.set( "Step Size", step );
  LOCA_PL.set( "Tol", tol );

  // Create an output manager to handle the I/O from the solver
  LOCA_OM = Teuchos::rcp( new Anasazi::OutputManager<double>() );
  LOCA_OM->SetVerbosity( debug );  

  // Create a sorting manager to handle the sorting of eigenvalues 
  LOCASort =
    Teuchos::rcp(new Anasazi::LOCASort<double, MV, OP>(which, 
						       cayleyPole,
						       cayleyZero));
#endif
}

LOCA::Eigensolver::AnasaziStrategy::~AnasaziStrategy() 
{
}

NOX::Abstract::Group::ReturnType
LOCA::Eigensolver::AnasaziStrategy::computeEigenvalues(
						  NOX::Abstract::Group& group)
{
#ifdef HAVE_LOCA_ANASAZI
  if (globalData->locaUtils->doPrint(Utils::StepperIteration)) {
    cout << "\n" << globalData->locaUtils->fill(64,'=')
         << "\nAnasazi Eigensolver starting with block size " << blksz
	 << "\n" << endl;
  }

  // Get reference to solution vector to clone
  const NOX::Abstract::Vector& xVector = group.getX();

  // Create the operator and initial vector
  LOCA::AnasaziOperator::Manager anasaziOperator(*eigenParams, *solverParams, 
						 group);
  Teuchos::RefCountPtr<Anasazi::LOCAMatrix> Amat = 
    Teuchos::rcp( new Anasazi::LOCAMatrix(anasaziOperator) );
  Teuchos::RefCountPtr<Anasazi::LOCAMultiVec> ivec =
    Teuchos::rcp( new Anasazi::LOCAMultiVec(xVector, blksz) );
  ivec->MvRandom();

  // Create an instance of the eigenproblem
  Teuchos::RefCountPtr<Anasazi::BasicEigenproblem<double, MV, OP> > 
    LOCAProblem =
    Teuchos::rcp( new Anasazi::BasicEigenproblem<double, MV, OP>(Amat, ivec) );

  // Set the number of eigenvalues requested
  LOCAProblem->SetNEV( nev );
 
  // Inform the eigenproblem that you are finishing passing it information
  assert( LOCAProblem->SetProblem() == 0 );

  // Initialize the solver
  Anasazi::BlockKrylovSchur<double, MV, OP> 
    LOCABlockKrylovSchur(LOCAProblem, LOCASort, LOCA_OM, LOCA_PL);

  // Solve the problem to the specified tolerance
  LOCABlockKrylovSchur.solve();

  // Look at the solutions once if debug=0
  if (globalData->locaUtils->doPrint(Utils::StepperIteration)) 
    if (debug == 0) 
      LOCABlockKrylovSchur.currentStatus();
  
  // Obtain the eigenvalues / eigenvectors
  // If the matrix is non-symmetric the real part of the eigenvalues are
  // stored in the first "narn" entries of "evals", the imaginary part is 
  // stored in the second "narn" entries of "evals".
  int narn = LOCABlockKrylovSchur.GetKrylovFactLength();
  std::vector<double> evals( *(LOCABlockKrylovSchur.GetRitzValues()) );

  if (globalData->locaUtils->doPrint(Utils::StepperIteration)) 
    cout << "Untransformed eigenvalues (since the operator was " 
	 << anasaziOperator.label() << ")" <<endl;
  
  // Obtain the eigenvectors
  // The real part is stored in the first "nev" vectors and the imaginary 
  // in the second "nev" vectors.
  Anasazi::LOCAMultiVec 
    evecs( *dynamic_cast<Anasazi::LOCAMultiVec *>(LOCAProblem->GetEvecs().
						    get()) );

  // Real & imaginary components of Rayleigh quotient
  double rq_r, rq_i;

  for (int i=0; i<nev; i++) {

    // Un-transform eigenvalues
    anasaziOperator.transformEigenvalue(evals[i], evals[narn+i]);

    // Compute Rayleigh quotient
    anasaziOperator.rayleighQuotient(evecs.GetNOXVector(i),
				     evecs.GetNOXVector(nev+i), 
				     rq_r, rq_i);

    // Print out untransformed eigenvalues and Rayleigh quotient residual
    if (globalData->locaUtils->doPrint(Utils::StepperIteration)) {
      cout << "Eigenvalue " << i << " : " 
	   << globalData->locaUtils->sci(evals[i]) <<"  "
	   << globalData->locaUtils->sci(evals[narn+i]) << " i    :  RQresid "
	   << globalData->locaUtils->sci(fabs(evals[i] - rq_r)) << "  "
	   << globalData->locaUtils->sci(fabs(evals[narn+i] - rq_i)) 
	   << " i" << endl;
    }

  }  

  // Print out remaining eigenvalue approximations from nev to 
  // final arnoldi size
  if (globalData->locaUtils->doPrint(Utils::StepperIteration) && narn>nev) {
    cout << "~~~~~~~ remaining eigenvalue approximations ~~~~~~~~~~~~" 
	 << endl;
  }
  for (int i=nev; i<narn; i++) {

      // Un-transform eigenvalues
      anasaziOperator.transformEigenvalue(evals[i], evals[narn+i]);

      if (globalData->locaUtils->doPrint(Utils::StepperIteration) && 
	  narn>nev) {
	cout << "Eigenvalue " << i << " : " 
	     << globalData->locaUtils->sci(evals[i]) << "  "
	     << globalData->locaUtils->sci(evals[narn+i]) << " i" <<endl;
      }

  }

  // Save eigenvectors/eigenvalues
//   if (saveEV > 0)
//     saveEigenVectors(saveEV, evals, evecs);

  if (globalData->locaUtils->doPrint(Utils::StepperIteration)) {
    cout << "\nAnasazi Eigensolver finished.\n" 
         << globalData->locaUtils->fill(64,'=') << "\n" << endl;
  }
#else
  if (globalData->locaUtils->doPrint(Utils::StepperIteration)) {
    cout << endl 
	 << "Warning: LOCA::Eigensolver::AnasaziStrategy::computeEigenvalues:"
	 << endl
         << "Anasazi Eigensolver requested but not compiled in!" << endl;
  }
  return NOX::Abstract::Group::Ok;
#endif

  return NOX::Abstract::Group::Ok;
}
