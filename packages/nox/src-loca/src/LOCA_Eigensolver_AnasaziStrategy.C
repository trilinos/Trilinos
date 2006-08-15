// $Id$
// $Source$

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
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

#include "NOX_Abstract_MultiVector.H"
#include "Teuchos_ParameterList.hpp"
#include "LOCA_Parameter_SublistParser.H"
#include "LOCA_GlobalData.H"
#include "NOX_Utils.H"
#include "LOCA_Factory.H"
#include "LOCA_Eigensolver_AnasaziStrategy.H"
#include "LOCA_EigenvalueSort_Strategies.H"

#ifdef HAVE_LOCA_ANASAZI
#include "Anasazi_LOCA_MultiVecTraits.H"
#include "Anasazi_LOCA_OperatorTraits.H"
#include "AnasaziBlockKrylovSchur.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#endif

LOCA::Eigensolver::AnasaziStrategy::AnasaziStrategy(
	const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
	const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams_,
	const Teuchos::RefCountPtr<Teuchos::ParameterList>& eigParams) :
  globalData(global_data),
  topParams(topParams_),
  eigenParams(eigParams),
  solverParams(),
  blksz(1),
  length(30),
  nev(4),
  tol(1.0e-7),
  step(1),
  restart(1),
  debug(1),
  which("LM"),
  saveEV(0),
  cayleyPole(0.0),
  cayleyZero(0.0),
  LOCA_PL(),
  LOCA_OM(),
  LOCASort()
{
  solverParams = topParams->getSublist("Linear Solver");

  // Get values out of parameter list
  blksz = eigenParams->get("Block Size", 1);
  length = eigenParams->get("Arnoldi Size", 30); 
  nev = eigenParams->get("NEV", 4);   
  tol = eigenParams->get("Tol", 1.0e-7);
  step = eigenParams->get("Convergence Check", 1);    
  restart = eigenParams->get("Restarts",1);
  debug = eigenParams->get("Debug Level",1);

  // Set default on which eigenvalues are of interest.  
  // Different defaults for Cayley and others.
  if (eigenParams->get("Operator","Jacobian Inverse") == "Cayley")
    which = eigenParams->get("Sorting Order","CA");   
  else  
    which = eigenParams->get("Sorting Order","LM");

  saveEV = eigenParams->get("Save Eigenvectors", 0);
  cayleyPole = eigenParams->get("Cayley Pole",0.0);
  cayleyZero = eigenParams->get("Cayley Zero",0.0);

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
  Teuchos::RefCountPtr<LOCA::EigenvalueSort::AbstractStrategy> sortingStrategy
    = globalData->locaFactory->createEigenvalueSortStrategy(topParams,
							    eigenParams);
  LOCASort =
    Teuchos::rcp(new Anasazi::LOCASort(sortingStrategy));
}

LOCA::Eigensolver::AnasaziStrategy::~AnasaziStrategy() 
{
}

NOX::Abstract::Group::ReturnType
LOCA::Eigensolver::AnasaziStrategy::computeEigenvalues(
		 NOX::Abstract::Group& group,
		 Teuchos::RefCountPtr< std::vector<double> >& evals_r,
		 Teuchos::RefCountPtr< std::vector<double> >& evals_i,
		 Teuchos::RefCountPtr< NOX::Abstract::MultiVector >& evecs_r,
		 Teuchos::RefCountPtr< NOX::Abstract::MultiVector >& evecs_i)
{
  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperIteration)) {
    globalData->locaUtils->out() << "\n" << 
      globalData->locaUtils->fill(64,'=') << 
      "\nAnasazi Eigensolver starting with block size " << blksz << 
      "\n" << std::endl;
  }

  // Get reference to solution vector to clone
  const NOX::Abstract::Vector& xVector = group.getX();

  // Create the operator and initial vector
  Teuchos::RefCountPtr<LOCA::AnasaziOperator::AbstractStrategy> anasaziOp
    = globalData->locaFactory->createAnasaziOperatorStrategy(
						   topParams, 
						   eigenParams,
						   solverParams,
						   Teuchos::rcp(&group,false));
  Teuchos::RefCountPtr<MV> ivec = xVector.createMultiVector(blksz);
  ivec->random();

  // Create an instance of the eigenproblem
  Teuchos::RefCountPtr<Anasazi::BasicEigenproblem<double, MV, OP> > 
    LOCAProblem =
    Teuchos::rcp( new Anasazi::BasicEigenproblem<double, MV, OP>(anasaziOp, 
								 ivec) );

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
  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperIteration)) 
    if (debug == 0) 
      LOCABlockKrylovSchur.currentStatus();
  
  // Obtain the eigenvalues / eigenvectors
  // If the matrix is non-symmetric the real part of the eigenvalues are
  // stored in the first "narn" entries of "evals", the imaginary part is 
  // stored in the second "narn" entries of "evals".
  int narn = LOCABlockKrylovSchur.GetKrylovFactLength();
  Teuchos::RefCountPtr< const std::vector<double> > evals = 
    LOCABlockKrylovSchur.GetRitzValues();
  // Copy first narn values into evals_r
  evals_r = 
    Teuchos::rcp(new std::vector<double>(evals->begin(),evals->begin()+narn));
  // Copy second narn values in evals_i
  evals_i = 
    Teuchos::rcp(new std::vector<double>(evals->begin()+narn, evals->end()));

  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperIteration)) 
     globalData->locaUtils->out() << 
       "Untransformed eigenvalues (since the operator was " << 
       anasaziOp->label() << ")" << std::endl;
  
  // Obtain the eigenvectors
  // The real part is stored in the first "nev" vectors and the imaginary 
  // in the second "nev" vectors.
  Teuchos::RefCountPtr<MV> evecs = LOCAProblem->GetEvecs();
  std::vector<int> index_r(nev);
  std::vector<int> index_i(nev);
  for (int i=0; i<nev; i++) {
    index_r[i] = i;
    index_i[i] = nev+i;
  }
  evecs_r = evecs->subCopy(index_r);
  evecs_i = evecs->subCopy(index_i);

  // Real & imaginary components of Rayleigh quotient
  double rq_r, rq_i;

  for (int i=0; i<nev; i++) {

    // Un-transform eigenvalues
    anasaziOp->transformEigenvalue((*evals_r)[i], (*evals_i)[i]);

    // Compute Rayleigh quotient
    anasaziOp->rayleighQuotient((*evecs_r)[i], (*evecs_i)[i], rq_r, rq_i);

    // Print out untransformed eigenvalues and Rayleigh quotient residual
    if (globalData->locaUtils->isPrintType(NOX::Utils::StepperIteration)) {
       globalData->locaUtils->out() << "Eigenvalue " << i << " : " << 
	 globalData->locaUtils->sciformat((*evals_r)[i]) << "  " << 
	 globalData->locaUtils->sciformat((*evals_i)[i]) << 
	 " i    :  RQresid " << 
	 globalData->locaUtils->sciformat(fabs((*evals_r)[i] - rq_r)) << 
	 "  " << 
	 globalData->locaUtils->sciformat(fabs((*evals_i)[i] - rq_i)) << 
	 " i" << std::endl;
    }

  }  

  // Print out remaining eigenvalue approximations from nev to 
  // final arnoldi size
  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperIteration) && 
      narn>nev) {
    globalData->locaUtils->out() << 
      "~~~~~~~ remaining eigenvalue approximations ~~~~~~~~~~~~" << std::endl;
  }
  for (int i=nev; i<narn; i++) {

      // Un-transform eigenvalues
    anasaziOp->transformEigenvalue((*evals_r)[i], (*evals_i)[i]);

    if (globalData->locaUtils->isPrintType(NOX::Utils::StepperIteration) && 
	narn>nev) {
	globalData->locaUtils->out() << 
	  "Eigenvalue " << i << " : " << 
	  globalData->locaUtils->sciformat((*evals_r)[i]) << "  " << 
	  globalData->locaUtils->sciformat((*evals_i)[i]) << " i" << std::endl;
      }

  }

  // Save eigenvectors/eigenvalues
//   if (saveEV > 0)
//     saveEigenVectors(saveEV, evals, evecs);

  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperIteration)) {
    globalData->locaUtils->out() << 
      "\nAnasazi Eigensolver finished.\n" << 
      globalData->locaUtils->fill(64,'=') << "\n" << std::endl;
  }

  return NOX::Abstract::Group::Ok;
}

