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
#include "LOCA_ErrorCheck.H"
#include "NOX_Utils.H"
#include "LOCA_Factory.H"
#include "LOCA_Eigensolver_AnasaziStrategy.H"
#include "LOCA_EigenvalueSort_Strategies.H"

#ifdef HAVE_LOCA_ANASAZI
#include "Anasazi_LOCA_MultiVecTraits.H"
#include "Anasazi_LOCA_OperatorTraits.H"
#include "AnasaziBlockKrylovSchurSolMgr.hpp"
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
  nev(4),
  isSymmetric(false),
  locaSort()
{
  solverParams = topParams->getSublist("Linear Solver");

  // Get values out of parameter list
  blksz = eigenParams->get("Block Size", 1);
  nev = eigenParams->get("Num Eigenvalues", 4);
  isSymmetric = eigenParams->get("Symmetric", false);

  // Set more reasonable defaults
  eigenParams->get("Convergence Tolerance", 1.0e-7);
  eigenParams->get("Maximum Restarts", 1);
  eigenParams->get("Num Blocks", 30);
  eigenParams->get("Step Size", 1);
  if (!eigenParams->isParameter("Verbosity"))
    eigenParams->set("Verbosity",  
		     Anasazi::Errors + 
		     Anasazi::Warnings +
		     Anasazi::FinalSummary);
		   
  // Create a sorting manager to handle the sorting of eigenvalues 
  Teuchos::RefCountPtr<LOCA::EigenvalueSort::AbstractStrategy> sortingStrategy
    = globalData->locaFactory->createEigenvalueSortStrategy(topParams,
							    eigenParams);
  locaSort =
    Teuchos::rcp(new Anasazi::LOCASort(globalData, sortingStrategy));
  eigenParams->set( "Sort Manager", locaSort );
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
  LOCAProblem->setNEV( nev );

  // Set symmetry
  LOCAProblem->setHermitian(isSymmetric);
 
  // Inform the eigenproblem that you are finishing passing it information
  //assert( LOCAProblem->setProblem() == 0 );
  LOCAProblem->setProblem();

  // Initialize the solver
  Anasazi::BlockKrylovSchurSolMgr<double, MV, OP> 
    LOCABlockKrylovSchur(LOCAProblem, *eigenParams); // Need to pass in sorter

  // Solve the problem to the specified tolerance
  Anasazi::ReturnType returnCode = LOCABlockKrylovSchur.solve();
  
  // Obtain the eigenvalues / eigenvectors
  const Anasazi::Eigensolution<double,MV>& anasaziSolution = 
    LOCAProblem->getSolution();
  evals_r = 
    Teuchos::rcp(new std::vector<double>(nev));
  evals_i = 
    Teuchos::rcp(new std::vector<double>(nev));
  for (int i=0; i<nev; i++) {
    (*evals_r)[i] = anasaziSolution.Evals[i].realpart;
    (*evals_i)[i] = anasaziSolution.Evals[i].imagpart;
  }

  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperIteration)) 
     globalData->locaUtils->out() << 
       "Untransformed eigenvalues (since the operator was " << 
       anasaziOp->label() << ")" << std::endl;
  
  // Obtain the eigenvectors
  Teuchos::RefCountPtr<MV> evecs = anasaziSolution.Evecs;
  evecs_r = evecs->clone(nev);
  evecs_i = evecs->clone(nev);
  for (int i=0; i<nev; i++) {

    // Eigenvalue is real
    if (anasaziSolution.index[i] == 0) {
      (*evecs_r)[i] = (*evecs)[i];
      (*evecs_i)[i].init(0.0);
    }
    
    // Complex conjugate pair.  We must have i<nev-1 for this to be true
    else if (anasaziSolution.index[i] == 1) {
      (*evecs_r)[i] = (*evecs)[i];
      (*evecs_i)[i] = (*evecs)[i+1];
    }

    // Previous complex conjugate pair.  We must have i>0 for this to be true
    // Take conjugate of imaginary part
    else if (anasaziSolution.index[i] == -1) {
      (*evecs_r)[i] = (*evecs)[i-1];
      (*evecs_i)[i].update(-1.0, (*evecs)[i], 0.0);
    }
    else {
      string func = "LOCA::Eigensolver::AnasaziStrategy::computeEigenvalues()";
      globalData->locaErrorCheck->throwError(func, "Unknown anasazi index");
    }
  }

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

  // Print out remaining eigenvalue approximations
  std::vector<Anasazi::Value<double> > ritzValues = 
    LOCABlockKrylovSchur.getRitzValues();
  int numRitz = ritzValues.size();
  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperIteration) && 
      numRitz>nev) {
    globalData->locaUtils->out() << 
      "~~~~~~~ remaining eigenvalue approximations ~~~~~~~~~~~~" << std::endl;
  }
  for (int i=nev; i<numRitz; i++) {

      // Un-transform eigenvalues
    anasaziOp->transformEigenvalue(ritzValues[i].realpart, 
				   ritzValues[i].imagpart);

    if (globalData->locaUtils->isPrintType(NOX::Utils::StepperIteration)) {
	globalData->locaUtils->out() << 
	  "Eigenvalue " << i << " : " << 
	  globalData->locaUtils->sciformat(ritzValues[i].realpart) << "  " << 
	  globalData->locaUtils->sciformat(ritzValues[i].imagpart) << " i" <<
	  std::endl;
      }

  }

  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperIteration)) {
    globalData->locaUtils->out() << 
      "\nAnasazi Eigensolver finished.\n" << 
      globalData->locaUtils->fill(64,'=') << "\n" << std::endl;
  }

  if (returnCode == Anasazi::Converged)
    return NOX::Abstract::Group::Ok;
  else
    return NOX::Abstract::Group::NotConverged;
}

