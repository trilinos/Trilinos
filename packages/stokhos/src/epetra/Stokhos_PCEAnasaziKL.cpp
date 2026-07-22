// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Stokhos_PCEAnasaziKL.hpp"
#ifdef HAVE_STOKHOS_ANASAZI

#include "AnasaziBlockKrylovSchurSolMgr.hpp"
#include "AnasaziBasicSort.hpp"

Stokhos::PCEAnasaziKL::
PCEAnasaziKL(const Stokhos::VectorOrthogPoly<Epetra_Vector>& X_poly,
	     int num_KL_) :
  covOp(Teuchos::rcp(new Stokhos::PCECovarianceOp(X_poly))),
  num_KL(num_KL_)
{ 
}

Stokhos::PCEAnasaziKL::
PCEAnasaziKL(const Teuchos::RCP<const EpetraExt::BlockVector>& X,
	     const Stokhos::OrthogPolyBasis<int,double>& basis,
	     int num_KL_) :
  covOp(Teuchos::rcp(new Stokhos::PCECovarianceOp(X, basis))),
  num_KL(num_KL_)
{
}

Stokhos::PCEAnasaziKL::
PCEAnasaziKL(const Teuchos::RCP<const Epetra_MultiVector>& X,
	     const Stokhos::OrthogPolyBasis<int,double>& basis,
	     int num_KL_) :
  covOp(Teuchos::rcp(new Stokhos::PCECovarianceOp(X, basis))),
  num_KL(num_KL_)
{
}

Teuchos::ParameterList
Stokhos::PCEAnasaziKL::
getDefaultParams() const
{
  Teuchos::ParameterList params;

  params.set("Verbosity",  
	     Anasazi::FinalSummary + 
	     //Anasazi::TimingDetails + 
	     Anasazi::Errors + 
	     Anasazi::Warnings);
  params.set("Which", "LM");  
  params.set("Block Size", 1);
  params.set("Num Blocks", 3*num_KL);
  params.set("Step Size", 5);
  params.set("Maximum Restarts", 1);
  params.set("Convergence Tolerance", 1e-12);

  return params;
}

bool
Stokhos::PCEAnasaziKL::
computeKL(Teuchos::ParameterList& anasazi_params)
{
  // Create an Epetra_MultiVector for an initial vector to start the solver.
  // Note:  This needs to have the same number of columns as the blocksize.
  Teuchos::RCP<Epetra_MultiVector> ivec = 
    Teuchos::rcp(new Epetra_MultiVector(covOp->CoeffMap(), 
					anasazi_params.get<int>("Block Size")));
  ivec->SetSeed(1);
  ivec->Random();
  
  // Create the eigenproblem.
  anasazi_problem =
    Teuchos::rcp(new Anasazi::BasicEigenproblem<ScalarType,MV,OP>(covOp, ivec));
  
  // Inform the eigenproblem that the operator A is symmetric
  anasazi_problem->setHermitian(true); 
  
  // Set the number of eigenvalues requested
  anasazi_problem->setNEV(num_KL);
  
  // Inform the eigenproblem that you are finishing passing it information
  anasazi_problem->setProblem();
  
  // Initialize the Block Arnoldi solver
  Anasazi::BlockKrylovSchurSolMgr<ScalarType,MV,OP> solverMgr(anasazi_problem, 
							      anasazi_params);
  
  // Solve the problem to the specified tolerances or length
  Anasazi::ReturnType returnCode = solverMgr.solve();

  // Check convergence
  bool result = true;
  if (returnCode != Anasazi::Converged) {
    result = false;
  }

  // Get solution
  sol = anasazi_problem->getSolution();

  return result;
}

Teuchos::Array<double>
Stokhos::PCEAnasaziKL::
getEigenvalues() const
{
  Teuchos::Array<double> evals(num_KL);
  for (int i=0; i<num_KL; i++)
    evals[i] = std::abs(sol.Evals[i].realpart);
  return evals;
}

Teuchos::RCP<Epetra_MultiVector>
Stokhos::PCEAnasaziKL::
getEigenvectors() const
{
  return sol.Evecs;
}

#endif // HAVE_STOKHOS_ANASAZI
