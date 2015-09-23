// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
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
