// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright 2004 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER
//  This example computes the eigenvalues of smallest magnitude of the
//  discretized 2D Laplacian operator using the block Davidson method.
//  This problem shows the construction of a shifted eigenproblem that
//  targets the smallest eigenvalues around a certain value (sigma).
//  This operator is discretized using linear finite elements and constructed
//  as an Epetra matrix, then passed shifted using EpetraExt utilities.

// Include autoconfigured header
#include "AnasaziConfigDefs.hpp"

// Include header for block Davidson solver
#include "AnasaziBlockDavidsonSolMgr.hpp"

// Include header to define basic eigenproblem Ax = \lambda*Bx
#include "AnasaziBasicEigenproblem.hpp"

// Include header to provide Anasazi with Epetra adapters
#include "AnasaziEpetraAdapter.hpp"

// Include header to provide basic Anasazi output manager
#include "AnasaziBasicOutputManager.hpp"

// Include header for Epetra compressed-row storage matrix and linear problem
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"

// Include header for Teuchos serial dense matrix
#include "Teuchos_SerialDenseMatrix.hpp"

// Include header for the problem definition
#include "ModeLaplace2DQ2.h"

// Include EpetraExt MatrixMatrix helpers.
#include "EpetraExt_MatrixMatrix.h"

// Include selected communicator class and map required by Epetra objects
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"

// ****************************************************************************
// BEGIN MAIN ROUTINE
// ****************************************************************************

int main(int argc, char *argv[]) {
  int i;

#ifdef EPETRA_MPI
  // Initialize MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  int MyPID = Comm.MyPID();

  // Create an Anasazi output manager
  //
  Anasazi::BasicOutputManager<double> printer;
  printer.stream(Anasazi::Errors) << Anasazi::Anasazi_Version() << std::endl << std::endl;

  // Number of dimension of the domain
  int space_dim = 2;

  // Size of each of the dimensions of the domain
  std::vector<double> brick_dim( space_dim );
  brick_dim[0] = 1.0;
  brick_dim[1] = 1.0;

  // Number of elements in each of the dimensions of the domain
  std::vector<int> elements( space_dim );
  elements[0] = 10;
  elements[1] = 10;

  // Create problem
  Teuchos::RCP<ModalProblem> testCase = Teuchos::rcp( new ModeLaplace2DQ2(Comm, brick_dim[0], elements[0], brick_dim[1], elements[1]) );

  // Get the stiffness and mass matrices
  Teuchos::RCP<Epetra_CrsMatrix> K = Teuchos::rcp( const_cast<Epetra_CrsMatrix *>(testCase->getStiffness()), false );
  Teuchos::RCP<Epetra_CrsMatrix> M = Teuchos::rcp( const_cast<Epetra_CrsMatrix *>(testCase->getMass()), false );

  // Create the shifted system K - sigma * M.

  double sigma = 1.0;
  Teuchos::RCP<Epetra_CrsMatrix> Kshift = Teuchos::rcp( new Epetra_CrsMatrix( *K ) );

  int addErr = EpetraExt::MatrixMatrix::Add( *M, false, -sigma, *Kshift, 1.0 );
  if (addErr != 0) {
    printer.print(Anasazi::Errors,"EpetraExt::MatrixMatrix::Add returned with error.\n");
#ifdef HAVE_MPI
    MPI_Finalize() ;
#endif
    return -1;
  }

  //
  // ************************************
  // Start the block Davidson iteration
  // ************************************
  //
  //  Variables used for the Block Davidson Method
  //
  const int    nev         = 4;
  const int    blockSize   = 4;
  const int    numBlocks   = 5;
  const int    maxRestarts = 200;
  const double tol         = 1.0e-8;
  std::string which = "SM";
  int verbosity = Anasazi::Errors + Anasazi::Warnings + Anasazi::FinalSummary;
  //
  // Create parameter list to pass into solver
  //
  Teuchos::ParameterList MyPL;
  MyPL.set( "Verbosity", verbosity );
  MyPL.set( "Which", which );
  MyPL.set( "Block Size", blockSize );
  MyPL.set( "Num Blocks", numBlocks );
  MyPL.set( "Maximum Restarts", maxRestarts );
  MyPL.set( "Convergence Tolerance", tol );
  MyPL.set( "Use Locking", true );
  MyPL.set( "Locking Tolerance", tol/10 );

  typedef Epetra_MultiVector MV;
  typedef Epetra_Operator OP;
  typedef Anasazi::MultiVecTraits<double, MV> MVT;
  // typedef Anasazi::OperatorTraits<double, MV, OP> OPT; // unused

  // Create an Epetra_MultiVector for an initial vector to start the solver.
  // Note:  This needs to have the same number of columns as the blocksize.
  Teuchos::RCP<Epetra_MultiVector> ivec = Teuchos::rcp( new Epetra_MultiVector(K->Map(), blockSize) );
  MVT::MvRandom( *ivec );

  Teuchos::RCP<Anasazi::BasicEigenproblem<double,MV,OP> > MyProblem =
    Teuchos::rcp( new Anasazi::BasicEigenproblem<double,MV,OP>(Kshift, M, ivec) );

  // Inform the eigenproblem that the matrix pencil (K,M) is symmetric
  MyProblem->setHermitian(true);

  // Set the number of eigenvalues requested
  MyProblem->setNEV( nev );

  // Inform the eigenproblem that you are finished passing it information
  bool boolret = MyProblem->setProblem();
  if (boolret != true) {
    printer.print(Anasazi::Errors,"Anasazi::BasicEigenproblem::setProblem() returned an error.\n");
#ifdef HAVE_MPI
    MPI_Finalize() ;
#endif
    return -1;
  }

  // Initialize the Block Davidson solver
  Anasazi::BlockDavidsonSolMgr<double, MV, OP> MySolverMgr(MyProblem, MyPL);

  // Solve the problem to the specified tolerances or length
  Anasazi::ReturnType returnCode = MySolverMgr.solve();
  if (returnCode != Anasazi::Converged && MyPID==0) {
    std::cout << "Anasazi::EigensolverMgr::solve() returned unconverged." << std::endl;
  }

  // Get the eigenvalues and eigenvectors from the eigenproblem
  Anasazi::Eigensolution<double,MV> sol = MyProblem->getSolution();
  std::vector<Anasazi::Value<double> > evals = sol.Evals;
  Teuchos::RCP<MV> evecs = sol.Evecs;
  int numev = sol.numVecs;

  if (numev > 0) {

    // Undo shift transformation; computed eigenvalues are real
    std::vector<double> compEvals(numev);
    for (i=0; i<numev; ++i) {
      compEvals[i] = evals[i].realpart + sigma;
    }

    //************************************
    // Compute residuals, just for funsies
    //************************************
    //
    std::vector<double> normR(sol.numVecs);

    Teuchos::SerialDenseMatrix<int,double> T(sol.numVecs, sol.numVecs);
    Epetra_MultiVector Kvec( K->OperatorDomainMap(), evecs->NumVectors() );
    Epetra_MultiVector Mvec( M->OperatorDomainMap(), evecs->NumVectors() );
    T.putScalar(0.0);
    for (i=0; i<sol.numVecs; i++) {
      T(i,i) = compEvals[i];
    }
    K->Apply( *evecs, Kvec );
    M->Apply( *evecs, Mvec );
    MVT::MvTimesMatAddMv( -1.0, Mvec, T, 1.0, Kvec );
    MVT::MvNorm( Kvec, normR );

    //************************************
    // Print the results
    //************************************
    //
    std::ostringstream os;
    os.setf(std::ios_base::right, std::ios_base::adjustfield);
    os<<"Solver manager returned " << (returnCode == Anasazi::Converged ? "converged." : "unconverged.") << std::endl;
    os<<std::endl;
    os<<"------------------------------------------------------"<<std::endl;
    os<<std::setw(16)<<"Eigenvalue"
      <<std::setw(18)<<"Direct Residual"
      <<std::endl;
    os<<"------------------------------------------------------"<<std::endl;
    for (i=0; i<sol.numVecs; i++) {
      os<<std::setw(16)<<compEvals[i]
        <<std::setw(18)<<normR[i]/compEvals[i]
        <<std::endl;
    }
    os<<"------------------------------------------------------"<<std::endl;
    printer.print(Anasazi::Errors,os.str());

  }

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

  return 0;
}

