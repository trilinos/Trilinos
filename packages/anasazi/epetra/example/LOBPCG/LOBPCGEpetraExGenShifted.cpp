// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
/// \example LOBPCGEpetraExGenShifted.cpp
/// \brief Use LOBPCG with Epetra, with shifted eigenvalue problem
///
/// This example computes the eigenvalues of largest magnitude of the
/// discretized 2-D Laplacian operator, using Anasazi's implementation
/// of the LOBPCG method.  This problem constructs a shifted
/// eigenproblem that targets the smallest eigenvalues around a
/// certain value (sigma).  This operator is discretized using linear
/// finite elements and constructed as an Epetra matrix, then passed
/// shifted using EpetraExt utilities.

// Include autoconfigured header
#include "AnasaziConfigDefs.hpp"

// Include header for LOBPCG solver
#include "AnasaziLOBPCGSolMgr.hpp"

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

#include "Teuchos_StandardCatchMacros.hpp"


int
main (int argc, char *argv[])
{
#ifdef EPETRA_MPI
  // Initialize MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  bool success = false;
  try {
    int i;
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
      throw -1;
    }

    //
    // ************************************
    // Start the LOBPCG iteration
    // ************************************
    //
    //  Variables used for the LOBPCG Method
    //
    const int    nev         = 10;
    const int    blockSize   = 5;
    const int    maxIters    = 500;
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
    MyPL.set( "Maximum Iterations", maxIters );
    MyPL.set( "Convergence Tolerance", tol );
    MyPL.set( "Full Ortho", true );
    MyPL.set( "Use Locking", true );
    MyPL.set( "Locking Tolerance", tol/10 );

    typedef Epetra_MultiVector MV;
    typedef Epetra_Operator OP;
    typedef Anasazi::MultiVecTraits<double, MV> MVT;
    // typedef Anasazi::OperatorTraits<double, MV, OP> OPT;

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
    if (!boolret) {
      printer.print(Anasazi::Errors,"Anasazi::BasicEigenproblem::setProblem() returned an error.\n");
      throw -1;
    }

    // Initialize the LOBPCG solver
    Anasazi::LOBPCGSolMgr<double, MV, OP> MySolverMgr(MyProblem, MyPL);

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

    success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}
