// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
/// \example LOBPCGEpetraExGen.cpp
/// \brief Use LOBPCG with Epetra, for a generalized eigenvalue problem.
///
/// This example computes the eigenvalues of largest magnitude of an
/// generalized eigenvalue problem, using Anasazi's implementation of
/// the LOBPCG method, with Epetra linear algebra.
///
/// The test problem claims to come from research described in the
/// following report: " A comparison of algorithms for modal analysis
/// in the absence of a sparse direct method", P. Arbenz, R. Lehoucq,
/// and U. Hetmaniuk, Sandia National Laboratories, Technical report
/// SAND2003-1028J.

#include "AnasaziConfigDefs.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziLOBPCGSolMgr.hpp"
#include "AnasaziBasicOutputManager.hpp"
#include "AnasaziEpetraAdapter.hpp"
#include "Epetra_CrsMatrix.h"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"

#include "ModeLaplace2DQ2.h"

using namespace Anasazi;

int main(int argc, char *argv[]) {

#ifdef HAVE_MPI
  // Initialize MPI
  //
  MPI_Init(&argc,&argv);
#endif

  bool success = false;
  try {
    // Create an Epetra communicator
    //
#ifdef HAVE_MPI
    Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm Comm;
#endif

    // Create an Anasazi output manager
    //
    BasicOutputManager<double> printer;
    printer.stream(Errors) << Anasazi_Version() << std::endl << std::endl;

    // Get the sorting std::string from the command line
    //
    std::string which("SM");
    Teuchos::CommandLineProcessor cmdp(false,true);
    cmdp.setOption("sort",&which,"Targetted eigenvalues (SM or LM).");
    if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
      throw -1;
    }


    typedef Epetra_MultiVector MV;
    typedef Epetra_Operator OP;
    typedef MultiVecTraits<double, Epetra_MultiVector> MVT;

    // Number of dimension of the domain
    const int space_dim = 2;

    // Size of each of the dimensions of the domain
    std::vector<double> brick_dim( space_dim );
    brick_dim[0] = 1.0;
    brick_dim[1] = 1.0;

    // Number of elements in each of the dimensions of the domain
    std::vector<int> elements( space_dim );
    elements[0] = 10;
    elements[1] = 10;

    // Create problem
    Teuchos::RCP<ModalProblem> testCase =
      Teuchos::rcp( new ModeLaplace2DQ2(Comm, brick_dim[0], elements[0], brick_dim[1], elements[1]) );

    // Get the stiffness and mass matrices
    Teuchos::RCP<Epetra_CrsMatrix> K = Teuchos::rcp( const_cast<Epetra_CrsMatrix *>(testCase->getStiffness()), false );
    Teuchos::RCP<Epetra_CrsMatrix> M = Teuchos::rcp( const_cast<Epetra_CrsMatrix *>(testCase->getMass()), false );

    // Eigensolver parameters
    int nev = 10;
    int blockSize = 5;
    int maxIters = 500;
    double tol = 1.0e-8;
    int verbosity = Anasazi::Errors + Anasazi::Warnings + Anasazi::FinalSummary;

    Teuchos::RCP<Epetra_MultiVector> ivec = Teuchos::rcp( new Epetra_MultiVector(K->OperatorDomainMap(), blockSize) );
    ivec->Random();

    // Create the eigenproblem.
    Teuchos::RCP<BasicEigenproblem<double, MV, OP> > MyProblem =
      Teuchos::rcp( new BasicEigenproblem<double, MV, OP>(K, M, ivec) );

    // Inform the eigenproblem that the operator A is symmetric
    MyProblem->setHermitian(true);

    // Set the number of eigenvalues requested
    MyProblem->setNEV( nev );

    // Inform the eigenproblem that you are finishing passing it information
    bool boolret = MyProblem->setProblem();
    if (boolret != true) {
      printer.print(Errors,"Anasazi::BasicEigenproblem::setProblem() returned an error.\n");
      throw -1;
    }

    // Create parameter list to pass into the solver manager
    //
    Teuchos::ParameterList MyPL;
    MyPL.set( "Which", which );
    MyPL.set( "Block Size", blockSize );
    MyPL.set( "Maximum Iterations", maxIters );
    MyPL.set( "Convergence Tolerance", tol );
    MyPL.set( "Full Ortho", true );
    MyPL.set( "Use Locking", true );
    MyPL.set( "Verbosity", verbosity );
    //
    // Create the solver manager
    LOBPCGSolMgr<double, MV, OP> MySolverMan(MyProblem, MyPL);

    // Solve the problem
    //
    ReturnType returnCode = MySolverMan.solve();

    // Get the eigenvalues and eigenvectors from the eigenproblem
    //
    Eigensolution<double,MV> sol = MyProblem->getSolution();
    std::vector<Value<double> > evals = sol.Evals;
    Teuchos::RCP<MV> evecs = sol.Evecs;

    // Compute residuals.
    //
    std::vector<double> normR(sol.numVecs);
    if (sol.numVecs > 0) {
      Teuchos::SerialDenseMatrix<int,double> T(sol.numVecs, sol.numVecs);
      Epetra_MultiVector Kvec( K->OperatorDomainMap(), evecs->NumVectors() );
      Epetra_MultiVector Mvec( M->OperatorDomainMap(), evecs->NumVectors() );
      T.putScalar(0.0);
      for (int i=0; i<sol.numVecs; i++) {
        T(i,i) = evals[i].realpart;
      }
      K->Apply( *evecs, Kvec );
      M->Apply( *evecs, Mvec );
      MVT::MvTimesMatAddMv( -1.0, Mvec, T, 1.0, Kvec );
      MVT::MvNorm( Kvec, normR );
    }

    // Print the results
    //
    std::ostringstream os;
    os.setf(std::ios_base::right, std::ios_base::adjustfield);
    os<<"Solver manager returned " << (returnCode == Converged ? "converged." : "unconverged.") << std::endl;
    os<<std::endl;
    os<<"------------------------------------------------------"<<std::endl;
    os<<std::setw(16)<<"Eigenvalue"
      <<std::setw(18)<<"Direct Residual"
      <<std::endl;
    os<<"------------------------------------------------------"<<std::endl;
    for (int i=0; i<sol.numVecs; i++) {
      os<<std::setw(16)<<evals[i].realpart
        <<std::setw(18)<<normR[i]/evals[i].realpart
        <<std::endl;
    }
    os<<"------------------------------------------------------"<<std::endl;
    printer.print(Errors,os.str());

    success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}
