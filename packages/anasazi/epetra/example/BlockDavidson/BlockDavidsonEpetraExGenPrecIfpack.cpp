// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#include "AnasaziConfigDefs.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziBlockDavidsonSolMgr.hpp"
#include "AnasaziBasicOutputManager.hpp"
#include "AnasaziEpetraAdapter.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_InvOperator.h"
#include "Teuchos_CommandLineProcessor.hpp"

// Include header for Ifpack incomplete Cholesky preconditioner
#include "Ifpack.h"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"

#include "ModeLaplace2DQ2.h"


using namespace Anasazi;


//*****************************************************************************
// Begin main routine
//*****************************************************************************
int main(int argc, char *argv[]) {

#ifdef HAVE_MPI
  // Initialize MPI
  MPI_Init(&argc,&argv);
#endif

  // Create an Epetra communicator
#ifdef HAVE_MPI
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  //************************************
  // Get the parameters from the command line
  //************************************
  //
  int    nev       = 10;
  int    blockSize = 10;
  int    numBlocks   = 4;
  int    maxRestarts = 100;
  double tol       = 1.0e-8;
  int numElements = 10;
  bool verbose = true;
  std::string which("SM");
  bool usePrec = true;
  double prec_dropTol = 1e-4;
  int prec_lofill = 0;
  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("nev",&nev,"Number of eigenpairs to compted.");
  cmdp.setOption("blockSize",&blockSize,"Block size.");
  cmdp.setOption("numBlocks",&numBlocks,"Number of blocks in basis.");
  cmdp.setOption("maxRestarts",&maxRestarts,"Maximum number of restarts.");
  cmdp.setOption("tol",&tol,"Relative convergence tolerance.");
  cmdp.setOption("numElements",&numElements,"Number of elements in the discretization.");
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("sort",&which,"Targetted eigenvalues (SM or LM).");
  cmdp.setOption("usePrec","noPrec",&usePrec,"Use Ifpack for preconditioning.");
  cmdp.setOption("prec_dropTol",&prec_dropTol,"Preconditioner: drop tolerance.");
  cmdp.setOption("prec_lofill",&prec_lofill,"Preconditioner: level of fill.");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return -1;
  }

  //************************************
  // Create an Anasazi output manager
  //************************************
  //
  // Set verbosity level
  int verbosity = Anasazi::Errors + Anasazi::Warnings;
  if (verbose) {
    verbosity += Anasazi::FinalSummary;
  }
  BasicOutputManager<double> printer(verbosity);
  printer.stream(Errors) << Anasazi_Version() << std::endl << std::endl;

  //************************************
  // Some useful typedefs
  //************************************
  //
  typedef Epetra_MultiVector MV;
  typedef Epetra_Operator OP;
  typedef MultiVecTraits<double, Epetra_MultiVector> MVT;

  //************************************
  // Create the problem matrices
  //************************************
  //
  printer.stream(Errors) << "Generating problem matrices..." << std::flush;
  // Number of dimension of the domain
  const int space_dim = 2;
  // Size of each of the dimensions of the domain
  std::vector<double> brick_dim( space_dim );
  brick_dim[0] = 1.0;
  brick_dim[1] = 1.0;
  // Number of elements in each of the dimensions of the domain
  std::vector<int> elements( space_dim );
  elements[0] = numElements;
  elements[1] = numElements;
  // Create problem
  Teuchos::RCP<ModalProblem> testCase =
    Teuchos::rcp( new ModeLaplace2DQ2(Comm, brick_dim[0], elements[0], brick_dim[1], elements[1]) );
  // Get the stiffness and mass matrices
  Teuchos::RCP<Epetra_CrsMatrix> K = Teuchos::rcp( const_cast<Epetra_CrsMatrix *>(testCase->getStiffness()), false );
  Teuchos::RCP<Epetra_CrsMatrix> M = Teuchos::rcp( const_cast<Epetra_CrsMatrix *>(testCase->getMass()), false );
  // tell the user that we're done
  printer.stream(Errors) << " done." << std::endl;


  //************************************
  // Select the Preconditioner
  //************************************
  //
  Teuchos::RCP<Ifpack_Preconditioner> prec;
  Teuchos::RCP<Epetra_Operator> PrecOp;
  if (usePrec) {
    printer.stream(Errors) << "Constructing Incomplete Cholesky preconditioner..." << std::flush;
    Ifpack precFactory;
    // additive-Schwartz incomplete Cholesky with thresholding; see IFPACK documentation
    std::string precType = "IC stand-alone";
    int overlapLevel = 0;
    prec = Teuchos::rcp( precFactory.Create(precType,K.get(),overlapLevel) );
    // parameters for preconditioner
    Teuchos::ParameterList precParams;
    precParams.set("fact: drop tolerance",prec_dropTol);
    precParams.set("fact: level-of-fill",prec_lofill);
    IFPACK_CHK_ERR(prec->SetParameters(precParams));
    IFPACK_CHK_ERR(prec->Initialize());
    IFPACK_CHK_ERR(prec->Compute());
    //
    printer.stream(Errors)
      << " done." << std::endl;
    // encapsulate this preconditioner into a IFPACKPrecOp class
    PrecOp = Teuchos::rcp( new Epetra_InvOperator(&*prec) );
  }


  //************************************
  // Call the BlockDavidson solver manager
  //***********************************
  // Create an Epetra_MultiVector for an initial vector to start the solver.
  // Note:  This needs to have the same number of columns as the blocksize.
  //
  Teuchos::RCP<Epetra_MultiVector> ivec
    = Teuchos::rcp( new Epetra_MultiVector(K->OperatorDomainMap(), blockSize) );
  ivec->Random();
  // Create the eigenproblem.
  //
  Teuchos::RCP<BasicEigenproblem<double, MV, OP> > MyProblem =
    Teuchos::rcp( new BasicEigenproblem<double, MV, OP>(K, M, ivec) );
  // Inform the eigenproblem that the operator K is symmetric
  //
  MyProblem->setHermitian(true);
  // Pass the preconditioner to the eigenproblem
  //
  if (usePrec) {
    MyProblem->setPrec(PrecOp);
  }
  // Set the number of eigenvalues requested
  //
  MyProblem->setNEV( nev );
  // Inform the eigenproblem that you are finishing passing it information
  //
  bool boolret = MyProblem->setProblem();
  if (boolret != true) {
    printer.print(Errors,"Anasazi::BasicEigenproblem::setProblem() returned an error.\n");
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return -1;
  }


  //************************************
  // Create parameter list to pass into the solver manager
  //************************************
  //
  Teuchos::ParameterList MyPL;
  MyPL.set( "Verbosity", verbosity );
  MyPL.set( "Which", which );
  MyPL.set( "Block Size", blockSize );
  MyPL.set( "Num Blocks", numBlocks );
  MyPL.set( "Maximum Restarts", maxRestarts );
  MyPL.set( "Convergence Tolerance", tol );
  //
  // Create the solver manager
  BlockDavidsonSolMgr<double, MV, OP> MySolverMan(MyProblem, MyPL);


  //************************************
  // Solve the problem
  //************************************
  //
  printer.stream(Errors) << "Solving eigenvalue problem..." << std::endl;
  ReturnType returnCode = MySolverMan.solve();
  // print some precond info
  if (usePrec) {
    printer.stream(FinalSummary) << *prec << std::endl;
  }


  //************************************
  // Get the eigenvalues and eigenvectors from the eigenproblem
  //************************************
  //
  Eigensolution<double,MV> sol = MyProblem->getSolution();
  std::vector<Value<double> > evals = sol.Evals;
  Teuchos::RCP<MV> evecs = sol.Evecs;


  //************************************
  // Compute residuals, just for funsies
  //************************************
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


  //************************************
  // Print the results
  //************************************
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

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return 0;
}

