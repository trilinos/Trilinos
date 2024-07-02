// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
// This test compares the Anasazi solvers against ARPACK. The eigenproblems
// used are from the ARPACK examples: SYM, NONSYM, and COMPLEX
// See ARPACK_Operators.hpp and examlpesdesc for more information.

#include "AnasaziConfigDefs.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziBasicSort.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "AnasaziMVOPTester.hpp"
#include "AnasaziBasicOutputManager.hpp"
#include "AnasaziFactory.hpp"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

// templated multivector 
#include "MyMultiVec.hpp"
// ARPACK test problems
#include "ARPACK_Operators.hpp"

using namespace Teuchos;

int main(int argc, char *argv[]) 
{
#ifdef HAVE_MPI
  // Initialize MPI
  MPI_Init(&argc,&argv);
#endif

  typedef double ST;
  typedef ScalarTraits<ST>                   SCT;
  typedef SCT::magnitudeType                  MT;
  typedef Anasazi::MultiVec<ST>               MV;
  typedef Anasazi::Operator<ST>               OP;
  typedef Anasazi::MultiVecTraits<ST,MV>     MVT;
  typedef Anasazi::OperatorTraits<ST,MV,OP>  OPT;


  bool testFailed;
  bool verbose = false;
  std::string which("auto");
  int nx = 10;
  std::string problem("SDRV1");
  int ncv = 10;
  int nev = 5;
  MT tol = 1.0e-10;

  CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("sort",&which,"Targetted eigenvalues (auto, SR, LR, SI, LI, SM, LM).");
  cmdp.setOption("problem",&problem,"Problem to solve.");
  cmdp.setOption("nx",&nx,"Number of interior elements.");
  cmdp.setOption("ncv",&ncv,"Number of Arnoldi basis vectors.");
  cmdp.setOption("nev",&nev,"Number of Ritz values requested.");
  cmdp.setOption("tol",&tol,"Convergence tolerance.");
  if (cmdp.parse(argc,argv) != CommandLineProcessor::PARSE_SUCCESSFUL) {
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return -1;
  }


  // Create default output manager 
  RCP<Anasazi::OutputManager<ST> > MyOM = rcp( new Anasazi::BasicOutputManager<ST>() );
  // Set verbosity level
  int verbosity = Anasazi::Errors;
  if (verbose) {
    verbosity = Anasazi::Warnings + Anasazi::FinalSummary + Anasazi::TimingDetails;
  }
  MyOM->setVerbosity(verbosity);

  // print greeting
  MyOM->stream(Anasazi::Warnings) << Anasazi::Anasazi_Version() << std::endl << std::endl;

  // Eigensolver parameters
  int dim = nx*nx;
  int maxRestarts = 500;

  // Create initial vectors
  RCP<MV> ivec = rcp( new MyMultiVec<ST>(dim,1) );
  ivec->MvRandom();

  // Create matrices
  RCP< ARPACK_Example<ST> > prob;
  RCP<const OP> A, M, Op, B;

  prob = GetARPACKExample<ST>(problem,dim);
  if (!prob.get()) {
    MyOM->stream(Anasazi::Warnings)
      << "Invalid driver name. Try something like ""ndrv3"" or ""sdrv2""." << std::endl
      << "End Result: TEST FAILED" << std::endl;	
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return -1;
  }
  A = prob->getA();
  B = prob->getB();
  M = prob->getM();
  Op = prob->getOp();

  // determine sort
  which = prob->getSort();

  // test multivector and operators
  bool ierr;
  ierr = Anasazi::TestMultiVecTraits<ST,MV>(MyOM,ivec);
  MyOM->print(Anasazi::Warnings,"Testing MultiVector... ");
  if (ierr == true) {
    MyOM->print(Anasazi::Warnings,"PASSED TestMultiVecTraits()\n");
  } else {
    MyOM->print(Anasazi::Warnings,"FAILED TestMultiVecTraits()\n");
  }
  ierr = Anasazi::TestOperatorTraits<ST,MV,OP>(MyOM,ivec,A);
  MyOM->print(Anasazi::Warnings,"Testing OP... ");
  if (ierr == true) {
    MyOM->print(Anasazi::Warnings,"PASSED TestOperatorTraits()\n");
  } else {
    MyOM->print(Anasazi::Warnings,"FAILED TestOperatorTraits()\n");
  }
  ierr = Anasazi::TestOperatorTraits<ST,MV,OP>(MyOM,ivec,M);
  MyOM->print(Anasazi::Warnings,"Testing M... ");
  if (ierr == true) {
    MyOM->print(Anasazi::Warnings,"PASSED TestOperatorTraits()\n");
  } else {
    MyOM->print(Anasazi::Warnings,"FAILED TestOperatorTraits()\n");
  }

  // Create eigenproblem
  RCP<Anasazi::Eigenproblem<ST,MV,OP> > MyProblem = rcp( new Anasazi::BasicEigenproblem<ST,MV,OP>(Op, B, ivec) );
  //
  // Inform the eigenproblem if the operator is symmetric
  MyProblem->setHermitian(prob->isHerm());
  //
  // set the number of eigenvalues requested
  MyProblem->setNEV( nev );
  //
  // Inform the eigenproblem that you are done passing it information
  if (MyProblem->setProblem() != true) {
    MyOM->stream(Anasazi::Warnings)
      << "Anasazi::BasicEigenproblem::setProblem() had an error." << std::endl
      << "End Result: TEST FAILED" << std::endl;	
#ifdef HAVE_MPI
    MPI_Finalize() ;
#endif
    return -1;
  }

  // Create the solver manager
  Teuchos::ParameterList MyPL;
  MyPL.set( "Which", which );
  MyPL.set( "Block Size", 1 );
  MyPL.set( "Num Blocks", ncv );
  MyPL.set( "Maximum Restarts", maxRestarts );
  MyPL.set( "Orthogonalization", "DGKS" );
  MyPL.set( "Verbosity", verbosity );
  MyPL.set( "Convergence Tolerance", tol );
  RCP<Anasazi::SolverManager<ST,MV,OP>> BKS =
          Anasazi::Factory::create("BLOCK_KRYLOV_SCHUR", MyProblem, MyPL);

  // Solve the problem to the specified tolerances or length
  Anasazi::ReturnType ret = BKS->solve();
  testFailed = false;
  if (ret != Anasazi::Converged) {
    testFailed = true; 
  }

  // Get the eigenvalues and eigenvectors from the eigenproblem
  Anasazi::Eigensolution<ST,MV> sol = MyProblem->getSolution();
  RCP<MV> evecs = sol.Evecs;
  std::vector<Anasazi::Value<ST> > evals = sol.Evals;
  std::vector<int> index = sol.index;

  // Perform spectral transform on eigenvalues, if we used the 
  // spectral xformed operator (Op)
  // prob->xformeval(evals);

  /*
  // Compute the direct residual
  std::vector<MT> normV( nevecs );
  SerialDenseMatrix<int,ST> L(nevecs,nevecs);
  for (int i=0; i<nevecs; i++) {
    L(i,i) = (*evals)[i];
  }
  RCP<MV > Avecs = MVT::Clone( *evecs, nevecs );
  RCP<MV > Mvecs = MVT::Clone( *evecs, nevecs );

  OPT::Apply( *A, *evecs, *Avecs );
  OPT::Apply( *M, *evecs, *Mvecs );
  // Compute A*evecs - M*evecs*L
  MVT::MvTimesMatAddMv( -ONE, *Mvecs, L, ONE, *Avecs );
  MVT::MvNorm( *Avecs, normV );

  // check residuals
  for (int i=0; i<nevecs; i++) {
    if ( (*evals)[i] != SCT::zero() ) {
      normV[i] = SCT::magnitude(normV[i]/(*evals)[i]);
    }
    if ( normV[i] > ((MT)2.0)*tol ) {
      testFailed = true;
    }
  }

  {
    stringstream os;
    //      28,5,22
    os << "Back transformed eigenvalues     Relative Residual Norm" << std::endl
       << "-------------------------------------------------------" << std::endl;
    for (int i=0; i<nev; i++) {
      os.setf(ios::scientific, ios::floatfield);  
      os.precision(10);
      os << std::setw(28) << std::right << (*evals)[i] 
         << "     "
         << std::setw(22) << std::right << normV[i] 
         << std::endl;
    }
    MyOM->print(Anasazi::Warnings,os.str());
  }
  */

  // Exit
#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

  if (testFailed) {
    MyOM->print(Anasazi::Warnings,"End Result: TEST FAILED\n");
    return -1;
  }
  //
  // Default return value
  //
  MyOM->print(Anasazi::Warnings,"End Result: TEST PASSED\n");
  return 0;
}	
