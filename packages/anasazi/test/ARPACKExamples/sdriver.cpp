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
#include "AnasaziLOBPCG.hpp"
#include "AnasaziBlockDavidson.hpp"
#include "AnasaziBlockKrylovSchur.hpp"
#include "AnasaziBasicSort.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "AnasaziMVOPTester.hpp"
#include "AnasaziBasicOutputManager.hpp"

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
  int info = 0;

#ifdef HAVE_MPI 
  // Initialize MPI
  MPI_Init(&argc,&argv);
#endif

  bool testFailed;
  bool verbose = 0;
  std::string which("auto");
  int nx = 10;
  std::string problem("SDRV1");
  bool isherm;
  std::string solver("auto");

  CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("sort",&which,"Targetted eigenvalues (auto, SR, LR, SI, LI, SM, LM).");
  cmdp.setOption("nx",&nx,"Number of interior elements.");
  cmdp.setOption("problem",&problem,"Problem to solve.");
  cmdp.setOption("solver",&solver,"Eigensolver to use (LOBPCG, BKS, BD, auto)");
  if (cmdp.parse(argc,argv) != CommandLineProcessor::PARSE_SUCCESSFUL) {
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return -1;
  }

  typedef float ST;
  typedef ScalarTraits<ST>                   SCT;
  typedef SCT::magnitudeType                  MT;
  typedef Anasazi::MultiVec<ST>               MV;
  typedef Anasazi::Operator<ST>               OP;
  typedef Anasazi::MultiVecTraits<ST,MV>     MVT;
  typedef Anasazi::OperatorTraits<ST,MV,OP>  OPT;
  ST ONE  = SCT::one();


  // Create default output manager 
  RCP<Anasazi::OutputManager<ST> > MyOM 
    = rcp( new Anasazi::BasicOutputManager<ST>() );
  // Set verbosity level
  if (verbose) {
    MyOM->setVerbosity( Anasazi::Warning + Anasazi::FinalSummary );
  }

  // print greeting
  MyOM->stream(Anasazi::Warning) << Anasazi::Anasazi_Version() << endl << endl;

  Anasazi::ReturnType returnCode = Anasazi::Ok;  

  // Eigensolver parameters
  int dim = nx*nx;
  int nev = 4;
  int blockSize = 4;
  int maxIters = 500;
  int maxBlocks = 5;
  MT tol = 1.0e-4;

  // Create parameter list to pass into solver
  ParameterList MyPL;
  MyPL.set( "Block Size", blockSize );
  MyPL.set( "Max Iters", maxIters );
  MyPL.set( "Tol", tol );
  MyPL.set( "Max Blocks", maxBlocks );
  MyPL.set( "Max Restarts", maxIters );

  // Create initial vectors
  RCP<MV> ivec = rcp( new MyMultiVec<ST>(dim,blockSize) );
  ivec->MvRandom();

  // Create matrices
  RCP< ARPACK_Example<ST> > prob;
  RCP<const OP> A, M, Op, B;

  prob = GetARPACKExample<ST>(problem,dim);
  if (!prob.get()) {
    MyOM->stream(Anasazi::Warning)
      << "Invalid driver name. Try something like ""ndrv3"" or ""sdrv2""." << endl
      << "End Result: TEST FAILED" << endl;	
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return -1;
  }
  A = prob->getA();
  B = prob->getB();
  M = prob->getM();
  Op = prob->getOp();
  isherm = prob->isHerm();

  // determine solver
  if (solver == "auto") {
    if (isherm) {
      // solver = "LOBPCG";
      solver = "BKS";
    }
    else {
      solver = "BKS";
    }
  }

  // determine sort
  if (which == "auto") {
    if (solver == "LOBPCG" || solver == "BD") {
      which = "SM";
    }
    else {
      which = prob->getSort();
    }
  }

  // test multivector and operators
  int ierr;
  ierr = Anasazi::TestMultiVecTraits<ST,MV>(MyOM,ivec);
  MyOM->print(Anasazi::Warning,"Testing MultiVector... ");
  if (ierr == Anasazi::Ok) {      
    MyOM->print(Anasazi::Warning,"PASSED TestMultiVecTraits()\n");
  } else {
    MyOM->print(Anasazi::Warning,"FAILED TestMultiVecTraits()\n");
  }
  ierr = Anasazi::TestOperatorTraits<ST,MV,OP>(MyOM,ivec,A);
  MyOM->print(Anasazi::Warning,"Testing OP... ");
  if (ierr == Anasazi::Ok) {
    MyOM->print(Anasazi::Warning,"PASSED TestOperatorTraits()\n");
  } else {
    MyOM->print(Anasazi::Warning,"FAILED TestOperatorTraits()\n");
  }
  ierr = Anasazi::TestOperatorTraits<ST,MV,OP>(MyOM,ivec,M);
  MyOM->print(Anasazi::Warning,"Testing M... ");
  if (ierr == Anasazi::Ok) {
    MyOM->print(Anasazi::Warning,"PASSED TestOperatorTraits()\n");
  } else {
    MyOM->print(Anasazi::Warning,"FAILED TestOperatorTraits()\n");
  }

  // Create the sort manager
  RCP<Anasazi::BasicSort<ST,MV,OP> > MySM = 
     rcp( new Anasazi::BasicSort<ST,MV,OP>(which) );

  // Create eigenproblem
  RCP<Anasazi::BasicEigenproblem<ST,MV,OP> > MyProblem;
  if (solver == "LOBPCG" || solver == "BD") {
    MyProblem = rcp( new Anasazi::BasicEigenproblem<ST,MV,OP>(A, M, ivec) );
  }
  else {
    MyProblem = rcp( new Anasazi::BasicEigenproblem<ST,MV,OP>(Op, B, ivec) );
  }
  // Inform the eigenproblem if the operator is symmetric
  MyProblem->SetSymmetric(isherm);

  // Set the number of eigenvalues requested and the blocksize the solver should use
  MyProblem->SetNEV( nev );

  // Inform the eigenproblem that you are done passing it information
  info = MyProblem->SetProblem();
  if (info) {
    MyOM->stream(Anasazi::Warning)
      << "Anasazi::BasicEigenproblem::SetProblem() returned with code : "<< info << endl
      << "End Result: TEST FAILED" << endl;	
#ifdef HAVE_MPI
    MPI_Finalize() ;
#endif
    return -1;
  }

  // Create the eigensolver 
  RCP< Anasazi::Eigensolver<ST,MV,OP> > MySolver;
  if (solver == "LOBPCG") {
    MySolver = rcp( new Anasazi::LOBPCG<ST,MV,OP>(MyProblem, MySM, MyOM, MyPL));
  }
  else if (solver == "BKS") {
    MyPL.set( "Block Size", 1 );
    MyPL.set( "Max Blocks", 20 );
    MySolver = rcp( new Anasazi::BlockKrylovSchur<ST,MV,OP>(MyProblem, MySM, MyOM, MyPL));
  }
  else if (solver == "BD") {
    MySolver = rcp( new Anasazi::BlockDavidson<ST,MV,OP>(MyProblem, MySM, MyOM, MyPL));
  }
  else {
    MyOM->stream(Anasazi::Warning)
      << "Invalid solver: " << solver << endl
      << "End Result: TEST FAILED" << endl;	
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return -1;
  }
  MyOM->stream(Anasazi::Warning) << "Using solver: " << solver << endl;

  // Solve the problem to the specified tolerances or length
  returnCode = MySolver->solve();
  testFailed = false;
  if (returnCode != Anasazi::Ok) {
    testFailed = true; 
  }

  // Get the eigenvalues and eigenvectors from the eigenproblem
  RCP<std::vector<ST> > evals = MyProblem->GetEvals();
  RCP<MV > evecs = MyProblem->GetEvecs();
  int nevecs = MVT::GetNumberVecs(*evecs);

  // Perform spectral transform on eigenvalues, if we used the 
  // spectral xformed operator (Op)
  if (solver == "BKS") {
    prob->xformeval(*evals);
  }

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
    os << "Back transformed eigenvalues     Relative Residual Norm" << endl
       << "-------------------------------------------------------" << endl;
    for (int i=0; i<nev; i++) {
      os.setf(ios::scientific, ios::floatfield);  
      os.precision(10);
      os << std::setw(28) << std::right << (*evals)[i] 
         << "     "
         << std::setw(22) << std::right << normV[i] 
         << endl;
    }
    MyOM->print(Anasazi::Warning,os.str());
  }

  // Exit
#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

  if (testFailed) {
    MyOM->print(Anasazi::Warning,"End Result: TEST FAILED\n");
    return -1;
  }
  //
  // Default return value
  //
  MyOM->print(Anasazi::Warning,"End Result: TEST PASSED\n");
  return 0;

}	
