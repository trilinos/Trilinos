// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
//  This test is for the LOBPCG solver
//
#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"

#include "AnasaziTpetraAdapter.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziBasicOutputManager.hpp"
#include "AnasaziLOBPCGSolMgr.hpp"
#include "AnasaziSVQBOrthoManager.hpp"
#include "AnasaziBasicSort.hpp"
#include "AnasaziStatusTestMaxIters.hpp"
#include "AnasaziSolverUtils.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>

#include <ModeLaplace1DQ1.hpp>

using namespace Teuchos;
using namespace Anasazi;

typedef double                                      ST;
typedef double                           MT;
typedef Tpetra::MultiVector<ST>                     MV;
typedef Tpetra::Operator<ST>                        OP;
typedef MV::global_ordinal_type                     GO;
typedef MV::local_ordinal_type                      LO;
typedef MV::node_type                             Node;

class get_out : public std::logic_error {
  public: get_out(const std::string &whatarg) : std::logic_error(whatarg) {}
};

void checks( RCP<LOBPCG<ST,MV,OP> > solver, int blocksize, bool fullortho,
             RCP<Eigenproblem<ST,MV,OP> > problem,
             RCP<MatOrthoManager<ST,MV,OP> > ortho,
             SolverUtils<ST,MV,OP> &msutils)
{
  typedef MultiVecTraits<ST,MV>     MVT;
  typedef OperatorTraits<ST,MV,OP>  OPT;
  //typedef ScalarTraits<ST>          SCT; // unused
  LOBPCGState<ST,MV> state = solver->getState();

  TEUCHOS_TEST_FOR_EXCEPTION(MVT::GetNumberVecs(*state.X)  != solver->getBlockSize(),get_out,"blockSize() does not match allocated size for X");
  TEUCHOS_TEST_FOR_EXCEPTION(MVT::GetNumberVecs(*state.KX) != solver->getBlockSize(),get_out,"blockSize() does not match allocated size for KX");
  TEUCHOS_TEST_FOR_EXCEPTION(MVT::GetNumberVecs(*state.R)  != solver->getBlockSize(),get_out,"blockSize() does not match allocated size for R");
  TEUCHOS_TEST_FOR_EXCEPTION(MVT::GetNumberVecs(*state.H)  != solver->getBlockSize(),get_out,"blockSize() does not match allocated size for H");
  TEUCHOS_TEST_FOR_EXCEPTION(MVT::GetNumberVecs(*state.KH) != solver->getBlockSize(),get_out,"blockSize() does not match allocated size for KH");
  TEUCHOS_TEST_FOR_EXCEPTION(MVT::GetNumberVecs(*state.P)  != solver->getBlockSize(),get_out,"blockSize() does not match allocated size for P");
  TEUCHOS_TEST_FOR_EXCEPTION(MVT::GetNumberVecs(*state.KP) != solver->getBlockSize(),get_out,"blockSize() does not match allocated size for KP");
  if (solver->getProblem().getM() != null) {
    TEUCHOS_TEST_FOR_EXCEPTION(MVT::GetNumberVecs(*state.MX) != solver->getBlockSize(),get_out,"blockSize() does not match allocated size for MX");
    TEUCHOS_TEST_FOR_EXCEPTION(MVT::GetNumberVecs(*state.MH) != solver->getBlockSize(),get_out,"blockSize() does not match allocated size for MH");
    TEUCHOS_TEST_FOR_EXCEPTION(MVT::GetNumberVecs(*state.MP) != solver->getBlockSize(),get_out,"blockSize() does not match allocated size for MP");
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(state.MX != null,get_out,"MX should null; problem has no M matrix");
    TEUCHOS_TEST_FOR_EXCEPTION(state.MH != null,get_out,"MH should null; problem has no M matrix");
    TEUCHOS_TEST_FOR_EXCEPTION(state.MP != null,get_out,"MP should null; problem has no M matrix");
  }
  TEUCHOS_TEST_FOR_EXCEPTION(solver->getBlockSize() != blocksize, get_out,"Solver block size does not match specified block size.");
  TEUCHOS_TEST_FOR_EXCEPTION(solver->getFullOrtho() != fullortho, get_out,"Solver full ortho does not match specified state.");
  TEUCHOS_TEST_FOR_EXCEPTION(&solver->getProblem() != problem.get(),get_out,"getProblem() did not return the submitted problem.");
  TEUCHOS_TEST_FOR_EXCEPTION(solver->getMaxSubspaceDim() != 3*blocksize,get_out,"LOBPCG::getMaxSubspaceDim() should always be 3*blocksize");

  if (solver->isInitialized())
  {
    TEUCHOS_TEST_FOR_EXCEPTION(solver->getResNorms().size() != (unsigned int)blocksize,get_out,"getResNorms.size() does not match block size.");
    TEUCHOS_TEST_FOR_EXCEPTION(solver->getRes2Norms().size() != (unsigned int)blocksize,get_out,"getRes2Norms.size() does not match block size.");
    TEUCHOS_TEST_FOR_EXCEPTION(solver->getRitzRes2Norms().size() != (unsigned int)blocksize,get_out,"getRitzRes2Norms.size() does not match block size.");
    // check residual norms for consitency
    std::vector<double> RN2 = solver->getRes2Norms();
    std::vector<double> RR2 = solver->getRitzRes2Norms();
    for (int i=0; i<blocksize; i++) {
      TEUCHOS_TEST_FOR_EXCEPTION(RN2[i] != RR2[i],get_out,"getRitzRes2Norms() values do not match getRes2Norms() values.");
    }
    // check residuals
    RCP<const MV> evecs = state.X;
    RCP<MV> Kevecs, Mevecs;
    Kevecs = MVT::Clone(*evecs,blocksize);
    OPT::Apply(*problem->getOperator(),*evecs,*Kevecs);
    if (problem->getM() == null) {
      Mevecs = MVT::CloneCopy(*evecs);
    }
    else {
      Mevecs = MVT::Clone(*evecs,blocksize);
      OPT::Apply(*problem->getM(),*evecs,*Mevecs);
    }
    std::vector<Value<ST> > theta = solver->getRitzValues();
    TEUCHOS_TEST_FOR_EXCEPTION(theta.size() != (unsigned int)solver->getCurSubspaceDim(),get_out,"getRitzValues() has incorrect size.");
    SerialDenseMatrix<int,ST> T(blocksize,blocksize);
    for (int i=0; i<blocksize; i++) T(i,i) = theta[i].realpart;
    // LOBPCG computes residuals like R = K*X - M*X*T
    MVT::MvTimesMatAddMv(-1.0,*Mevecs,T,1.0,*Kevecs);
    MT error = msutils.errorEquality(*Kevecs,*state.R);
    // residuals from LOBPCG should be exact; we will cut a little slack
    TEUCHOS_TEST_FOR_EXCEPTION(error > 1e-14,get_out,"Residuals from solver did not match eigenvectors.");

    // check eigenvalues
    // X should be ritz vectors; they should diagonalize K to produce the current eigenvalues
    MT ninf = T.normInf();
    OPT::Apply(*problem->getOperator(),*evecs,*Kevecs);
    MVT::MvTransMv(1.0,*evecs,*Kevecs,T);
    for (int i=0; i<blocksize; i++) T(i,i) -= theta[i].realpart;
    error = T.normFrobenius() / ninf;
    TEUCHOS_TEST_FOR_EXCEPTION(error > 1e-14,get_out,"Ritz values don't match eigenvectors.");

    if (solver->hasP() && solver->getFullOrtho() == true) {
      // this should be small
      error = ortho->orthogError(*state.X,*state.P);
      TEUCHOS_TEST_FOR_EXCEPTION(error >= 1e-14,get_out,"FullOrtho enabled and P not orthogonal to X.");
    }
  }
  else {
    // not initialized
    TEUCHOS_TEST_FOR_EXCEPTION(solver->hasP() != false,get_out,"In unitialized state, hasP() should be false.");
    TEUCHOS_TEST_FOR_EXCEPTION(solver->getCurSubspaceDim() != 0,get_out,"In unitialized state, getCurSubspaceDim() should be 0.");
  }
}

void testsolver( RCP<BasicEigenproblem<ST,MV,OP> > problem,
                 RCP< OutputManager<ST> > printer,
                 RCP< MatOrthoManager<ST,MV,OP> > ortho,
                 RCP< SortManager<MT> > sorter,
                 ParameterList &pls,bool invalid=false)
{
  // create a status tester to run for one iteration
  RCP< StatusTest<ST,MV,OP> > tester = rcp( new StatusTestMaxIters<ST,MV,OP>(1) );

  // create the solver
  RCP< LOBPCG<ST,MV,OP> > solver;
  try {
    solver = rcp( new LOBPCG<ST,MV,OP>(problem,sorter,printer,tester,ortho,pls) );
    TEUCHOS_TEST_FOR_EXCEPTION(invalid, get_out, "Initializing with invalid parameters failed to throw exception.")
  }
  catch (const std::invalid_argument &ia) {
    TEUCHOS_TEST_FOR_EXCEPTION(!invalid, get_out, "Initializing with valid parameters unexpectadly threw exception.");
    // caught expected exception
    return;
  }


  const int  blocksize = pls.get<int>("Block Size");
  const bool fullortho = pls.get<bool>("Full Ortho");

  SolverUtils<ST,MV,OP> msutils;

  // solver should be uninitialized
  TEUCHOS_TEST_FOR_EXCEPTION(solver->isInitialized() != false,get_out,"Solver should be un-initialized after instantiation.");
  TEUCHOS_TEST_FOR_EXCEPTION(solver->getNumIters() != 0,get_out,"Number of iterations after initialization should be zero after init.")
  TEUCHOS_TEST_FOR_EXCEPTION(solver->hasP() != false,get_out,"Uninitialized solver should not have valid search directions.");
  TEUCHOS_TEST_FOR_EXCEPTION(solver->getAuxVecs().size() != 0,get_out,"getAuxVecs() should return empty.");
  checks(solver,blocksize,fullortho,problem,ortho,msutils);

  // initialize solver and perform checks
  solver->initialize();
  std::vector<Value<ST> > vals1 = solver->getRitzValues();
  vals1.resize(blocksize);
  MT sum1 = 0.0;
  for (int i=0; i<blocksize; i++) sum1 += vals1[i].realpart;
  TEUCHOS_TEST_FOR_EXCEPTION(solver->isInitialized() != true,get_out,"Solver should be initialized after call to initialize().");
  TEUCHOS_TEST_FOR_EXCEPTION(solver->getNumIters() != 0,get_out,"Number of iterations should be zero.")
  TEUCHOS_TEST_FOR_EXCEPTION(solver->hasP() != false,get_out,"Solver should not have valid P.");
  TEUCHOS_TEST_FOR_EXCEPTION(solver->getAuxVecs().size() != 0,get_out,"getAuxVecs() should return empty.");
  TEUCHOS_TEST_FOR_EXCEPTION(solver->getCurSubspaceDim() != blocksize,get_out,"after init, getCurSubspaceDim() should be blocksize.");
  checks(solver,blocksize,fullortho,problem,ortho,msutils);

  // call iterate(); solver should perform exactly one iteration and return; status test should be passed
  solver->iterate();
  std::vector<Value<ST> > vals2 = solver->getRitzValues();
  vals2.resize(blocksize);
  MT sum2 = 0.0;
  for (int i=0; i<blocksize; i++) sum2 += vals2[i].realpart;
  TEUCHOS_TEST_FOR_EXCEPTION(tester->getStatus() != Passed,get_out,"Solver returned from iterate() but getStatus() not Passed.");
  TEUCHOS_TEST_FOR_EXCEPTION(sum2 < sum1,get_out,"LOBPCG set to ascent method; sum of eigenvalues should increase.");
  TEUCHOS_TEST_FOR_EXCEPTION(solver->isInitialized() != true,get_out,"Solver should be initialized after call to initialize().");
  TEUCHOS_TEST_FOR_EXCEPTION(solver->getNumIters() != 1,get_out,"Number of iterations should be zero.")
  TEUCHOS_TEST_FOR_EXCEPTION(solver->hasP() != true,get_out,"Solver should have valid P.");
  TEUCHOS_TEST_FOR_EXCEPTION(solver->getAuxVecs().size() != 0,get_out,"getAuxVecs() should return empty.");
  TEUCHOS_TEST_FOR_EXCEPTION(solver->getCurSubspaceDim() != 2*blocksize,get_out,"after one step, getCurSubspaceDim() should be 2*blocksize.");
  checks(solver,blocksize,fullortho,problem,ortho,msutils);

  // reset numiters, call iterate(); solver should perform exactly one iteration and return; status test should be passed
  solver->resetNumIters();
  TEUCHOS_TEST_FOR_EXCEPTION(solver->getNumIters() != 0,get_out,"Number of iterations should be zero after resetNumIters().")
  solver->iterate();
  std::vector<Value<ST> > vals3 = solver->getRitzValues();
  vals3.resize(blocksize);
  MT sum3 = 0.0;
  for (int i=0; i<blocksize; i++) sum3 += vals3[i].realpart;
  TEUCHOS_TEST_FOR_EXCEPTION(tester->getStatus() != Passed,get_out,"Solver returned from iterate() but getStatus() not Passed.");
  TEUCHOS_TEST_FOR_EXCEPTION(sum3 < sum2,get_out,"LOBPCG set to ascent method; sum of eigenvalues should increase.");
  TEUCHOS_TEST_FOR_EXCEPTION(solver->isInitialized() != true,get_out,"Solver should be initialized after call to initialize().");
  TEUCHOS_TEST_FOR_EXCEPTION(solver->getNumIters() != 1,get_out,"Number of iterations should be zero.")
  TEUCHOS_TEST_FOR_EXCEPTION(solver->hasP() != true,get_out,"Solver should have valid P.");
  TEUCHOS_TEST_FOR_EXCEPTION(solver->getAuxVecs().size() != 0,get_out,"getAuxVecs() should return empty.");
  TEUCHOS_TEST_FOR_EXCEPTION(solver->getCurSubspaceDim() != 3*blocksize,get_out,"after two steps, getCurSubspaceDim() should be 3*blocksize.");
  checks(solver,blocksize,fullortho,problem,ortho,msutils);

  // call setBlockSize with current blocksize and ensure that it wasn't resized and that it is still initialized
  solver->setBlockSize(blocksize);
  TEUCHOS_TEST_FOR_EXCEPTION(!solver->isInitialized(),get_out,"After trivial setBlockSize(), solver should still be initialized.");
  TEUCHOS_TEST_FOR_EXCEPTION(solver->getBlockSize() != blocksize,get_out,"After trivial setBlockSize(), solver should have same block size.");

  // decrease block size and see the difference
  solver->setBlockSize(blocksize-1);
  TEUCHOS_TEST_FOR_EXCEPTION(!solver->isInitialized(),get_out,"After decreasing setBlockSize(), solver should be initialized.");
  TEUCHOS_TEST_FOR_EXCEPTION(solver->getBlockSize() != blocksize-1,get_out,"After setBlockSize(), new block size was not in effect.");
  TEUCHOS_TEST_FOR_EXCEPTION(solver->getMaxSubspaceDim() != 3*solver->getBlockSize(),get_out,"After setBlockSize(), getMaxSubspaceDim() should be changed.");

  // increase block size and see the difference
  solver->setBlockSize(blocksize);
  TEUCHOS_TEST_FOR_EXCEPTION(solver->isInitialized(),get_out,"After increasing setBlockSize(), solver should be uninitialized.");
  TEUCHOS_TEST_FOR_EXCEPTION(solver->getBlockSize() != blocksize,get_out,"After setBlockSize(), new block size was not in effect.");
  TEUCHOS_TEST_FOR_EXCEPTION(solver->getCurSubspaceDim() != 0,get_out,"Unitialized solver should have getCurSubspaceDim() of 0.");
  TEUCHOS_TEST_FOR_EXCEPTION(solver->getMaxSubspaceDim() != 3*solver->getBlockSize(),get_out,"After setBlockSize(), getMaxSubspaceDim() should be changed.");
}

int main(int argc, char *argv[])
{
  Tpetra::ScopeGuard tpetraScope (&argc,&argv);
  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm ();

  bool testFailed;
  bool verbose = false;

  CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("debug","quiet",&verbose,"Print messages and results.");
  if (cmdp.parse(argc,argv) != CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }

  // create the output manager
  RCP< OutputManager<ST> > printer = rcp( new BasicOutputManager<ST>() );

  if (verbose) printer->stream(Errors) << Anasazi_Version() << std::endl << std::endl;

  const int veclength = 99;

  //  Problem information
  int space_dim = 1;
  std::vector<double> brick_dim( space_dim );
  brick_dim[0] = 1.0;
  std::vector<int> elements( space_dim );
  elements[0] = veclength+1;

  // Create problem
  ModeLaplace1DQ1<ST,LO,GO,Node> testCase( comm, brick_dim[0], elements[0]);
  //
  // Get the stiffness and mass matrices
  RCP<const Tpetra::CrsMatrix<ST,LO,GO,Node> > K = testCase.getStiffness();
  RCP<const Tpetra::CrsMatrix<ST,LO,GO,Node> > M = testCase.getMass();
  //
  // Create the initial vectors
  const int nev = 4;
  RCP<MV> ivec = rcp (new MV (K->getDomainMap(), nev));
  ivec->randomize();
  //
  // Create eigenproblem: one standard and one generalized
  RCP<BasicEigenproblem<ST,MV,OP> > probstd = rcp( new BasicEigenproblem<ST, MV, OP>(K, ivec) );
  RCP<BasicEigenproblem<ST,MV,OP> > probgen = rcp( new BasicEigenproblem<ST, MV, OP>(K, M, ivec) );
  //
  // Inform the eigenproblem that the operator A is symmetric
  probstd->setHermitian(true);
  probgen->setHermitian(true);
  //
  // Set the number of eigenvalues requested
  probstd->setNEV( nev );
  probgen->setNEV( nev );
  //
  // Inform the eigenproblem that you are finishing passing it information
  if ( probstd->setProblem() != true || probgen->setProblem() != true ) {
    if (verbose) {
      printer->stream(Errors) << "Anasazi::BasicEigenproblem::SetProblem() returned with error." << std::endl
                              << "End Result: TEST FAILED" << std::endl;
    }
    return -1;
  }

  // create the orthogonalization managers: one standard and one M-based
  RCP< MatOrthoManager<ST,MV,OP> > orthostd = rcp( new SVQBOrthoManager<ST,MV,OP>() );
  RCP< MatOrthoManager<ST,MV,OP> > orthogen = rcp( new SVQBOrthoManager<ST,MV,OP>(M) );
  // create the sort manager
  RCP< SortManager<MT> > sorter = rcp( new BasicSort<MT>("LM") );
  // create the parameter list specifying blocksize > nev and full orthogonalization
  ParameterList pls;

  // begin testing
  testFailed = false;

  try
  {

    if (verbose) printer->stream(Errors) << "Testing solver(default,default) with standard eigenproblem..." << std::endl;
    testsolver(probstd,printer,orthostd,sorter,pls);
    if (verbose) printer->stream(Errors) << "Testing solver(default,default) with generalized eigenproblem..." << std::endl;
    testsolver(probgen,printer,orthogen,sorter,pls);

    pls.set<int>("Block Size",nev);
    pls.set<bool>("Full Ortho",false);
    if (verbose) printer->stream(Errors) << "Testing solver(nev,false) with standard eigenproblem..." << std::endl;
    testsolver(probstd,printer,orthostd,sorter,pls);
    TEUCHOS_TEST_FOR_EXCEPTION(pls.getEntryPtr("Block Size")->isUsed() == false, get_out, "Solver did not consume parameter \"Block Size\".");
    TEUCHOS_TEST_FOR_EXCEPTION(pls.getEntryPtr("Full Ortho")->isUsed() == false, get_out, "Solver did not consume parameter \"Full Ortho\".");
    pls.set<bool>("Full Ortho",true);
    if (verbose) printer->stream(Errors) << "Testing solver(nev,true) with standard eigenproblem..." << std::endl;
    testsolver(probstd,printer,orthostd,sorter,pls);
    pls.set<bool>("Full Ortho",false);
    if (verbose) printer->stream(Errors) << "Testing solver(nev,false) with generalized eigenproblem..." << std::endl;
    testsolver(probgen,printer,orthogen,sorter,pls);
    pls.set<bool>("Full Ortho",true);
    if (verbose) printer->stream(Errors) << "Testing solver(nev,true) with generalized eigenproblem..." << std::endl;
    testsolver(probgen,printer,orthogen,sorter,pls);

    pls.set<int>("Block Size",2*nev);
    pls.set<bool>("Full Ortho",false);
    if (verbose) printer->stream(Errors) << "Testing solver(2*nev,false) with standard eigenproblem..." << std::endl;
    testsolver(probstd,printer,orthostd,sorter,pls);
    pls.set<bool>("Full Ortho",true);
    if (verbose) printer->stream(Errors) << "Testing solver(2*nev,true) with standard eigenproblem..." << std::endl;
    testsolver(probstd,printer,orthostd,sorter,pls);
    pls.set<bool>("Full Ortho",false);
    if (verbose) printer->stream(Errors) << "Testing solver(2*nev,false) with generalized eigenproblem..." << std::endl;
    testsolver(probgen,printer,orthogen,sorter,pls);
    pls.set<bool>("Full Ortho",true);
    if (verbose) printer->stream(Errors) << "Testing solver(2*nev,true) with generalized eigenproblem..." << std::endl;
    testsolver(probgen,printer,orthogen,sorter,pls);

    pls.set<int>("Block Size",nev/2);
    pls.set<bool>("Full Ortho",false);
    if (verbose) printer->stream(Errors) << "Testing solver(nev/2,false) with standard eigenproblem..." << std::endl;
    testsolver(probstd,printer,orthostd,sorter,pls);
    pls.set<bool>("Full Ortho",true);
    if (verbose) printer->stream(Errors) << "Testing solver(nev/2,true) with standard eigenproblem..." << std::endl;
    testsolver(probstd,printer,orthostd,sorter,pls);
    pls.set<bool>("Full Ortho",false);
    if (verbose) printer->stream(Errors) << "Testing solver(nev/2,false) with generalized eigenproblem..." << std::endl;
    testsolver(probgen,printer,orthogen,sorter,pls);
    pls.set<bool>("Full Ortho",true);
    if (verbose) printer->stream(Errors) << "Testing solver(nev/2,true) with generalized eigenproblem..." << std::endl;
    testsolver(probgen,printer,orthogen,sorter,pls);

    // try with an invalid block size
    pls.set<int>("Block Size",0);
    if (verbose) printer->stream(Errors) << "Testing solver(0) with standard eigenproblem..." << std::endl;
    testsolver(probstd,printer,orthostd,sorter,pls,true);

    // try with a too-large block size
    pls.set<int>("Num Blocks",veclength+1);
    if (verbose) printer->stream(Errors) << "Testing solver(toomany) with standard eigenproblem..." << std::endl;
    testsolver(probstd,printer,orthostd,sorter,pls,true);

    // try with an unset problem
    // setHermitian will mark the problem as unset
    probstd->setHermitian(false);
    if (verbose) printer->stream(Errors) << "Testing solver with unset eigenproblem..." << std::endl;
    testsolver(probstd,printer,orthostd,sorter,pls,true);

    // set the problem, and try with a non-Hermitian problem
    if ( probstd->setProblem() != true ) {
      if (verbose) {
        printer->stream(Errors) << "Anasazi::BasicEigenproblem::SetProblem() returned with error." << std::endl
                                << "End Result: TEST FAILED" << std::endl;
      }
      return -1;
    }
    if (verbose) printer->stream(Errors) << "Testing solver with non-Hermitian eigenproblem..." << std::endl;
    testsolver(probstd,printer,orthostd,sorter,pls,true);
    // fix it now
    probstd->setHermitian(true);
    probstd->setProblem();

    // create a dummy status tester
    RCP< StatusTest<ST,MV,OP> > dumtester = rcp( new StatusTestMaxIters<ST,MV,OP>(1) );

    // try with a null problem
    if (verbose) printer->stream(Errors) << "Testing solver with null eigenproblem..." << std::endl;
    try {
      RCP< LOBPCG<ST,MV,OP> > solver
        = rcp( new LOBPCG<ST,MV,OP>(Teuchos::null,sorter,printer,dumtester,orthostd,pls) );
      TEUCHOS_TEST_FOR_EXCEPTION(true,get_out,"Initializing with invalid parameters failed to throw exception.");
    }
    catch (const std::invalid_argument &ia) {
      // caught expected exception
    }

    // try with a null sortman
    if (verbose) printer->stream(Errors) << "Testing solver with null sort manager..." << std::endl;
    try {
      RCP< LOBPCG<ST,MV,OP> > solver
        = rcp( new LOBPCG<ST,MV,OP>(probstd,Teuchos::null,printer,dumtester,orthostd,pls) );
      TEUCHOS_TEST_FOR_EXCEPTION(true,get_out,"Initializing with invalid parameters failed to throw exception.");
    }
    catch (const std::invalid_argument &ia) {
      // caught expected exception
    }

    // try with a output man problem
    if (verbose) printer->stream(Errors) << "Testing solver with null output manager..." << std::endl;
    try {
      RCP< LOBPCG<ST,MV,OP> > solver
        = rcp( new LOBPCG<ST,MV,OP>(probstd,sorter,Teuchos::null,dumtester,orthostd,pls) );
      TEUCHOS_TEST_FOR_EXCEPTION(true,get_out,"Initializing with invalid parameters failed to throw exception.");
    }
    catch (const std::invalid_argument &ia) {
      // caught expected exception
    }

    // try with a null status test
    if (verbose) printer->stream(Errors) << "Testing solver with null status test..." << std::endl;
    try {
      RCP< LOBPCG<ST,MV,OP> > solver
        = rcp( new LOBPCG<ST,MV,OP>(probstd,sorter,printer,Teuchos::null,orthostd,pls) );
      TEUCHOS_TEST_FOR_EXCEPTION(true,get_out,"Initializing with invalid parameters failed to throw exception.");
    }
    catch (const std::invalid_argument &ia) {
      // caught expected exception
    }

    // try with a null orthoman
    if (verbose) printer->stream(Errors) << "Testing solver with null ortho manager..." << std::endl;
    try {
      RCP< LOBPCG<ST,MV,OP> > solver
        = rcp( new LOBPCG<ST,MV,OP>(probstd,sorter,printer,dumtester,Teuchos::null,pls) );
      TEUCHOS_TEST_FOR_EXCEPTION(true,get_out,"Initializing with invalid parameters failed to throw exception.");
    }
    catch (const std::invalid_argument &ia) {
      // caught expected exception
    }
  }
  catch (const get_out &go) {
    printer->stream(Errors) << "Test failed: " << go.what() << std::endl;
    testFailed = true;
  }
  catch (const std::exception &e) {
    printer->stream(Errors) << "Caught unexpected exception: " << e.what() << std::endl;
    testFailed = true;
  }

  if (testFailed) {
    if (verbose) {
      printer->stream(Errors) << std::endl << "End Result: TEST FAILED" << std::endl;
    }
    return -1;
  }
  //
  // Default return value
  //
  if (verbose) {
    printer->stream(Errors) << std::endl << "End Result: TEST PASSED" << std::endl;
  }
  return 0;

}
