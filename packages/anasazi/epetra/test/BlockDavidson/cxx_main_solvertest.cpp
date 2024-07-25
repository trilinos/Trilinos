// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
//  This test is for the BlockDavidson solver
//
#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"

#include "AnasaziEpetraAdapter.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"

#include "AnasaziBlockDavidson.hpp"

#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziBasicOutputManager.hpp"
#include "AnasaziSVQBOrthoManager.hpp"
#include "AnasaziBasicSort.hpp"
#include "AnasaziStatusTestMaxIters.hpp"
#include "AnasaziSolverUtils.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif

#include "ModeLaplace1DQ1.h"

namespace { // (anonymous)

using namespace Anasazi;
using namespace Teuchos;
typedef double                                  ScalarType;
typedef ScalarTraits<ScalarType>::magnitudeType MagnitudeType;
typedef ScalarTraits<MagnitudeType>             MT;
typedef Epetra_MultiVector                 MV;
typedef Epetra_Operator                    OP;
typedef MultiVecTraits<ScalarType,MV>     MVT;
typedef OperatorTraits<ScalarType,MV,OP>  OPT;

class get_out : public std::logic_error {
  public: get_out(const std::string &whatarg) : std::logic_error(whatarg) {}
};

void checks( RCP<BlockDavidson<ScalarType,MV,OP> > solver, int blocksize, int numblocks,
             RCP<Eigenproblem<ScalarType,MV,OP> > problem,
             RCP<MatOrthoManager<ScalarType,MV,OP> > ortho,
             SolverUtils<ScalarType,MV,OP> &msutils) {
  BlockDavidsonState<ScalarType,MV> state = solver->getState();

  TEUCHOS_TEST_FOR_EXCEPTION(MVT::GetNumberVecs(*state.V)  != solver->getMaxSubspaceDim(),get_out,     "getMaxSubspaceDim() does not match allocated size for V.");
  TEUCHOS_TEST_FOR_EXCEPTION(MVT::GetNumberVecs(*state.X)  != solver->getBlockSize(),get_out,"blockSize() does not match allocated size for X.");
  TEUCHOS_TEST_FOR_EXCEPTION(MVT::GetNumberVecs(*state.KX) != solver->getBlockSize(),get_out,"blockSize() does not match allocated size for KX.");
  TEUCHOS_TEST_FOR_EXCEPTION(MVT::GetNumberVecs(*solver->getRitzVectors()) != solver->getBlockSize(),get_out,"blockSize() does not match getRitzVectors().");
  TEUCHOS_TEST_FOR_EXCEPTION(state.T->size() != (unsigned int)solver->getCurSubspaceDim(),get_out,"state.T->size() does not match getCurSubspaceDim().");
  if (solver->getProblem().getM() != null) {
    TEUCHOS_TEST_FOR_EXCEPTION(MVT::GetNumberVecs(*state.MX) != solver->getBlockSize(),get_out,"blockSize() does not match allocated size for MX.");
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(state.MX != null,get_out,"MX should null; problem has no M matrix.");
  }
  TEUCHOS_TEST_FOR_EXCEPTION(MVT::GetNumberVecs(*state.R)  != solver->getBlockSize(),get_out,"blockSize() does not match allocated size for R.");
  TEUCHOS_TEST_FOR_EXCEPTION(solver->getBlockSize() != blocksize, get_out,"Solver block size does not match specified block size.");
  TEUCHOS_TEST_FOR_EXCEPTION(solver->getMaxSubspaceDim()/solver->getBlockSize() != numblocks, get_out, "Solver num blaocks does not match specified num blocks.");
  TEUCHOS_TEST_FOR_EXCEPTION(&solver->getProblem() != problem.get(),get_out,"getProblem() did not return the submitted problem.");
  TEUCHOS_TEST_FOR_EXCEPTION(solver->getMaxSubspaceDim() != numblocks*blocksize,get_out,"BlockDavidson::getMaxSubspaceDim() does not match numblocks*blocksize.");

  if (solver->isInitialized())
  {
    TEUCHOS_TEST_FOR_EXCEPTION(solver->getResNorms().size() != (unsigned int)blocksize,get_out,"getResNorms.size() does not match block size.");
    TEUCHOS_TEST_FOR_EXCEPTION(solver->getRes2Norms().size() != (unsigned int)blocksize,get_out,"getRes2Norms.size() does not match block size.");
    TEUCHOS_TEST_FOR_EXCEPTION(solver->getRitzRes2Norms().size() != (unsigned int)blocksize,get_out,"getRitzRes2Norms.size() does not match size.");
    // check residual norms for consitency
    std::vector<double> RN2 = solver->getRes2Norms();
    std::vector<double> RR2 = solver->getRitzRes2Norms();
    for (int i=0; i<blocksize; i++) {
      TEUCHOS_TEST_FOR_EXCEPTION(RN2[i] != RR2[i],get_out,"getRitzRes2Norms() values do not match getRes2Norms() values.");
    }
    // check ritz values
    std::vector<Value<ScalarType> > theta = solver->getRitzValues();
    TEUCHOS_TEST_FOR_EXCEPTION(theta.size() != (unsigned int)solver->getCurSubspaceDim(),get_out,"getRitzValues().size() does not match getCurSubspaceDim().");
    for (unsigned int i=0; i<theta.size(); i++) {
      TEUCHOS_TEST_FOR_EXCEPTION(theta[i].imagpart != MT::zero(),get_out,"getRitzValues() returned complex eigenvalues.");
    }
    // check ritz index
    std::vector<int> index = solver->getRitzIndex();
    TEUCHOS_TEST_FOR_EXCEPTION( index.size() != (unsigned int)solver->getCurSubspaceDim(), get_out, "Ritz index size not consistent with eigenvector size.");
    for (unsigned int i=0; i<index.size(); i++) {
      TEUCHOS_TEST_FOR_EXCEPTION(index[i] != 0,get_out,"Ritz index contained non-zeros.");
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
    SerialDenseMatrix<int,ScalarType> T(blocksize,blocksize);
    for (int i=0; i<blocksize; i++) T(i,i) = theta[i].realpart;
    // BlockDavidson computes residuals like R = K*X - M*X*T
    MVT::MvTimesMatAddMv(-1.0,*Mevecs,T,1.0,*Kevecs);
    MagnitudeType error = msutils.errorEquality(*Kevecs,*state.R);
    // residuals from BlockDavidson should be exact; we will cut a little slack
    TEUCHOS_TEST_FOR_EXCEPTION(error > 1e-14,get_out,"Residuals from solver did not match eigenvectors.");

    // check eigenvalues
    // X should be ritz vectors; they should diagonalize K to produce the current eigenvalues
    // this is the largest magnitude ritz value
    MagnitudeType ninf = T.normInf();
    OPT::Apply(*problem->getOperator(),*evecs,*Kevecs);
    MVT::MvTransMv(1.0,*evecs,*Kevecs,T);
    for (int i=0; i<blocksize; i++) T(i,i) -= theta[i].realpart;
    error = T.normFrobenius() / ninf;
    TEUCHOS_TEST_FOR_EXCEPTION(error > 1e-12,get_out,"Disagreement between Ritz vectors and Ritz values.");
  }
  else {
    // not initialized
    TEUCHOS_TEST_FOR_EXCEPTION(solver->getCurSubspaceDim() != 0,get_out,"In unitialized state, getCurSubspaceDim() should be 0.");
  }
}

void testsolver( RCP<BasicEigenproblem<ScalarType,MV,OP> > problem,
                 RCP< OutputManager<ScalarType> > printer,
                 RCP< MatOrthoManager<ScalarType,MV,OP> > ortho,
                 RCP< SortManager<MagnitudeType> > sorter,
                 ParameterList &pls,bool invalid,
                 BlockDavidsonState<ScalarType,MV> initstate, bool invalidinit)
{
  // create a status tester to run for two iterations
  RCP< StatusTest<ScalarType,MV,OP> > tester = rcp( new StatusTestMaxIters<ScalarType,MV,OP>(2) );

  const int blocksize = pls.get<int>("Block Size",problem->getNEV());
  const int numblocks = pls.get<int>("Num Blocks",2);

  // create the solver
  RCP< BlockDavidson<ScalarType,MV,OP> > solver;
  try {
    solver = rcp( new BlockDavidson<ScalarType,MV,OP>(problem,sorter,printer,tester,ortho,pls) );
    TEUCHOS_TEST_FOR_EXCEPTION(invalid, get_out, "Instantiating with invalid parameters failed to throw exception.")
  }
  catch (const std::invalid_argument &ia) {
    TEUCHOS_TEST_FOR_EXCEPTION(!invalid, get_out, "Instantiating with valid parameters unexpectadly threw exception.");
    // caught expected exception
    return;
  }

  SolverUtils<ScalarType,MV,OP> msutils;

  // solver should be uninitialized
  TEUCHOS_TEST_FOR_EXCEPTION(solver->isInitialized() != false,get_out,"Solver should be un-initialized after instantiation.");
  TEUCHOS_TEST_FOR_EXCEPTION(solver->getNumIters() != 0,get_out,"Number of iterations after initialization should be zero after init.")
  TEUCHOS_TEST_FOR_EXCEPTION(solver->getAuxVecs().size() != 0,get_out,"getAuxVecs() should return empty.");
  checks(solver,blocksize,numblocks,problem,ortho,msutils);

  int cursize  = MVT::GetNumberVecs(*problem->getInitVec());
  if (cursize < blocksize) cursize = blocksize;

  // initialize solver and perform checks
  try {
    solver->initialize(initstate);
    TEUCHOS_TEST_FOR_EXCEPTION(invalidinit, get_out, "Initializing with invalid data failed to throw exception.")
  }
  catch (const std::invalid_argument &ia) {
    TEUCHOS_TEST_FOR_EXCEPTION(!invalidinit, get_out, "Initializing with valid data unexpectadly threw exception.");
    // caught expected exception
    return;
  }

  TEUCHOS_TEST_FOR_EXCEPTION(solver->isInitialized() != true,get_out,"Solver should be initialized after call to initialize().");
  TEUCHOS_TEST_FOR_EXCEPTION(solver->getNumIters() != 0,get_out,"Number of iterations should be zero.")
  TEUCHOS_TEST_FOR_EXCEPTION(solver->getAuxVecs().size() != 0,get_out,"getAuxVecs() should return empty.");
  TEUCHOS_TEST_FOR_EXCEPTION(solver->getCurSubspaceDim() != cursize,get_out,"after init, getCurSubspaceDim() should be size of problem->getInitVec().");
  checks(solver,blocksize,numblocks,problem,ortho,msutils);

  while ( solver->getCurSubspaceDim() < solver->getMaxSubspaceDim() ) {

    // call iterate(); solver should perform at most two iterations and return; status test should be consistent
    solver->iterate();
    TEUCHOS_TEST_FOR_EXCEPTION(solver->isInitialized() != true,get_out,"Solver should be initialized after return from iterate().");
    TEUCHOS_TEST_FOR_EXCEPTION(solver->getNumIters() != 2 && tester->getStatus() == Passed,get_out,"Number of iterations not consistent with StatusTest return.");
    TEUCHOS_TEST_FOR_EXCEPTION(solver->getNumIters() == 2 && tester->getStatus() != Passed,get_out,"Number of iterations not consistent with StatusTest return.");
    TEUCHOS_TEST_FOR_EXCEPTION(solver->getCurSubspaceDim() != solver->getMaxSubspaceDim() && tester->getStatus() != Passed,get_out,"solver should not have returned from iterate().");
    TEUCHOS_TEST_FOR_EXCEPTION(solver->getAuxVecs().size() != 0,get_out,"getAuxVecs() should return empty.");
    TEUCHOS_TEST_FOR_EXCEPTION(solver->getCurSubspaceDim() != cursize+solver->getNumIters()*blocksize,get_out,"getCurSubspaceDim() did not grow as expected.");
    TEUCHOS_TEST_FOR_EXCEPTION(solver->getCurSubspaceDim() > solver->getMaxSubspaceDim(),get_out,"impossibly large basis.");
    checks(solver,blocksize,numblocks,problem,ortho,msutils);
    // record current size
    cursize = solver->getCurSubspaceDim();

    // reset numiters and check
    solver->resetNumIters();
    TEUCHOS_TEST_FOR_EXCEPTION(solver->getNumIters() != 0,get_out,"Number of iterations should be zero after resetNumIters().")
  }

  // call setSize with current blocksize,numblocks and ensure that it wasn't resized and that it is still initialized
  solver->setSize(blocksize,numblocks);
  TEUCHOS_TEST_FOR_EXCEPTION(!solver->isInitialized(),get_out,"After trivial setSize(), solver should still be initialized.");
  TEUCHOS_TEST_FOR_EXCEPTION(solver->getBlockSize() != blocksize,get_out,"After trivial setSize(), solver should have same block size.");
  TEUCHOS_TEST_FOR_EXCEPTION(solver->getMaxSubspaceDim()/blocksize != numblocks,get_out,"After trivial setSize(), solver should have same num blocks.");

  // change block size and see the difference
  solver->setBlockSize(blocksize+1);
  TEUCHOS_TEST_FOR_EXCEPTION(solver->isInitialized(),get_out,"After setBlockSize(), solver should be uninitialized.");
  TEUCHOS_TEST_FOR_EXCEPTION(solver->getCurSubspaceDim() != 0,get_out,"After setBlocksize(): Uninitialized solver should have getCurSubspaceDim() == 0.");
  TEUCHOS_TEST_FOR_EXCEPTION(solver->getBlockSize() != blocksize+1,get_out,"After setBlockSize(), new block size was not in effect.");
  TEUCHOS_TEST_FOR_EXCEPTION(solver->getMaxSubspaceDim()/solver->getBlockSize() != numblocks,get_out,"After setBlockSize(), num blocks should not have changed.");
  // call setSize and see the difference
  solver->initialize();
  solver->setSize(blocksize,numblocks+1);
  TEUCHOS_TEST_FOR_EXCEPTION(solver->isInitialized(),get_out,"After setSize(), solver should be uninitialized.");
  TEUCHOS_TEST_FOR_EXCEPTION(solver->getCurSubspaceDim() != 0,get_out,"After setSize(): Uninitialized solver should have getCurSubspaceDim() == 0.");
  TEUCHOS_TEST_FOR_EXCEPTION(solver->getBlockSize() != blocksize,get_out,"After setSize(), new block size was not in effect.");
  TEUCHOS_TEST_FOR_EXCEPTION(solver->getMaxSubspaceDim()/solver->getBlockSize() != numblocks+1,get_out,"After setSize(), new num blocks was not in effect.");
}

} // namespace (anonymous)

int
main (int argc, char *argv[])
{
  using namespace Anasazi;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  typedef double                                  ScalarType;
  typedef ScalarTraits<ScalarType>::magnitudeType MagnitudeType;
  typedef Epetra_MultiVector                 MV;
  typedef Epetra_Operator                    OP;
  typedef MultiVecTraits<ScalarType,MV>     MVT;

#ifdef HAVE_MPI
  // Initialize MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  bool testFailed;
  bool verbose = false;
  bool debug = false;

  Teuchos::CommandLineProcessor cmdp (false, true);
  cmdp.setOption ("verbose", "quiet", &verbose, "Print messages and results.");
  cmdp.setOption ("debug", "nodebug", &debug, "Print debugging output from iteration.");
  if (cmdp.parse (argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
#ifdef HAVE_MPI
    MPI_Finalize ();
#endif
    return -1;
  }
  if (debug) verbose = true;

  // create the output manager
  int verbosity = Anasazi::Errors;
  if (debug) {
    verbosity += Anasazi::Debug;
  }
  RCP< OutputManager<ScalarType> > printer =
    rcp( new BasicOutputManager<ScalarType>( verbosity ) );

  printer->stream(Debug) << Anasazi_Version() << std::endl << std::endl;

  //  Problem information
  const int veclength = 99;
  int space_dim = 1;
  std::vector<double> brick_dim( space_dim );
  brick_dim[0] = 1.0;
  std::vector<int> elements( space_dim );
  elements[0] = veclength+1;

  // Create problem
  RCP<ModalProblem> testCase = rcp( new ModeLaplace1DQ1(Comm, brick_dim[0], elements[0]) );
  //
  // Get the stiffness and mass matrices
  RCP<const Epetra_CrsMatrix> K = rcp( const_cast<Epetra_CrsMatrix *>(testCase->getStiffness()), false );
  RCP<const Epetra_CrsMatrix> M = rcp( const_cast<Epetra_CrsMatrix *>(testCase->getMass()), false );
  //
  // Create the initial vectors
  // For BlockDavidson, this will dictate the initial size of the basis after a call to initialize()
  const int nev = 4;
  RCP<Epetra_MultiVector> ivec = rcp( new Epetra_MultiVector(K->OperatorDomainMap(), nev) );
  ivec->Random();
  //
  // Create eigenproblem: one standard and one generalized
  RCP<BasicEigenproblem<ScalarType,MV,OP> > probstd = rcp( new BasicEigenproblem<ScalarType, MV, OP>(K, ivec) );
  RCP<BasicEigenproblem<ScalarType,MV,OP> > probgen = rcp( new BasicEigenproblem<ScalarType, MV, OP>(K, M, ivec) );
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
#ifdef HAVE_MPI
    MPI_Finalize() ;
#endif
    return -1;
  }

  // create the orthogonalization managers: one standard and one M-based
  RCP< MatOrthoManager<ScalarType,MV,OP> > orthostd = rcp( new SVQBOrthoManager<ScalarType,MV,OP>() );
  RCP< MatOrthoManager<ScalarType,MV,OP> > orthogen = rcp( new SVQBOrthoManager<ScalarType,MV,OP>(M) );
  // create the sort manager
  RCP< SortManager<MagnitudeType> > sorter = rcp( new BasicSort<MagnitudeType>("LR") );
  // create the parameter list specifying blocksize > nev and full orthogonalization
  ParameterList pls;

  // begin testing
  testFailed = false;

  try
  {
    BlockDavidsonState<ScalarType,MV> istate;

    // try with default args
    if (verbose) printer->stream(Errors) << "Testing solver(default,default) with standard eigenproblem..." << std::endl;
    testsolver(probstd,printer,orthostd,sorter,pls,false,istate,false);
    if (verbose) printer->stream(Errors) << "Testing solver(default,default) with generalized eigenproblem..." << std::endl;
    testsolver(probgen,printer,orthogen,sorter,pls,false,istate,false);

    pls.set<int>("Num Blocks",4);

    // try with blocksize == getInitVec() size
    pls.set<int>("Block Size",nev);
    if (verbose) printer->stream(Errors) << "Testing solver(nev,4) with standard eigenproblem..." << std::endl;
    testsolver(probstd,printer,orthostd,sorter,pls,false,istate,false);
    TEUCHOS_TEST_FOR_EXCEPTION(pls.getEntryPtr("Block Size")->isUsed() == false, get_out, "Solver did not consume parameter \"Block Size\".");
    TEUCHOS_TEST_FOR_EXCEPTION(pls.getEntryPtr("Num Blocks")->isUsed() == false, get_out, "Solver did not consume parameter \"Num Blocks\".");
    if (verbose) printer->stream(Errors) << "Testing solver(nev,4) with generalized eigenproblem..." << std::endl;
    testsolver(probgen,printer,orthogen,sorter,pls,false,istate,false);

    // try with getInitVec() too small
    pls.set<int>("Block Size",2*nev);
    if (verbose) printer->stream(Errors) << "Testing solver(2*nev,4) with standard eigenproblem..." << std::endl;
    testsolver(probstd,printer,orthostd,sorter,pls,false,istate,false);
    if (verbose) printer->stream(Errors) << "Testing solver(2*nev,4) with generalized eigenproblem..." << std::endl;
    testsolver(probgen,printer,orthogen,sorter,pls,false,istate,false);

    // try with getInitVec() == two blocks
    pls.set<int>("Block Size",nev/2);
    if (verbose) printer->stream(Errors) << "Testing solver(nev/2,4) with standard eigenproblem..." << std::endl;
    testsolver(probstd,printer,orthostd,sorter,pls,false,istate,false);
    if (verbose) printer->stream(Errors) << "Testing solver(nev/2,4) with generalized eigenproblem..." << std::endl;
    testsolver(probgen,printer,orthogen,sorter,pls,false,istate,false);

    // try with a larger number of blocks; leave room for some expansion
    pls.set<int>("Block Size",nev);
    pls.set<int>("Num Blocks",15);
    if (verbose) printer->stream(Errors) << "Testing solver(nev,15) with standard eigenproblem..." << std::endl;
    testsolver(probstd,printer,orthostd,sorter,pls,false,istate,false);
    if (verbose) printer->stream(Errors) << "Testing solver(nev,15) with generalized eigenproblem..." << std::endl;
    testsolver(probgen,printer,orthogen,sorter,pls,false,istate,false);

    // try with a larger number of blocks+1
    pls.set<int>("Block Size",nev);
    pls.set<int>("Num Blocks",16);
    if (verbose) printer->stream(Errors) << "Testing solver(nev,16) with standard eigenproblem..." << std::endl;
    testsolver(probstd,printer,orthostd,sorter,pls,false,istate,false);
    if (verbose) printer->stream(Errors) << "Testing solver(nev,16) with generalized eigenproblem..." << std::endl;
    testsolver(probgen,printer,orthogen,sorter,pls,false,istate,false);

    // try with an invalid blocksize
    pls.set<int>("Block Size",0);
    pls.set<int>("Num Blocks",4);
    if (verbose) printer->stream(Errors) << "Testing solver(0,4) with standard eigenproblem..." << std::endl;
    testsolver(probstd,printer,orthostd,sorter,pls,true,istate,false);

    // try with an invalid numblocks
    pls.set<int>("Block Size",4);
    pls.set<int>("Num Blocks",1);
    if (verbose) printer->stream(Errors) << "Testing solver(4,1) with standard eigenproblem..." << std::endl;
    testsolver(probstd,printer,orthostd,sorter,pls,true,istate,false);

    // try with a too-large subspace
    pls.set<int>("Block Size",4);
    pls.set<int>("Num Blocks",veclength/4+1);
    if (verbose) printer->stream(Errors) << "Testing solver(4,toomany) with standard eigenproblem..." << std::endl;
    testsolver(probstd,printer,orthostd,sorter,pls,true,istate,false);

    // try with a non-Hermitian problem
    probstd->setHermitian(false);
    probstd->setProblem();
    if (verbose) printer->stream(Errors) << "Testing solver with non-Hermitian eigenproblem..." << std::endl;
    testsolver(probstd,printer,orthostd,sorter,pls,true,istate,false);
    // fix it now
    probstd->setHermitian(true);
    probstd->setProblem();

    // try with an unset problem
    // setHermitian will mark the problem as unset
    probstd->setHermitian(false);
    if (verbose) printer->stream(Errors) << "Testing solver with unset eigenproblem..." << std::endl;
    testsolver(probstd,printer,orthostd,sorter,pls,true,istate,false);
    // restore to hermitian, set problem
    probstd->setHermitian(true);
    probstd->setProblem();

    // try with a too-small initial basis
    if (verbose) printer->stream(Errors) << "Initializing solver with too-small basis..." << std::endl;
    pls.set("Block Size",4);
    pls.set("Num Blocks",2);
    istate.V  = MVT::Clone(*ivec,3);
    istate.KK = rcp( new SerialDenseMatrix<int,ScalarType>(3,3) );
    testsolver(probstd,printer,orthostd,sorter,pls,false,istate,true);

    // try with a non-integral initial basis
    if (verbose) printer->stream(Errors) << "Initializing solver with oddly sized basis..." << std::endl;
    pls.set("Block Size",4);
    pls.set("Num Blocks",2);
    istate.V  = MVT::Clone(*ivec,5);
    istate.KK = rcp( new SerialDenseMatrix<int,ScalarType>(5,5) );
    testsolver(probstd,printer,orthostd,sorter,pls,false,istate,true);

    // try with a inconsistent KK and V
    if (verbose) printer->stream(Errors) << "Initializing solver with inconsistent KK and V..." << std::endl;
    pls.set("Block Size",4);
    pls.set("Num Blocks",2);
    istate.V  = MVT::Clone(*ivec,4);
    istate.KK = rcp( new SerialDenseMatrix<int,ScalarType>(3,3) );
    testsolver(probstd,printer,orthostd,sorter,pls,false,istate,true);

    // create a dummy status tester
    RCP< StatusTest<ScalarType,MV,OP> > dumtester = rcp( new StatusTestMaxIters<ScalarType,MV,OP>(1) );

    // try with a null problem
    if (verbose) printer->stream(Errors) << "Testing solver with null eigenproblem..." << std::endl;
    try {
      RCP< BlockDavidson<ScalarType,MV,OP> > solver
        = rcp( new BlockDavidson<ScalarType,MV,OP>(Teuchos::null,sorter,printer,dumtester,orthostd,pls) );
      TEUCHOS_TEST_FOR_EXCEPTION(true,get_out,"Instantiating with invalid parameters failed to throw exception.");
    }
    catch (const std::invalid_argument &ia) {
      // caught expected exception
    }

    // try with a null sortman
    if (verbose) printer->stream(Errors) << "Testing solver with null sort manager..." << std::endl;
    try {
      RCP< BlockDavidson<ScalarType,MV,OP> > solver
        = rcp( new BlockDavidson<ScalarType,MV,OP>(probstd,Teuchos::null,printer,dumtester,orthostd,pls) );
      TEUCHOS_TEST_FOR_EXCEPTION(true,get_out,"Instantiating with invalid parameters failed to throw exception.");
    }
    catch (const std::invalid_argument &ia) {
      // caught expected exception
    }

    // try with a output man problem
    if (verbose) printer->stream(Errors) << "Testing solver with null output manager..." << std::endl;
    try {
      RCP< BlockDavidson<ScalarType,MV,OP> > solver
        = rcp( new BlockDavidson<ScalarType,MV,OP>(probstd,sorter,Teuchos::null,dumtester,orthostd,pls) );
      TEUCHOS_TEST_FOR_EXCEPTION(true,get_out,"Instantiating with invalid parameters failed to throw exception.");
    }
    catch (const std::invalid_argument &ia) {
      // caught expected exception
    }

    // try with a null status test
    if (verbose) printer->stream(Errors) << "Testing solver with null status test..." << std::endl;
    try {
      RCP< BlockDavidson<ScalarType,MV,OP> > solver
        = rcp( new BlockDavidson<ScalarType,MV,OP>(probstd,sorter,printer,Teuchos::null,orthostd,pls) );
      TEUCHOS_TEST_FOR_EXCEPTION(true,get_out,"Instantiating with invalid parameters failed to throw exception.");
    }
    catch (const std::invalid_argument &ia) {
      // caught expected exception
    }

    // try with a null orthoman
    if (verbose) printer->stream(Errors) << "Testing solver with null ortho manager..." << std::endl;
    try {
      RCP< BlockDavidson<ScalarType,MV,OP> > solver
        = rcp( new BlockDavidson<ScalarType,MV,OP>(probstd,sorter,printer,dumtester,Teuchos::null,pls) );
      TEUCHOS_TEST_FOR_EXCEPTION(true,get_out,"Instantiating with invalid parameters failed to throw exception.");
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

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

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
