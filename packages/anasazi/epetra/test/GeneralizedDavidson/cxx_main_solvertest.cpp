// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
//  This test is for the generalized Davidsoneigensolver
//
#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"

#include "AnasaziEpetraAdapter.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"

#include "AnasaziGeneralizedDavidson.hpp"

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

using namespace Teuchos;
using namespace Anasazi;

typedef double                            ScalarType;
typedef ScalarTraits<ScalarType>          ScalarTypeTraits;
typedef ScalarTypeTraits::magnitudeType   MagnitudeType;
typedef Epetra_MultiVector                MV;
typedef Epetra_Operator                   OP;
typedef MultiVecTraits<ScalarType,MV>     MVTraits;
typedef OperatorTraits<ScalarType,MV,OP>  OpTraits;

class get_out : public std::logic_error {
  public: get_out(const std::string &whatarg) : std::logic_error(whatarg) {}
};

void checks( RCP<GeneralizedDavidson<ScalarType,MV,OP> > solver, int blocksize, int maxdim,
             RCP<Eigenproblem<ScalarType,MV,OP> > problem,
             RCP<MatOrthoManager<ScalarType,MV,OP> > ortho) {
  GeneralizedDavidsonState<ScalarType,MV> state = solver->getState();

  TEUCHOS_TEST_FOR_EXCEPTION(MVTraits::GetNumberVecs(*state.V) != solver->getMaxSubspaceDim(),get_out,"getMaxSubspaceDim() does not match allocated size for V");

  TEUCHOS_TEST_FOR_EXCEPTION(&solver->getProblem() != problem.get(),get_out,"getProblem() did not return the submitted problem.");

  if (solver->isInitialized())
  {
      // Generalized Davidson block size is equal to or one greater than user specified block size
      // Because GeneralizedDavidson block size is variable, this check only applied to an initialized solver
      TEUCHOS_TEST_FOR_EXCEPTION(solver->getBlockSize() != blocksize && solver->getBlockSize() != blocksize+1, get_out,"Solver block size does not match specified block size.");

    std::vector<Anasazi::Value<ScalarType> > ritzValues = solver->getRitzValues();

    // check Ritz residuals
    std::vector<MagnitudeType> ritzResids = solver->getRitzRes2Norms();

    // get Ritz index
    std::vector<int> ritzIndex = solver->getRitzIndex();

    // get Ritz vector
    RCP<const MV> ritzVectors = solver->getRitzVectors();

    int numRitzVecs = MVTraits::GetNumberVecs(*ritzVectors);

    RCP<MV> tmpVecs = MVTraits::Clone( *ritzVectors, numRitzVecs );

    // Compute Ritz residuals like R = A*X - B*X*T
    Teuchos::SerialDenseMatrix<int,ScalarType> T(numRitzVecs,numRitzVecs);
    Teuchos::RCP<MV> ritzResiduals = MVTraits::Clone( *ritzVectors, numRitzVecs );
    for (int i=0; i<T.numRows(); i++) T(i,i) = ritzValues[i].realpart;
    OpTraits::Apply( *(problem->getA()), *ritzVectors, *ritzResiduals );
    if( problem->getM() != Teuchos::null )
    {
        OpTraits::Apply( *(problem->getM()), *ritzVectors, *tmpVecs );
    }
    else
    {
        std::vector<int> inds(numRitzVecs);
        for( int i=0; i<numRitzVecs; ++i ) inds[i]=i;
        MVTraits::SetBlock( *ritzVectors, inds, *tmpVecs );
    }
    MVTraits::MvTimesMatAddMv(-1.0,*tmpVecs,T,1.0,*ritzResiduals);

    // Compute the norm of the Ritz residual vectors
    std::vector<MagnitudeType> ritzVecNrm( numRitzVecs );
    MVTraits::MvNorm( *ritzVectors, ritzVecNrm );
    MagnitudeType error;
    for (int i=0; i<numRitzVecs; i++) {
      error = Teuchos::ScalarTraits<MagnitudeType>::magnitude( ritzVecNrm[i] - 1.0 );
      TEUCHOS_TEST_FOR_EXCEPTION(error > 1e-14,get_out,"Ritz vectors are not normalized.");
    }

    std::vector<MagnitudeType> ritzResNrm( MVTraits::GetNumberVecs( *ritzResiduals ) );
    MVTraits::MvNorm( *ritzResiduals, ritzResNrm );
    for (int i=0; i<(int)ritzResNrm.size(); i++) {
      error = Teuchos::ScalarTraits<MagnitudeType>::magnitude( ritzResids[i] - ritzResNrm[i] );
      TEUCHOS_TEST_FOR_EXCEPTION(error > 1e-12,get_out,"Ritz residuals from iteration do not compare to those computed.");
    }
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
                 GeneralizedDavidsonState<ScalarType,MV> initstate, bool invalidinit)
{
  // create a status tester
  RCP< StatusTest<ScalarType,MV,OP> > tester = rcp( new StatusTestMaxIters<ScalarType,MV,OP>(1) );

  // create the solver
  RCP< GeneralizedDavidson<ScalarType,MV,OP> > solver;
  try {
    solver = rcp( new GeneralizedDavidson<ScalarType,MV,OP>(problem,sorter,printer,tester,ortho,pls) );
    TEUCHOS_TEST_FOR_EXCEPTION(invalid, get_out, "Instantiating with invalid parameters failed to throw exception.")
  }
  catch (const std::invalid_argument &ia) {
    if (!invalid) {
      printer->stream(Warnings) << "Error thrown at instantiation: " << ia.what() << std::endl;
    }
    TEUCHOS_TEST_FOR_EXCEPTION(!invalid, get_out, "Instantiating with valid parameters unexpectadly threw exception.");

    // caught expected exception
    return;
  }

  const int  blocksize = pls.get<int>("Block Size");
  const int  maxdim = pls.get<int>("Maximum Subspace Dimension");

  SolverUtils<ScalarType,MV,OP> msutils;

  // solver should be uninitialized
  TEUCHOS_TEST_FOR_EXCEPTION(solver->isInitialized() != false,get_out,"Solver should be un-initialized after instantiation.");
  TEUCHOS_TEST_FOR_EXCEPTION(solver->getNumIters() != 0,get_out,"Number of iterations after initialization should be zero after init.")
  TEUCHOS_TEST_FOR_EXCEPTION(solver->getAuxVecs().size() != 0,get_out,"getAuxVecs() should return empty.");
  checks(solver,blocksize,maxdim,problem,ortho);

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
  TEUCHOS_TEST_FOR_EXCEPTION(solver->getCurSubspaceDim() != blocksize,get_out,"after init, getCurSubspaceDim() should be equal to block size.");
  checks(solver,blocksize,maxdim,problem,ortho);

  // call iterate(); solver should perform exactly one iteration and return; status test should be passed
  solver->iterate();
  TEUCHOS_TEST_FOR_EXCEPTION(tester->getStatus() != Passed,get_out,"Solver returned from iterate() but getStatus() not Passed.");
  TEUCHOS_TEST_FOR_EXCEPTION(solver->isInitialized() != true,get_out,"Solver should be initialized after call to initialize().");
  TEUCHOS_TEST_FOR_EXCEPTION(solver->getNumIters() != 1,get_out,"Number of iterations should be one.")
  TEUCHOS_TEST_FOR_EXCEPTION(solver->getAuxVecs().size() != 0,get_out,"getAuxVecs() should return empty.");
  TEUCHOS_TEST_FOR_EXCEPTION(solver->getCurSubspaceDim() != 2*blocksize,get_out,"after one step, getCurSubspaceDim() should be 2*blocksize.");
  checks(solver,blocksize,maxdim,problem,ortho);

  // reset numiters, call iterate(); solver should perform exactly one iteration and return; status test should be passed
  solver->resetNumIters();
  TEUCHOS_TEST_FOR_EXCEPTION(solver->getNumIters() != 0,get_out,"Number of iterations should be zero after resetNumIters().")
  solver->iterate();
  TEUCHOS_TEST_FOR_EXCEPTION(tester->getStatus() != Passed,get_out,"Solver returned from iterate() but getStatus() not Passed.");
  TEUCHOS_TEST_FOR_EXCEPTION(solver->isInitialized() != true,get_out,"Solver should be initialized after call to initialize().");
  TEUCHOS_TEST_FOR_EXCEPTION(solver->getNumIters() != 0,get_out,"Number of iterations should be zero.")
  TEUCHOS_TEST_FOR_EXCEPTION(solver->getAuxVecs().size() != 0,get_out,"getAuxVecs() should return empty.");
  TEUCHOS_TEST_FOR_EXCEPTION(solver->getCurSubspaceDim() != 2*blocksize,get_out,"after two steps, getCurSubspaceDim() should be 2*blocksize.");
  checks(solver,blocksize,maxdim,problem,ortho);
}

int main(int argc, char *argv[])
{

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

  CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("debug","nodebug",&debug,"Print debugging output from iteration.");
  if (cmdp.parse(argc,argv) != CommandLineProcessor::PARSE_SUCCESSFUL) {
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return -1;
  }
  if (debug) verbose = true;

  // create the output manager
  int verbosity = Anasazi::Errors;
  if (verbose) {
    verbosity += Anasazi::Warnings;
  }
  if (debug) {
    verbosity += Anasazi::Debug;
  }
  RCP< OutputManager<ScalarType> > printer =
    rcp( new BasicOutputManager<ScalarType>( verbosity ) );

  printer->stream(Debug) << Anasazi_Version() << std::endl;

  //  Problem information
  int space_dim = 1;
  std::vector<double> brick_dim( space_dim );
  brick_dim[0] = 1.0;
  std::vector<int> elements( space_dim );
  elements[0] = 100+1;

  // Create problem
  RCP<ModalProblem> testCase = rcp( new ModeLaplace1DQ1(Comm, brick_dim[0], elements[0]) );
  //
  // Get the stiffness and mass matrices
  RCP<const Epetra_CrsMatrix> K = rcp( const_cast<Epetra_CrsMatrix *>(testCase->getStiffness()), false );
  RCP<const Epetra_CrsMatrix> M = rcp( const_cast<Epetra_CrsMatrix *>(testCase->getMass()), false );

  //
  // Create the initial vectors
  const int nev = 5;
  RCP<Epetra_MultiVector> ivec = rcp( new Epetra_MultiVector(K->OperatorDomainMap(), nev) );
  ivec->Random();
  //
  // Create eigenproblem: one standard and one generalized
  RCP<BasicEigenproblem<ScalarType,MV,OP> > probstd = rcp( new BasicEigenproblem<ScalarType, MV, OP>() );
  probstd->setA(K);
  probstd->setInitVec(ivec);
  RCP<BasicEigenproblem<ScalarType,MV,OP> > probgen = rcp( new BasicEigenproblem<ScalarType, MV, OP>() );
  probgen->setA(K);
  probgen->setM(M);
  probgen->setInitVec(ivec);
  //
  // Inform the eigenproblem that the operator A is not symmetric (even though it is)
  probstd->setHermitian(false);
  probgen->setHermitian(false);
  //
  // Set the number of eigenvalues requested
  probstd->setNEV( nev );
  probgen->setNEV( nev );
  //
  // Inform the eigenproblem that you are finishing passing it information
  if ( probstd->setProblem() != true || probgen->setProblem() != true ) {
    printer->stream(Warnings) << "Anasazi::BasicEigenproblem::SetProblem() returned with error." << std::endl
      << "End Result: TEST FAILED" << std::endl;
#ifdef HAVE_MPI
    MPI_Finalize() ;
#endif
    return -1;
  }

  // create the orthogonalization managers: one standard and one M-based
  RCP< MatOrthoManager<ScalarType,MV,OP> > orthostd = rcp( new SVQBOrthoManager<ScalarType,MV,OP>() );
  RCP< MatOrthoManager<ScalarType,MV,OP> > orthogen = rcp( new SVQBOrthoManager<ScalarType,MV,OP>() );
  // create the sort manager
  RCP< SortManager<MagnitudeType> > sorter = rcp( new BasicSort<MagnitudeType>("LM") );
  // create the parameter list specifying blocksize > nev and full orthogonalization
  ParameterList pls;

  // begin testing
  testFailed = false;

  try
  {
    GeneralizedDavidsonState<ScalarType,MV> istate;

    pls.set<int>("Block Size",nev);
    pls.set<int>("Maximum Subspace Dimension",3*nev);
    pls.set<int>("Number of Ritz Vectors",nev);
    printer->stream(Warnings) << "Testing solver(nev,3*nev) with standard eigenproblem..." << std::endl;
    testsolver(probstd,printer,orthostd,sorter,pls,false,istate,false);
    pls.set<int>("Maximum Subspace Dimension",3*nev);
    printer->stream(Warnings) << "Testing solver(nev,3*nev) with generalized eigenproblem..." << std::endl;
    testsolver(probgen,printer,orthogen,sorter,pls,false,istate,false);

    pls.set<int>("Block Size",nev);
    pls.set<int>("Maximum Subspace Dimension",4*nev);
    printer->stream(Warnings) << "Testing solver(nev,4*nev) with standard eigenproblem..." << std::endl;
    testsolver(probstd,printer,orthostd,sorter,pls,false,istate,false);
    pls.set<int>("Maximum Subspace Dimension",4*nev);
    printer->stream(Warnings) << "Testing solver(nev,4*nev) with generalized eigenproblem..." << std::endl;
    testsolver(probgen,printer,orthogen,sorter,pls,false,istate,false);

    pls.set<int>("Block Size",2*nev);
    pls.set<int>("Maximum Subspace Dimension",4*nev);
    printer->stream(Warnings) << "Testing solver(2*nev,4*nev) with standard eigenproblem..." << std::endl;
    testsolver(probstd,printer,orthostd,sorter,pls,false,istate,false);
    pls.set<int>("Maximum Subspace Dimension",4*nev);
    printer->stream(Warnings) << "Testing solver(2*nev,4*nev) with generalized eigenproblem..." << std::endl;
    testsolver(probgen,printer,orthogen,sorter,pls,false,istate,false);

    pls.set<int>("Block Size",nev/2);
    pls.set<int>("Maximum Subspace Dimension",3*nev);
    printer->stream(Warnings) << "Testing solver(nev/2,3*nev) with standard eigenproblem..." << std::endl;
    testsolver(probstd,printer,orthostd,sorter,pls,false,istate,false);
    pls.set<int>("Maximum Subspace Dimension",3*nev);
    printer->stream(Warnings) << "Testing solver(nev/2,3*nev) with generalized eigenproblem..." << std::endl;
    testsolver(probgen,printer,orthogen,sorter,pls,false,istate,false);

    pls.set<int>("Block Size",nev/2);
    pls.set<int>("Maximum Subspace Dimension",4*nev);
    printer->stream(Warnings) << "Testing solver(nev/2,4*nev) with standard eigenproblem..." << std::endl;
    testsolver(probstd,printer,orthostd,sorter,pls,false,istate,false);
    pls.set<int>("Maximum Subspace Dimension",4*nev);
    printer->stream(Warnings) << "Testing solver(nev/2,4*nev) with generalized eigenproblem..." << std::endl;
    testsolver(probgen,printer,orthogen,sorter,pls,false,istate,false);

    // try with an invalid blocksize
    pls.set<int>("Block Size",0);
    pls.set<int>("Maximum Subspace Dimension",4*nev);
    printer->stream(Warnings) << "Testing solver(0,4*nev) with standard eigenproblem..." << std::endl;
    testsolver(probstd,printer,orthostd,sorter,pls,true,istate,false);

    // try with an invalid maxdim
    pls.set<int>("Block Size",nev);
    pls.set<int>("Maximum Subspace Dimension",0);
    printer->stream(Warnings) << "Testing solver(nev,0) with standard eigenproblem..." << std::endl;
    testsolver(probstd,printer,orthostd,sorter,pls,true,istate,false);

    // try with an invalid maxdim: invalid because it must be greater than nev
    pls.set<int>("Block Size",nev);
    pls.set<int>("Maximum Subspace Dimension",nev);
    printer->stream(Warnings) << "Testing solver(nev,nev) with standard eigenproblem..." << std::endl;
    testsolver(probstd,printer,orthostd,sorter,pls,true,istate,false);

    // try with a too-large subspace
    probstd->setHermitian(true);
    probstd->setProblem();
    pls.set<int>("Maximum Subspace Dimension",100+1);
    printer->stream(Warnings) << "Testing solver(4,toomany,Hermitian) with standard eigenproblem..." << std::endl;
    testsolver(probstd,printer,orthostd,sorter,pls,true,istate,false);

    // try with an unset problem
    // setHermitian will mark the problem as unset
    probstd->setHermitian(true);
    printer->stream(Warnings) << "Testing solver with unset eigenproblem..." << std::endl;
    testsolver(probstd,printer,orthostd,sorter,pls,true,istate,false);
    // set problem: now hermitian
    probstd->setProblem();

    // create a dummy status tester
    RCP< StatusTest<ScalarType,MV,OP> > dumtester = rcp( new StatusTestMaxIters<ScalarType,MV,OP>(1) );

    // try with a null problem
    printer->stream(Warnings) << "Testing solver with null eigenproblem..." << std::endl;
    try {
      RCP< GeneralizedDavidson<ScalarType,MV,OP> > solver
        = rcp( new GeneralizedDavidson<ScalarType,MV,OP>(Teuchos::null,sorter,printer,dumtester,orthostd,pls) );
      TEUCHOS_TEST_FOR_EXCEPTION(true,get_out,"Instantiating with invalid parameters failed to throw exception.");
    }
    catch (const std::invalid_argument &ia) {
      // caught expected exception
    }

    // try with a null sortman
    printer->stream(Warnings) << "Testing solver with null sort manager..." << std::endl;
    try {
      RCP< GeneralizedDavidson<ScalarType,MV,OP> > solver
        = rcp( new GeneralizedDavidson<ScalarType,MV,OP>(probstd,Teuchos::null,printer,dumtester,orthostd,pls) );
      TEUCHOS_TEST_FOR_EXCEPTION(true,get_out,"Instantiating with invalid parameters failed to throw exception.");
    }
    catch (const std::invalid_argument &ia) {
      // caught expected exception
    }

    // try with a output man problem
    printer->stream(Warnings) << "Testing solver with null output manager..." << std::endl;
    try {
      RCP< GeneralizedDavidson<ScalarType,MV,OP> > solver
        = rcp( new GeneralizedDavidson<ScalarType,MV,OP>(probstd,sorter,Teuchos::null,dumtester,orthostd,pls) );
      TEUCHOS_TEST_FOR_EXCEPTION(true,get_out,"Instantiating with invalid parameters failed to throw exception.");
    }
    catch (const std::invalid_argument &ia) {
      // caught expected exception
    }

    // try with a null status test
    printer->stream(Warnings) << "Testing solver with null status test..." << std::endl;
    try {
      RCP< GeneralizedDavidson<ScalarType,MV,OP> > solver
        = rcp( new GeneralizedDavidson<ScalarType,MV,OP>(probstd,sorter,printer,Teuchos::null,orthostd,pls) );
      TEUCHOS_TEST_FOR_EXCEPTION(true,get_out,"Instantiating with invalid parameters failed to throw exception.");
    }
    catch (const std::invalid_argument &ia) {
      // caught expected exception
    }

    // try with a null orthoman
    printer->stream(Warnings) << "Testing solver with null ortho manager..." << std::endl;
    try {
      RCP< GeneralizedDavidson<ScalarType,MV,OP> > solver
        = rcp( new GeneralizedDavidson<ScalarType,MV,OP>(probstd,sorter,printer,dumtester,Teuchos::null,pls) );
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
    printer->stream(Warnings) << std::endl << "End Result: TEST FAILED" << std::endl;
    return -1;
  }
  //
  // Default return value
  //
  printer->stream(Warnings) << std::endl << "End Result: TEST PASSED" << std::endl;
  return 0;

}
