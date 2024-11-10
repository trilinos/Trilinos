// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
//  This test is for the GeneralizedDavidson solver using the Thyra interface
//  The Thyra objects will be extracted from Epetra objects using the
//  Epetra-Thyra interface.
//  A few steps of BlockPCG is used as a preconditioner.
//
#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"

#include "AnasaziEpetraAdapter.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_InvOperator.h"
#include "Epetra_Vector.h"

#include "AnasaziBasicSort.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziGeneralizedDavidsonSolMgr.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif

#include "Teuchos_StandardCatchMacros.hpp"

#ifdef HAVE_EPETRA_THYRA
#include "AnasaziThyraAdapter.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#endif

#include "ModeLaplace1DQ1.h"
#include "BlockPCGSolver.h"

using namespace Teuchos;

int main(int argc, char *argv[])
{
  using std::cout;
  using std::endl;
  bool boolret;
  int MyPID;

#ifdef HAVE_MPI
  // Initialize MPI
  MPI_Init(&argc,&argv);
#endif

  bool success = false;
  bool verbose = false;
  try {
#ifdef HAVE_MPI
    // Initialize MPI
    Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm Comm;
#endif

    MyPID = Comm.MyPID();

    bool debug = false;
    std::string which("LM");

    CommandLineProcessor cmdp(false,true);
    cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
    cmdp.setOption("debug","nodebug",&debug,"Print debugging information.");
    cmdp.setOption("sort",&which,"Targetted eigenvalues (SM or LM).");
    if (cmdp.parse(argc,argv) != CommandLineProcessor::PARSE_SUCCESSFUL) {
#ifdef HAVE_MPI
      MPI_Finalize();
#endif
      return -1;
    }
    if (debug) verbose = true;

    typedef double ScalarType;
    typedef ScalarTraits<ScalarType>                   SCT;
    typedef SCT::magnitudeType               MagnitudeType;
    typedef Thyra::MultiVectorBase<ScalarType>          MV;
    typedef Thyra::LinearOpBase<ScalarType>             OP;
    typedef Anasazi::MultiVecTraits<ScalarType,MV>     MVT;
    typedef Anasazi::OperatorTraits<ScalarType,MV,OP>  OPT;
    const ScalarType ONE  = SCT::one();

    if (verbose && MyPID == 0) {
      cout << Anasazi::Anasazi_Version() << endl << endl;
    }

    //  Problem information
    int space_dim = 1;
    std::vector<double> brick_dim( space_dim );
    brick_dim[0] = 1.0;
    std::vector<int> elements( space_dim );
    elements[0] = 100;

    // Create problem
    RCP<ModalProblem> testCase = rcp( new ModeLaplace1DQ1(Comm, brick_dim[0], elements[0]) );
    //
    // Get the stiffness and mass matrices
    RCP<Epetra_CrsMatrix> K = rcp( const_cast<Epetra_CrsMatrix *>(testCase->getStiffness()), false );
    RCP<Epetra_CrsMatrix> M = rcp( const_cast<Epetra_CrsMatrix *>(testCase->getMass()), false );
    //
    // Create solver for mass matrix
    // Note that accuracy of Davidson solution does NOT depend on how accurately the BlockPCG is solved
    const int maxIterCG = 10;
    const double tolCG = 1e-2;

    RCP<BlockPCGSolver> opStiffness = rcp( new BlockPCGSolver(Comm, M.get(), tolCG, maxIterCG, 0) );
    opStiffness->setPreconditioner( 0 );
    RCP<Epetra_Operator> invStiffness = rcp( new Epetra_InvOperator(opStiffness.get()) );

    // Create the initial vectors
    int blockSize = 3;
    RCP<Epetra_MultiVector> ivec = rcp( new Epetra_MultiVector(K->OperatorDomainMap(), blockSize) );
    ivec->Random();

    // Get a pointer to the Epetra_Map
    RCP<const Epetra_Map> Map =
      rcp( &K->OperatorDomainMap(), false );

    // create a Thyra::VectorSpaceBase
    RCP<const Thyra::VectorSpaceBase<double> > epetra_vs =
      Thyra::create_VectorSpace(Map);

    // create a MultiVectorBase (from the Epetra_MultiVector)
    RCP<Thyra::MultiVectorBase<double> > thyra_ivec =
      Thyra::create_MultiVector(ivec, epetra_vs);

    // Create Thyra LinearOpBase objects from the Epetra_Operator objects
    RCP<const Thyra::LinearOpBase<double> > thyra_K =
      Thyra::epetraLinearOp(K);
    RCP<const Thyra::LinearOpBase<double> > thyra_M =
      Thyra::epetraLinearOp(M);
    RCP<const Thyra::LinearOpBase<double> > thyra_Prec =
      Thyra::epetraLinearOp(invStiffness);

    // Create eigenproblem
    const int nev = 5;
    RCP<Anasazi::BasicEigenproblem<ScalarType,MV,OP> > problem =
      rcp( new Anasazi::BasicEigenproblem<ScalarType,MV,OP>() );
    problem->setA(thyra_K);
    problem->setM(thyra_M);
    problem->setPrec(thyra_Prec);
    problem->setInitVec(thyra_ivec);
    //
    // Set the number of eigenvalues requested
    problem->setNEV( nev );
    //
    // Inform the eigenproblem that you are done passing it information
    boolret = problem->setProblem();
    if (boolret != true) {
      if (verbose && MyPID == 0) {
        cout << "Anasazi::BasicEigenproblem::SetProblem() returned with error." << endl
          << "End Result: TEST FAILED" << endl;
      }
#ifdef HAVE_MPI
      MPI_Finalize() ;
#endif
      return -1;
    }


    // Set verbosity level
    int verbosity = Anasazi::Errors + Anasazi::Warnings;
    if (verbose) {
      verbosity += Anasazi::FinalSummary + Anasazi::TimingDetails;
    }
    if (debug) {
      verbosity += Anasazi::Debug;
    }


    // Eigensolver parameters
    int maxRestarts = 25;
    int maxDim = 50;
    MagnitudeType tol = 1e-6;
    //
    // Create parameter list to pass into the solver manager
    ParameterList MyPL;
    MyPL.set( "Verbosity", verbosity );
    MyPL.set( "Which", which );
    MyPL.set( "Maximum Subspace Dimension", maxDim );
    MyPL.set( "Block Size", blockSize );
    MyPL.set( "Maximum Restarts", maxRestarts );
    MyPL.set( "Convergence Tolerance", tol );
    //
    // Create the solver manager
    Anasazi::GeneralizedDavidsonSolMgr<ScalarType,MV,OP> MySolverMgr(problem, MyPL);
    //
    // Check that the parameters were all consumed
    if (MyPL.getEntryPtr("Verbosity")->isUsed() == false ||
        MyPL.getEntryPtr("Which")->isUsed() == false ||
        MyPL.getEntryPtr("Maximum Subspace Dimension")->isUsed() == false ||
        MyPL.getEntryPtr("Block Size")->isUsed() == false ||
        MyPL.getEntryPtr("Maximum Restarts")->isUsed() == false ||
        MyPL.getEntryPtr("Convergence Tolerance")->isUsed() == false) {
      if (verbose && MyPID==0) {
        cout << "Failure! Unused parameters: " << endl;
        MyPL.unused(cout);
      }
    }


    // Solve the problem to the specified tolerances or length
    Anasazi::ReturnType returnCode = MySolverMgr.solve();
    success = (returnCode == Anasazi::Converged);

    // Get the eigenvalues and eigenvectors from the eigenproblem
    Anasazi::Eigensolution<ScalarType,MV> sol = problem->getSolution();
    std::vector<Anasazi::Value<ScalarType> > evals = sol.Evals;
    RCP<MV> evecs = sol.Evecs;
    int numev = sol.numVecs;

    if (numev > 0) {

      std::ostringstream os;
      os.setf(std::ios::scientific, std::ios::floatfield);
      os.precision(6);

      // Compute the direct residual
      std::vector<ScalarType> normV( numev );
      SerialDenseMatrix<int,ScalarType> T(numev,numev);
      for (int i=0; i<numev; i++) {
        T(i,i) = evals[i].realpart;
      }
      RCP<MV> Mvecs = MVT::Clone( *evecs, numev ),
        Kvecs = MVT::Clone( *evecs, numev );
      OPT::Apply( *thyra_K, *evecs, *Kvecs );
      OPT::Apply( *thyra_M, *evecs, *Mvecs );
      MVT::MvTimesMatAddMv( -ONE, *Mvecs, T, ONE, *Kvecs );
      // compute 2-norm of residuals
      std::vector<MagnitudeType> resnorm(numev);
      MVT::MvNorm( *Kvecs, resnorm );

      os << "Direct residual norms computed in GeneralizedDavidson_test.exe" << endl
        << std::setw(20) << "Eigenvalue" << std::setw(20) << "Residual" << endl
        << "----------------------------------------" << endl;
      for (int i=0; i<numev; i++) {
        os << std::setw(20) << evals[i].realpart << std::setw(20) << resnorm[i] << endl;
        success &= (resnorm[i] < tol);
      }
      if (verbose && MyPID==0) {
        cout << endl << os.str() << endl;
      }
    }

    if (verbose && MyPID==0) {
      if (success)
        cout << "End Result: TEST PASSED" << endl;
      else
        cout << "End Result: TEST FAILED" << endl;
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}
