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
//
//  This test is for the LOBPCG solver using the Thyra interface
//  The Thyra objects will be extracted from Epetra objects using the
//  Epetra-Thyra interface.
//  Therefore, this test should yield identical results compared against
//  the Epetra-only LOBPCG solver test.
//
#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"

#include "AnasaziEpetraAdapter.hpp"
#include "Epetra_CrsMatrix.h"

#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziLOBPCGSolMgr.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif

#ifdef HAVE_EPETRA_THYRA
#include "AnasaziThyraAdapter.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#endif

#include "ModeLaplace1DQ1.h"

using namespace Teuchos;

int main(int argc, char *argv[]) 
{

  using Teuchos::rcp_implicit_cast;
  using std::cout;
  using std::endl;

  bool boolret;
  int MyPID;

#ifdef HAVE_MPI
  // Initialize MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif
  MyPID = Comm.MyPID();

  bool testFailed = false;
  bool verbose = false;
  bool debug = false;
  std::string filename("mhd1280b.cua");
  std::string which("LM");

  bool success = true;
  try {

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

#ifndef HAVE_EPETRA_THYRA
    if (verbose && MyPid == 0) {
      cout << "Please configure Anasazi with:" << endl;
      cout << "--enable-epetra-thyra" << endl;
      cout << "--enable-anasazi-thyra" << endl;
    }
    return 0;
#endif

    typedef double ScalarType;
    typedef ScalarTraits<ScalarType>                   SCT;
    typedef SCT::magnitudeType               MagnitudeType;
    typedef Thyra::MultiVectorBase<double>              MV;
    typedef Thyra::LinearOpBase<double>                 OP;
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
    // Create the initial vectors
    int blockSize = 5;
    //
    // Get a pointer to the Epetra_Map
    RCP<const Epetra_Map> Map =  rcp( &K->OperatorDomainMap(), false );
    //
    // create an epetra multivector
    RCP<Epetra_MultiVector> ivec = 
      rcp( new Epetra_MultiVector(K->OperatorDomainMap(), blockSize) );
    ivec->Random();

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

    // Create eigenproblem
    const int nev = 5;
    RCP<Anasazi::BasicEigenproblem<ScalarType,MV,OP> > problem =
      rcp( new Anasazi::BasicEigenproblem<ScalarType,MV,OP>(thyra_K,thyra_M,thyra_ivec) );
    //
    // Inform the eigenproblem that the operator K is symmetric
    problem->setHermitian(true);
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
    int maxIters = 450;
    MagnitudeType tol = 1.0e-6;
    //
    // Create parameter list to pass into the solver manager
    ParameterList MyPL;
    MyPL.set( "Verbosity", verbosity );
    MyPL.set( "Which", which );
    MyPL.set( "Block Size", blockSize );
    MyPL.set( "Maximum Iterations", maxIters );
    MyPL.set( "Convergence Tolerance", tol );
    MyPL.set( "Use Locking", true );
    MyPL.set( "Locking Tolerance", tol/10 );
    MyPL.set( "Full Ortho", true );
    //
    // Create the solver manager
    Anasazi::LOBPCGSolMgr<ScalarType,MV,OP> MySolverMan(problem, MyPL);
    // 
    // Check that the parameters were all consumed
    if (MyPL.getEntryPtr("Verbosity")->isUsed() == false ||
        MyPL.getEntryPtr("Which")->isUsed() == false ||
        MyPL.getEntryPtr("Block Size")->isUsed() == false ||
        MyPL.getEntryPtr("Maximum Iterations")->isUsed() == false ||
        MyPL.getEntryPtr("Convergence Tolerance")->isUsed() == false ||
        MyPL.getEntryPtr("Use Locking")->isUsed() == false ||
        MyPL.getEntryPtr("Locking Tolerance")->isUsed() == false ||
        MyPL.getEntryPtr("Full Ortho")->isUsed() == false) {
      if (verbose && MyPID==0) {
        cout << "Failure! Unused parameters: " << endl;
        MyPL.unused(cout);
      }
    }


    // Solve the problem to the specified tolerances or length
    Anasazi::ReturnType returnCode = MySolverMan.solve();
    if (returnCode != Anasazi::Converged) {
      testFailed = true;
    }

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
      // compute M-norm of residuals
      OPT::Apply( *thyra_M, *Kvecs, *Mvecs );
      MVT::MvDot( *Mvecs, *Kvecs, normV );

      os << "Direct residual norms computed in LOBPCGThyra_test.exe" << endl
        << std::setw(20) << "Eigenvalue" << std::setw(20) << "Residual(M)" << endl
        << "----------------------------------------" << endl;
      for (int i=0; i<numev; i++) {
        if ( SCT::magnitude(evals[i].realpart) != SCT::zero() ) {
          normV[i] = SCT::magnitude( SCT::squareroot( normV[i] ) / evals[i].realpart );
        }
        else {
          normV[i] = SCT::magnitude( SCT::squareroot( normV[i] ) );
        }
        os << std::setw(20) << evals[i].realpart << std::setw(20) << normV[i] << endl;
        if ( normV[i] > tol ) {
          testFailed = true;
        }
      }
      if (verbose && MyPID==0) {
        cout << endl << os.str() << endl;
      }

    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,cout,success);

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

  if (testFailed || success==false) {
    if (verbose && MyPID==0) {
      cout << "End Result: TEST FAILED" << endl;	
    }
    return -1;
  }
  //
  // Default return value
  //
  if (verbose && MyPID==0) {
    cout << "End Result: TEST PASSED" << endl;
  }
  return 0;

}
