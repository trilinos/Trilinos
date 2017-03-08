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
// This test is for RTR solving a generalized (Ax=Bxl) real Hermitian
// eigenvalue problem, using the RTRSolMgr solver manager.
//
// This test checks that the auxiliary vector functionality, intended to support
// checkpointing, work as intended. The test will create an eigenproblem and solve
// for 15 eigenpairs. Then it will solve the same problem for 15 eigenpairs via checkpointing
// three groups of 5.
//
#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"

#include "AnasaziEpetraAdapter.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"

#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziRTRSolMgr.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif

#include "ModeLaplace1DQ1.h"

using namespace Teuchos;


class get_out : public std::logic_error {
  public: get_out(const std::string &whatarg) : std::logic_error(whatarg) {}
};


int main(int argc, char *argv[]) 
{
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
  std::string which("LR");
  bool skinny = true;

  bool success = true;
  try {

    CommandLineProcessor cmdp(false,true);
    cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
    cmdp.setOption("debug","nodebug",&debug,"Print debugging information.");
    cmdp.setOption("skinny","hefty",&skinny,"Use a skinny (low-mem) or hefty (higher-mem) implementation of IRTR.");
    cmdp.setOption("sort",&which,"Targetted eigenvalues (SR or LR).");
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
    typedef Epetra_MultiVector                          MV;
    typedef Epetra_Operator                             OP;
    typedef Anasazi::MultiVecTraits<ScalarType,MV>     MVT;
    typedef Anasazi::OperatorTraits<ScalarType,MV,OP>  OPT;

    if (verbose && MyPID == 0) {
      std::cout << Anasazi::Anasazi_Version() << std::endl << std::endl;
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
    RCP<const Epetra_CrsMatrix> K = rcp( const_cast<Epetra_CrsMatrix *>(testCase->getStiffness()), false );
    RCP<const Epetra_CrsMatrix> M = rcp( const_cast<Epetra_CrsMatrix *>(testCase->getMass()), false );
    const int FIRST_BS = 5;
    const int SECOND_BS = 5;
    const int THIRD_BS = 5;
    const int NEV = FIRST_BS+SECOND_BS+THIRD_BS;


    ////////////////////////////////////////////////////////////////////////////////
    //
    // Shared parameters
    RCP<Anasazi::BasicEigenproblem<ScalarType,MV,OP> > problem;
    Anasazi::Eigensolution<ScalarType,MV> sol1, sol21, sol22, sol23;
    RCP<const MV> cpoint;
    RCP<MV> ev2 = rcp( new Epetra_MultiVector(K->OperatorDomainMap(), FIRST_BS+SECOND_BS+THIRD_BS) );
    Anasazi::ReturnType returnCode;
    //
    // Verbosity level
    int verbosity = Anasazi::Errors + Anasazi::Warnings;
    if (debug) {
      verbosity += Anasazi::Debug;
    }
    //
    // Eigensolver parameters
    int maxIters = 450;
    MagnitudeType tol = 1.0e-8;
    //
    // Create parameter list to pass into the solver managers
    ParameterList MyPL;
    MyPL.set( "Skinny Solver", skinny);
    MyPL.set( "Verbosity", verbosity );
    MyPL.set( "Which", which );
    MyPL.set( "Maximum Iterations", maxIters );
    MyPL.set( "Convergence Tolerance", tol );
    MyPL.set( "Use Locking", true );
    MyPL.set( "Locking Tolerance", tol/10 );


    try {

      ////////////////////////////////////////////////////////////////////////////////
      //
      // Build the first eigenproblem
      {
        RCP<Epetra_MultiVector> ivec = rcp( new Epetra_MultiVector(K->OperatorDomainMap(), FIRST_BS+SECOND_BS+THIRD_BS) );
        ivec->Random();
        problem = rcp( new Anasazi::BasicEigenproblem<ScalarType,MV,OP>(K,M,ivec) );
        problem->setHermitian(true);
        problem->setNEV( FIRST_BS+SECOND_BS+THIRD_BS );
        boolret = problem->setProblem();
        TEUCHOS_TEST_FOR_EXCEPTION(boolret != true,get_out,"Anasazi::BasicEigenproblem::SetProblem() returned with error.");
      }


      ////////////////////////////////////////////////////////////////////////////////
      //
      // Build the first solver manager
      {
        MyPL.set( "Block Size", FIRST_BS+SECOND_BS+THIRD_BS );
        Anasazi::RTRSolMgr<ScalarType,MV,OP> solverman1(problem, MyPL);
        returnCode = solverman1.solve();
        TEUCHOS_TEST_FOR_EXCEPTION(returnCode != Anasazi::Converged, get_out, "First problem was not fully solved.");
        sol1 = problem->getSolution();
      }


      ////////////////////////////////////////////////////////////////////////////////
      //
      // Build the second/1 eigenproblem
      {
        RCP<Epetra_MultiVector> ivec = rcp( new Epetra_MultiVector(K->OperatorDomainMap(), FIRST_BS ) );
        ivec->Random();
        problem->setInitVec(ivec);
        problem->setNEV( FIRST_BS );
        boolret = problem->setProblem();
        TEUCHOS_TEST_FOR_EXCEPTION(boolret != true, get_out, "Anasazi::BasicEigenproblem::SetProblem() returned with error.");
      }


      ////////////////////////////////////////////////////////////////////////////////
      //
      // Build the second/1 solver manager
      {
        MyPL.set( "Block Size", FIRST_BS );
        Anasazi::RTRSolMgr<ScalarType,MV,OP> solverman21(problem, MyPL);
        returnCode = solverman21.solve();
        TEUCHOS_TEST_FOR_EXCEPTION(returnCode != Anasazi::Converged, get_out, "Second/1 problem was not fully solved.");
        sol21 = problem->getSolution();
        std::vector<int> bsind1(FIRST_BS);
        for (int i=0; i<FIRST_BS; i++) bsind1[i] = i;
        MVT::SetBlock(*sol21.Evecs,bsind1,*ev2);
        cpoint = MVT::CloneView(*ev2,bsind1);
      }


      ////////////////////////////////////////////////////////////////////////////////
      //
      // Build the second/2 eigenproblem
      {
        RCP<Epetra_MultiVector> ivec = rcp( new Epetra_MultiVector(K->OperatorDomainMap(), SECOND_BS ) );
        ivec->Random();
        problem->setAuxVecs(cpoint);
        problem->setInitVec(ivec);
        problem->setNEV( SECOND_BS );
        boolret = problem->setProblem();
        TEUCHOS_TEST_FOR_EXCEPTION(boolret != true, get_out, "Anasazi::BasicEigenproblem::SetProblem() returned with error.");
      }


      ////////////////////////////////////////////////////////////////////////////////
      //
      // Build the second/2 solver manager
      {
        MyPL.set( "Block Size", SECOND_BS );
        Anasazi::RTRSolMgr<ScalarType,MV,OP> solverman22(problem, MyPL);
        returnCode = solverman22.solve();
        TEUCHOS_TEST_FOR_EXCEPTION(returnCode != Anasazi::Converged, get_out, "Second/2 problem was not fully solved." );
        sol22 = problem->getSolution();
        std::vector<int> bsind2(SECOND_BS);
        for (int i=0; i<SECOND_BS; i++) bsind2[i] = FIRST_BS+i;
        MVT::SetBlock(*sol22.Evecs,bsind2,*ev2);
        std::vector<int> bsind12(FIRST_BS+SECOND_BS);
        for (int i=0; i<FIRST_BS+SECOND_BS; i++) bsind12[i] = i;
        cpoint = MVT::CloneView(*ev2,bsind12);
      }


      ////////////////////////////////////////////////////////////////////////////////
      //
      // Build the second/3 eigenproblem
      {
        RCP<Epetra_MultiVector> ivec = rcp( new Epetra_MultiVector(K->OperatorDomainMap(), THIRD_BS ) );
        ivec->Random();
        problem->setAuxVecs(cpoint);
        problem->setInitVec(ivec);
        problem->setNEV( THIRD_BS );
        boolret = problem->setProblem();
        TEUCHOS_TEST_FOR_EXCEPTION(boolret != true, get_out, "Anasazi::BasicEigenproblem::SetProblem() returned with error.");
      }


      ////////////////////////////////////////////////////////////////////////////////
      //
      // Build the second/3 solver manager
      {
        MyPL.set( "Block Size", THIRD_BS );
        Anasazi::RTRSolMgr<ScalarType,MV,OP> solverman23(problem, MyPL);
        returnCode = solverman23.solve();
        TEUCHOS_TEST_FOR_EXCEPTION(returnCode != Anasazi::Converged, get_out, "Second/3 problem was not fully solved." );
        sol23 = problem->getSolution();
        std::vector<int> bsind3(THIRD_BS);
        for (int i=0; i<THIRD_BS; i++) bsind3[i] = FIRST_BS+SECOND_BS+i;
        MVT::SetBlock(*sol23.Evecs,bsind3,*ev2);
        cpoint = Teuchos::null;
      }
    }
    catch (const get_out &go) {
      if (verbose && MyPID==0) {
        std::cout << go.what() << std::endl;
      }
      testFailed = true;
    }

    if (testFailed == false) {
      std::cout.setf(std::ios::scientific, std::ios::floatfield);  
      std::cout.precision(6);
      //
      // check the checkpointed solution against the non-checkpointed solution
      //
      // First, check the eigenvalues
      //
      // unless we altogether missed an eigenvalue, then 
      // {sol21.Evals  sol22.Evals  sol23.Evals}  ==  {sol1.Evals}
      std::vector<Anasazi::Value<double> > Evals1 = sol1.Evals;
      std::vector<Anasazi::Value<double> > Evals2 = sol21.Evals;
      Evals2.insert(Evals2.end(),sol22.Evals.begin(),sol22.Evals.end());
      Evals2.insert(Evals2.end(),sol23.Evals.begin(),sol23.Evals.end());
      // compare the differences
      double maxd = 0;
      if (verbose && MyPID==0) {
        std::cout << std::setw(40) << "Computed Eigenvalues" << std::endl;
        std::cout << std::setw(20) << "Without c/p" << std::setw(20) << "With c/p" << std::setw(20) << "Rel. error" << std::endl;
        std::cout << "============================================================" << std::endl;
      }
      for (int i=0; i<NEV; i++) {
        double tmpd = SCT::magnitude((Evals1[i].realpart - Evals2[i].realpart)/Evals1[i].realpart);
        maxd = (tmpd > maxd ? tmpd : maxd);
        if (verbose && MyPID==0) {
          std::cout << std::setw(20) << Evals1[i].realpart << std::setw(20) << Evals2[i].realpart << std::setw(20) << tmpd << std::endl;
        }
      }
      if (maxd > tol) {
        testFailed = true;
      }
      if (verbose && MyPID==0) {
        std::cout << std::endl;
      }
      //
      // Second, check the eigenspaces
      RCP<MV> Mev1;
      if (M != null) {
        Mev1 = MVT::Clone(*sol1.Evecs,NEV);
        OPT::Apply(*M,*sol1.Evecs,*Mev1);
      }
      else {
        Mev1 = sol1.Evecs;
      }
      SerialDenseMatrix<int,double> vtv(NEV,NEV);
      MVT::MvTransMv(1.0,*Mev1,*ev2,vtv);
      for (int i=0; i<NEV; i++) vtv(i,i) = SCT::magnitude(vtv(i,i)) - 1.0;
      maxd = vtv.normFrobenius();
      if (verbose && MyPID==0) {
        std::cout << std::setw(20) << "|| EV1^H M EV2 - I ||_F" << std::endl;
        std::cout << std::setw(20) << maxd << std::endl;
        std::cout << std::endl;
      }
      if (maxd > SCT::squareroot(tol)*10) {
        testFailed = true;
      }
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cout,success);

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

  if (testFailed || success==false) {
    if (verbose && MyPID==0) {
      std::cout << "End Result: TEST FAILED" << std::endl;	
    }
    return -1;
  }
  //
  // Default return value
  //
  if (verbose && MyPID==0) {
    std::cout << "End Result: TEST PASSED" << std::endl;
  }
  return 0;

}	
