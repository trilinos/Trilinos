// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER
//
// This test is for LOBPCG solving a generalized (Ax=Bxl) real Hermitian
// eigenvalue problem, using the LOBPCGSolMgr solver manager.
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
#include "AnasaziLOBPCGSolMgr.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif

#include "ModeLaplace1DQ1.h"

using namespace std;
using namespace Teuchos;


class get_out : public logic_error {
  public: get_out(const string &whatarg) : logic_error(whatarg) {}
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
  string filename("mhd1280b.cua");
  string which("LM");

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
  typedef Epetra_MultiVector                          MV;
  typedef Epetra_Operator                             OP;
  typedef Anasazi::MultiVecTraits<ScalarType,MV>     MVT;
  typedef Anasazi::OperatorTraits<ScalarType,MV,OP>  OPT;

  if (verbose && MyPID == 0) {
    cout << Anasazi::Anasazi_Version() << endl << endl;
  }

  //  Problem information
  int space_dim = 1;
  vector<double> brick_dim( space_dim );
  brick_dim[0] = 1.0;
  vector<int> elements( space_dim );
  elements[0] = 100;

  // Create problem
  RefCountPtr<ModalProblem> testCase = rcp( new ModeLaplace1DQ1(Comm, brick_dim[0], elements[0]) );
  //
  // Get the stiffness and mass matrices
  RefCountPtr<Epetra_CrsMatrix> K = rcp( const_cast<Epetra_CrsMatrix *>(testCase->getStiffness()), false );
  RefCountPtr<Epetra_CrsMatrix> M = rcp( const_cast<Epetra_CrsMatrix *>(testCase->getMass()), false );
  const int FIRST_BS = 5;
  const int SECOND_BS = 5;
  const int THIRD_BS = 5;
  const int NEV = FIRST_BS+SECOND_BS+THIRD_BS;


  ////////////////////////////////////////////////////////////////////////////////
  //
  // Shared parameters
  RefCountPtr<Anasazi::BasicEigenproblem<ScalarType,MV,OP> > problem;
  Anasazi::Eigensolution<ScalarType,MV> sol1, sol21, sol22, sol23;
  RefCountPtr<MV> cpoint, ev2 = rcp( new Epetra_MultiVector(K->OperatorDomainMap(), FIRST_BS+SECOND_BS+THIRD_BS) );
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
      RefCountPtr<Epetra_MultiVector> ivec = rcp( new Epetra_MultiVector(K->OperatorDomainMap(), FIRST_BS+SECOND_BS+THIRD_BS) );
      ivec->Random();
      problem = rcp( new Anasazi::BasicEigenproblem<ScalarType,MV,OP>(K,M,ivec) );
      problem->setHermitian(true);
      problem->setNEV( FIRST_BS+SECOND_BS+THIRD_BS );
      boolret = problem->setProblem();
      TEST_FOR_EXCEPTION(boolret != true,get_out,"Anasazi::BasicEigenproblem::SetProblem() returned with error.");
    }


    ////////////////////////////////////////////////////////////////////////////////
    //
    // Build the first solver manager
    {
      MyPL.set( "Block Size", FIRST_BS+SECOND_BS+THIRD_BS );
      Anasazi::LOBPCGSolMgr<ScalarType,MV,OP> solverman1(problem, MyPL);
      returnCode = solverman1.solve();
      TEST_FOR_EXCEPTION(returnCode != Anasazi::Converged, get_out, "First problem was not fully solved.");
      sol1 = problem->getSolution();
    }


    ////////////////////////////////////////////////////////////////////////////////
    //
    // Build the second/1 eigenproblem
    {
      RefCountPtr<Epetra_MultiVector> ivec = rcp( new Epetra_MultiVector(K->OperatorDomainMap(), FIRST_BS ) );
      ivec->Random();
      problem->setInitVec(ivec);
      problem->setNEV( FIRST_BS );
      boolret = problem->setProblem();
      TEST_FOR_EXCEPTION(boolret != true, get_out, "Anasazi::BasicEigenproblem::SetProblem() returned with error.");
    }


    ////////////////////////////////////////////////////////////////////////////////
    //
    // Build the second/1 solver manager
    {
      MyPL.set( "Block Size", FIRST_BS );
      Anasazi::LOBPCGSolMgr<ScalarType,MV,OP> solverman21(problem, MyPL);
      returnCode = solverman21.solve();
      TEST_FOR_EXCEPTION(returnCode != Anasazi::Converged, get_out, "Second/1 problem was not fully solved.");
      sol21 = problem->getSolution();
      vector<int> bsind1(FIRST_BS);
      for (int i=0; i<FIRST_BS; i++) bsind1[i] = i;
      MVT::SetBlock(*sol21.Evecs,bsind1,*ev2);
      cpoint = MVT::CloneView(*ev2,bsind1);
    }


    ////////////////////////////////////////////////////////////////////////////////
    //
    // Build the second/2 eigenproblem
    {
      RefCountPtr<Epetra_MultiVector> ivec = rcp( new Epetra_MultiVector(K->OperatorDomainMap(), SECOND_BS ) );
      ivec->Random();
      problem->setAuxVecs(cpoint);
      problem->setInitVec(ivec);
      problem->setNEV( SECOND_BS );
      boolret = problem->setProblem();
      TEST_FOR_EXCEPTION(boolret != true, get_out, "Anasazi::BasicEigenproblem::SetProblem() returned with error.");
    }


    ////////////////////////////////////////////////////////////////////////////////
    //
    // Build the second/2 solver manager
    {
      MyPL.set( "Block Size", SECOND_BS );
      Anasazi::LOBPCGSolMgr<ScalarType,MV,OP> solverman22(problem, MyPL);
      returnCode = solverman22.solve();
      TEST_FOR_EXCEPTION(returnCode != Anasazi::Converged, get_out, "Second/2 problem was not fully solved." );
      sol22 = problem->getSolution();
      vector<int> bsind2(SECOND_BS);
      for (int i=0; i<SECOND_BS; i++) bsind2[i] = FIRST_BS+i;
      MVT::SetBlock(*sol22.Evecs,bsind2,*ev2);
      vector<int> bsind12(FIRST_BS+SECOND_BS);
      for (int i=0; i<FIRST_BS+SECOND_BS; i++) bsind12[i] = i;
      cpoint = MVT::CloneView(*ev2,bsind12);
    }


    ////////////////////////////////////////////////////////////////////////////////
    //
    // Build the second/3 eigenproblem
    {
      RefCountPtr<Epetra_MultiVector> ivec = rcp( new Epetra_MultiVector(K->OperatorDomainMap(), THIRD_BS ) );
      ivec->Random();
      problem->setAuxVecs(cpoint);
      problem->setInitVec(ivec);
      problem->setNEV( THIRD_BS );
      boolret = problem->setProblem();
      TEST_FOR_EXCEPTION(boolret != true, get_out, "Anasazi::BasicEigenproblem::SetProblem() returned with error.");
    }


    ////////////////////////////////////////////////////////////////////////////////
    //
    // Build the second/3 solver manager
    {
      MyPL.set( "Block Size", THIRD_BS );
      Anasazi::LOBPCGSolMgr<ScalarType,MV,OP> solverman23(problem, MyPL);
      returnCode = solverman23.solve();
      TEST_FOR_EXCEPTION(returnCode != Anasazi::Converged, get_out, "Second/3 problem was not fully solved." );
      sol23 = problem->getSolution();
      vector<int> bsind3(THIRD_BS);
      for (int i=0; i<THIRD_BS; i++) bsind3[i] = FIRST_BS+SECOND_BS+i;
      MVT::SetBlock(*sol23.Evecs,bsind3,*ev2);
      cpoint = Teuchos::null;
    }
  }
  catch (get_out go) {
    if (verbose && MyPID==0) {
      cout << go.what() << endl;
    }
    testFailed = true;
  }

  if (testFailed == false) {
    cout.setf(ios::scientific, ios::floatfield);  
    cout.precision(6);
    //
    // check the checkpointed solution against the non-checkpointed solution
    //
    // First, check the eigenvalues
    //
    // unless we altogether missed an eigenvalue, then 
    // {sol21.Evals  sol22.Evals  sol23.Evals}  ==  {sol1.Evals}
    vector<Anasazi::Value<double> > Evals1 = sol1.Evals;
    vector<Anasazi::Value<double> > Evals2 = sol21.Evals;
    Evals2.insert(Evals2.end(),sol22.Evals.begin(),sol22.Evals.end());
    Evals2.insert(Evals2.end(),sol23.Evals.begin(),sol23.Evals.end());
    // compare the differences
    double maxd = 0;
    if (verbose && MyPID==0) {
      cout << setw(40) << "Computed Eigenvalues" << endl;
      cout << setw(20) << "Without c/p" << setw(20) << "With c/p" << setw(20) << "Rel. error" << endl;
      cout << "============================================================" << endl;
    }
    for (int i=0; i<NEV; i++) {
      double tmpd = SCT::magnitude((Evals1[i].realpart - Evals2[i].realpart)/Evals1[i].realpart);
      maxd = (tmpd > maxd ? tmpd : maxd);
      if (verbose && MyPID==0) {
        cout << setw(20) << Evals1[i].realpart << setw(20) << Evals2[i].realpart << setw(20) << tmpd << endl;
      }
    }
    if (maxd > tol) {
      testFailed = true;
    }
    if (verbose && MyPID==0) {
      cout << endl;
    }
    //
    // Second, check the eigenspaces
    RefCountPtr<MV> Mev1;
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
      cout << setw(20) << "|| EV1^H M EV2 - I ||_F" << endl;
      cout << setw(20) << maxd << endl;
      cout << endl;
    }
    if (maxd > SCT::squareroot(tol)*10) {
      testFailed = true;
    }
  }

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

  if (testFailed) {
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
