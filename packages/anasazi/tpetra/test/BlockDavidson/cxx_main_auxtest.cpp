// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
// This test is for BlockDavidson solving a generalized (Ax=Bxl) real Hermitian
// eigenvalue problem, using the BlockDavidsonSolMgr solver manager.
//
// This test checks that the auxiliary vector functionality, intended to support
// checkpointing, work as intended. The test will create an eigenproblem and solve
// for 15 eigenpairs. Then it will solve the same problem for 15 eigenpairs via checkpointing
// three groups of 5.
//
#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"

#include "AnasaziTpetraAdapter.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziBlockDavidsonSolMgr.hpp"
#include <Teuchos_CommandLineProcessor.hpp>

#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>

#include <ModeLaplace1DQ1.hpp>

using namespace Teuchos;

class get_out : public std::logic_error {
  public: get_out(const string &whatarg) : std::logic_error(whatarg) {}
};


int main(int argc, char *argv[]) 
{
  using std::cout;
  using std::endl;
  using std::ios_base;
  using std::setw;

  bool boolret;

  Tpetra::ScopeGuard tpetraScope (&argc,&argv);
  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm ();
  int MyPID = comm->getRank();

  bool testFailed = false;
  bool verbose = false;
  bool debug = false;
  string which("LM");
  bool printOnAllProcs = false;
  string outputfn = "";

  CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("debug","nodebug",&debug,"Print debugging information.");
  cmdp.setOption("sort",&which,"Targetted eigenvalues (SM or LM).");
  cmdp.setOption("allprint","oneprint",&printOnAllProcs,"Print output on all processors.");
  cmdp.setOption("outputfn",&outputfn,"File template for printing, ""%d"" substitues for processor rank, blank implies standard out.");
  if (cmdp.parse(argc,argv) != CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }
  if (debug) verbose = true;

  typedef double                              ST;
  typedef Teuchos::ScalarTraits<ST>          SCT;
  typedef SCT::magnitudeType                  MT;
  typedef Tpetra::MultiVector<ST>             MV;
  typedef MV::global_ordinal_type             GO;
  typedef MV::local_ordinal_type              LO;
  typedef MV::node_type                     Node;
  typedef Tpetra::Operator<ST>                OP;
  typedef Anasazi::MultiVecTraits<ST,MV>     MVT;
  typedef Anasazi::OperatorTraits<ST,MV,OP>  OPT;

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
  ModeLaplace1DQ1<ST,LO,GO,Node> testCase( comm, brick_dim[0], elements[0]);
  //
  // Get the stiffness and mass matrices
  RCP<const Tpetra::CrsMatrix<ST,LO,GO,Node> > K = testCase.getStiffness();
  RCP<const Tpetra::CrsMatrix<ST,LO,GO,Node> > M = testCase.getMass();
  //
  const int FIRST_BS = 5;
  const int SECOND_BS = 5;
  const int THIRD_BS = 5;
  const int NEV = FIRST_BS+SECOND_BS+THIRD_BS;


  ////////////////////////////////////////////////////////////////////////////////
  //
  // Shared parameters
  RCP<Anasazi::BasicEigenproblem<ST,MV,OP> > problem;
  Anasazi::Eigensolution<ST,MV> sol1, sol21, sol22, sol23;
  RCP<const MV> cpoint;
  RCP<MV> ev2 = rcp( new MV (K->getDomainMap(), FIRST_BS+SECOND_BS+THIRD_BS) );
  Anasazi::ReturnType returnCode;
  //
  // Verbosity level
  int verbosity = Anasazi::Errors + Anasazi::Warnings;
  if (debug) {
    verbosity += Anasazi::Debug;
  }
  //
  // Eigensolver parameters
  int numBlocks = 4;
  int maxRestarts = 250;
  MT tol = 1.0e-8;
  //
  // Create parameter list to pass into the solver managers
  ParameterList MyPL;
  MyPL.set( "Verbosity", verbosity );
  MyPL.set( "Which", which );
  MyPL.set( "Num Blocks", numBlocks );
  MyPL.set( "Maximum Restarts", maxRestarts );
  MyPL.set( "Convergence Tolerance", tol );
  MyPL.set( "Use Locking", true );
  MyPL.set( "Locking Tolerance", tol/10 );

  try {

    ////////////////////////////////////////////////////////////////////////////////
    //
    // Build the first eigenproblem
    {
      RCP<MV> ivec = rcp( new MV(K->getDomainMap(), FIRST_BS+SECOND_BS+THIRD_BS) );
      ivec->randomize();
      problem = rcp( new Anasazi::BasicEigenproblem<ST,MV,OP>(K,M,ivec) );
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
      Anasazi::BlockDavidsonSolMgr<ST,MV,OP> solverman1(problem, MyPL);
      returnCode = solverman1.solve();
      TEUCHOS_TEST_FOR_EXCEPTION(returnCode != Anasazi::Converged, get_out, "First problem was not fully solved.");
      sol1 = problem->getSolution();
    }


    ////////////////////////////////////////////////////////////////////////////////
    //
    // Build the second/1 eigenproblem
    {
      RCP<MV> ivec = rcp( new MV(K->getDomainMap(), FIRST_BS ) );
      ivec->randomize();
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
      Anasazi::BlockDavidsonSolMgr<ST,MV,OP> solverman21(problem, MyPL);
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
      RCP<MV> ivec = rcp( new MV(K->getDomainMap(), SECOND_BS ) );
      ivec->randomize();
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
      Anasazi::BlockDavidsonSolMgr<ST,MV,OP> solverman22(problem, MyPL);
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
      RCP<MV> ivec = rcp( new MV(K->getDomainMap(), THIRD_BS ) );
      ivec->randomize();
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
      Anasazi::BlockDavidsonSolMgr<ST,MV,OP> solverman23(problem, MyPL);
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
      cout << go.what() << endl;
    }
    testFailed = true;
  }

  if (testFailed == false) {
    cout.setf(std::ios::scientific, std::ios::floatfield);  
    cout.precision(6);
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
      cout << std::setw(40) << "Computed Eigenvalues" << endl;
      cout << std::setw(20) << "Without c/p" << std::setw(20) << "With c/p" << std::setw(20) << "Rel. error" << endl;
      cout << "============================================================" << endl;
    }
    for (int i=0; i<NEV; i++) {
      double tmpd = SCT::magnitude((Evals1[i].realpart - Evals2[i].realpart)/Evals1[i].realpart);
      maxd = (tmpd > maxd ? tmpd : maxd);
      if (verbose && MyPID==0) {
        cout << std::setw(20) << Evals1[i].realpart << std::setw(20) << Evals2[i].realpart << std::setw(20) << tmpd << endl;
      }
    }
    // a little fudge room
    if (maxd > tol*10) {
      testFailed = true;
    }
    if (verbose && MyPID==0) {
      cout << endl;
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
      cout << std::setw(20) << "|| EV1^H M EV2 - I ||_F" << endl;
      cout << std::setw(20) << maxd << endl;
      cout << endl;
    }
    if (maxd > SCT::squareroot(tol)*10) {
      testFailed = true;
    }
  }

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
