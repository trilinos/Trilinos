// @HEADER
// ***********************************************************************
//
//                 Belos: Block Linear Solvers Package
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
// This driver reads a problem from a Harwell-Boeing (HB) file.
// Multiple right-hand-sides are created randomly.
// The initial guesses are all set to zero. 
//
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosOutputManager.hpp"
#include "BelosStatusTestMaxIters.hpp"
#include "BelosStatusTestMaxRestarts.hpp"
#include "BelosStatusTestResNorm.hpp"
#include "BelosStatusTestOutputter.hpp"
#include "BelosStatusTestCombo.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosPseudoBlockGmres.hpp"
#include "createEpetraProblem.hpp"
#include "Ifpack_IlukGraph.h"
#include "Ifpack_CrsRiluk.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_ParameterList.hpp"

int main(int argc, char *argv[]) {
  //
#ifdef EPETRA_MPI	
  MPI_Init(&argc,&argv);
  Belos::MPIFinalize mpiFinalize; // Will call finalize with *any* return
#endif	
  //
  typedef double                            ST;
  typedef Teuchos::ScalarTraits<ST>        SCT;
  typedef SCT::magnitudeType                MT;
  typedef Belos::MultiVec<ST>               MV;
  typedef Belos::Operator<ST>               OP;
  typedef Belos::MultiVecTraits<ST,MV>     MVT;
  typedef Belos::OperatorTraits<ST,MV,OP>  OPT;

  using Teuchos::ParameterList;
  using Teuchos::RefCountPtr;
  using Teuchos::rcp;

  bool verbose = false, proc_verbose = false;
  int frequency = -1;  // how often residuals are printed by solver
  int numrhs = 15;
  int blocksize = 3; // number of linear systems solved at one time
  int numrestarts = 15; // number of restarts allowed 
  int length = 25;
  std::string filename("orsirr1.hb");
  MT tol = 1.0e-5;  // relative residual tolerance

  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("frequency",&frequency,"Solvers frequency for printing residuals (#iters).");
  cmdp.setOption("filename",&filename,"Filename for Harwell-Boeing test matrix.");
  cmdp.setOption("tol",&tol,"Relative residual tolerance used by GMRES solver.");
  cmdp.setOption("num-rhs",&numrhs,"Number of right-hand sides to be solved for.");
  cmdp.setOption("blocksize",&blocksize,"Number of linear systems solved together.");
  cmdp.setOption("num-restarts",&numrestarts,"Number of restarts allowed for GMRES solver.");
  cmdp.setOption("subspace-size",&length,"Dimension of Krylov subspace used by GMRES.");  
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }
  if (!verbose)
    frequency = -1;  // reset frequency if test is not verbose
  //
  // Get the problem
  //
  int MyPID;
  RefCountPtr<Epetra_CrsMatrix> A;
  int return_val =Belos::createEpetraProblem(filename,NULL,&A,NULL,NULL,&MyPID);
  const Epetra_Map &Map = A->RowMap();
  if(return_val != 0) return return_val;
  proc_verbose = verbose && (MyPID==0); /* Only print on zero processor */
  //
  // *****Construct the Preconditioner*****
  //
  if (proc_verbose) cout << endl << endl;
  if (proc_verbose) cout << "Constructing ILU preconditioner" << endl;
  int Lfill = 2;
  // if (argc > 2) Lfill = atoi(argv[2]);
  if (proc_verbose) cout << "Using Lfill = " << Lfill << endl;
  int Overlap = 2;
  // if (argc > 3) Overlap = atoi(argv[3]);
  if (proc_verbose) cout << "Using Level Overlap = " << Overlap << endl;
  double Athresh = 0.0;
  // if (argc > 4) Athresh = atof(argv[4]);
  if (proc_verbose) cout << "Using Absolute Threshold Value of " << Athresh << endl;
  double Rthresh = 1.0;
  // if (argc >5) Rthresh = atof(argv[5]);
  if (proc_verbose) cout << "Using Relative Threshold Value of " << Rthresh << endl;
  //
  Teuchos::RefCountPtr<Ifpack_IlukGraph> ilukGraph;
  Teuchos::RefCountPtr<Ifpack_CrsRiluk> ilukFactors;
  //
  if (Lfill > -1) {
    ilukGraph = Teuchos::rcp(new Ifpack_IlukGraph(A->Graph(), Lfill, Overlap));
    assert(ilukGraph->ConstructFilledGraph()==0);
    ilukFactors = Teuchos::rcp(new Ifpack_CrsRiluk(*ilukGraph));
    int initerr = ilukFactors->InitValues(*A);
    if (initerr != 0) cout << "InitValues error = " << initerr;
    assert(ilukFactors->Factor() == 0);
  }
  //
  bool transA = false;
  double Cond_Est;
  ilukFactors->Condest(transA, Cond_Est);
  if (proc_verbose) {
    cout << "Condition number estimate for this preconditoner = " << Cond_Est << endl;
    cout << endl;
  }
  //
  // Solve using Belos
  //
  // Construct a Belos::Operator instance through the Epetra interface.
  //
  Belos::EpetraOp Amat( A );
  //
  // call the ctor for the preconditioning object
  //
  Belos::EpetraPrecOp Prec( ilukFactors );
  //
  // ********Other information used by block solver***********
  // *****************(can be user specified)******************
  //
  const int NumGlobalElements = Map.NumGlobalElements();
  int maxits = NumGlobalElements - 1; // maximum number of iterations to run
  //
  ParameterList My_PL;
  My_PL.set( "Length", length );
  //
  // *****Construct solution vector and random right-hand-sides *****
  //
  Belos::EpetraMultiVec soln(Map, numrhs);
  Belos::EpetraMultiVec rhs(Map, numrhs);
  rhs.MvRandom();
  Belos::LinearProblem<double,MV,OP>
  My_LP( rcp(&Amat,false), rcp(&soln,false), rcp(&rhs,false) );
  //My_LP.SetRightPrec( rcp(&Prec,false) );
  My_LP.SetLeftPrec( rcp(&Prec,false) );
  My_LP.SetBlockSize( blocksize );
  
  Belos::OutputManager<double> My_OM( MyPID );
  if (verbose)
    My_OM.SetVerbosity( Belos::Errors + Belos::Warnings 
			+ Belos::TimingDetails + Belos::FinalSummary );

  typedef Belos::StatusTestCombo<double,MV,OP>  StatusTestCombo_t;
  typedef Belos::StatusTestResNorm<double,MV,OP>  StatusTestResNorm_t;
  Belos::StatusTestMaxIters<double,MV,OP> test1( maxits );
  Belos::StatusTestMaxRestarts<double,MV,OP> test2( numrestarts );
  StatusTestCombo_t BasicTest( StatusTestCombo_t::OR, test1, test2 );
  StatusTestResNorm_t test3( tol );
  test3.DefineScaleForm( StatusTestResNorm_t::NormOfPrecInitRes, Belos::TwoNorm );
  Belos::StatusTestOutputter<ST,MV,OP> test4( frequency, false, "Native Residual: ||A*x-b||/||b||" );
  test4.set_resNormStatusTest( rcp(&test3, false) );
  test4.set_outputManager( rcp(&My_OM,false) );
  BasicTest.AddStatusTest( test4 );      
  StatusTestResNorm_t ExpTest( tol );
  ExpTest.DefineResForm( StatusTestResNorm_t::Explicit, Belos::TwoNorm ); 
  StatusTestCombo_t My_Test( StatusTestCombo_t::SEQ, BasicTest, ExpTest );  
  //
  // *******************************************************************
  // *************Start the block Gmres iteration*************************
  // *******************************************************************
  //
  Belos::PseudoBlockGmres<double,MV,OP>
    MyBlockGmres( rcp(&My_LP,false), rcp(&My_Test,false), rcp(&My_OM,false), rcp(&My_PL,false));
  //
  // **********Print out information about problem*******************
  //
  if (proc_verbose) {
    cout << endl << endl;
    cout << "Dimension of matrix: " << NumGlobalElements << endl;
    cout << "Number of right-hand sides: " << numrhs << endl;
    cout << "Number of linear systems solved together: " << blocksize << endl;
    cout << "Number of restarts allowed: " << numrestarts << endl;
    cout << "Length of block Arnoldi factorization: " << length << endl;
    cout << "Max number of Gmres iterations: " << maxits << endl; 
    cout << "Relative residual tolerance: " << tol << endl;
    cout << endl;
  }
  //
  //
  if (proc_verbose) {
    cout << endl << endl;
    cout << "Running Block Gmres -- please wait" << endl;
    cout << (numrhs+blocksize-1)/blocksize 
	 << " pass(es) through the solver required to solve for " << endl; 
    cout << numrhs << " right-hand side(s) -- using a block size of " << blocksize
	 << endl << endl;
  }

  //
  // Perform solve
  //
  MyBlockGmres.Solve();

  //
  // Compute actual residuals.
  //
  std::vector<double> actual_resids( numrhs );
  std::vector<double> rhs_norm( numrhs );
  Belos::EpetraMultiVec resid(Map, numrhs);
  OPT::Apply( Amat, soln, resid );
  MVT::MvAddMv( -1.0, resid, 1.0, rhs, resid ); 
  MVT::MvNorm( resid, &actual_resids );
  MVT::MvNorm( rhs, &rhs_norm );
  if (proc_verbose) {
    cout<< "---------- Actual Residuals (normalized) ----------"<<endl<<endl;
    for ( int i=0; i<numrhs; i++) {
      cout<<"Problem "<<i<<" : \t"<< actual_resids[i]/rhs_norm[i] <<endl;
    }
  }

  if (My_Test.GetStatus()!=Belos::Converged) {
	if (proc_verbose)
      		cout << "End Result: TEST FAILED" << endl;	
	return -1;
  }
  //
  // Default return value
  //
  if (proc_verbose)
    cout << "End Result: TEST PASSED" << endl;
  return 0;
  //
} // end test_bl_pgmres_hb.cpp
