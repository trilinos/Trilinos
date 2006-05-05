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
// The right-hand-side from the HB file is used instead of random vectors.
// The initial guesses are all set to zero. 
//
// NOTE: No preconditioner is used in this case. 
//
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosOutputManager.hpp"
#include "BelosStatusTestMaxIters.hpp"
#include "BelosStatusTestOutputter.hpp"
#include "BelosStatusTestResNorm.hpp"
#include "BelosStatusTestCombo.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosBlockCG.hpp"
#include "createEpetraProblem.hpp"
#include "Trilinos_Util.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Teuchos_CommandLineProcessor.hpp"

int main(int argc, char *argv[]) {
  //
#ifdef EPETRA_MPI	
  // Initialize MPI	
  MPI_Init(&argc,&argv); 	
  Belos::MPIFinalize mpiFinalize; // Will call finalize with *any* return
#endif
  //
  using Teuchos::RefCountPtr;
  using Teuchos::rcp;
  //
  // Get test parameters from command-line processor
  //  
  bool verbose = false;
  int frequency = -1;  // how often residuals are printed by solver
  int numrhs = 1;  // total number of right-hand sides to solve for
  int blockSize = 1;  // blocksize used by solver
  std::string filename("bcsstk14.hb");
  double tol = 1.0e-5;  // relative residual tolerance

  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("frequency",&frequency,"Solvers frequency for printing residuals (#iters).");
  cmdp.setOption("tol",&tol,"Relative residual tolerance used by CG solver.");
  cmdp.setOption("filename",&filename,"Filename for Harwell-Boeing test matrix.");
  cmdp.setOption("num-rhs",&numrhs,"Number of right-hand sides to be solved for.");
  cmdp.setOption("block-size",&blockSize,"Block size to be used by CG solver.");
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
  RefCountPtr<Epetra_MultiVector> B, X;
  int return_val =Belos::createEpetraProblem(filename,NULL,&A,&B,&X,&MyPID);
  if(return_val != 0) return return_val;
  const Epetra_Map &Map = A->RowMap();
  verbose &= (MyPID==0); /* Only print on the zero processor */
  //
  // Solve using Belos
  //
  typedef double                          ST;
  typedef Epetra_Operator                 OP;
  typedef Epetra_MultiVector              MV;
  typedef Belos::OperatorTraits<ST,MV,OP> OPT;
  typedef Belos::MultiVecTraits<ST,MV>    MVT;
  //
  // ********Other information used by block solver***********
  // *****************(can be user specified)******************
  //
  const int NumGlobalElements = B->GlobalLength();
  int maxits = NumGlobalElements/blockSize - 1; // maximum number of iterations to run
  //
  // *****Construct initial guess and random right-hand-sides *****
  //
  if (numrhs != 1) {
    X = rcp( new Epetra_MultiVector( Map, numrhs ) );
    MVT::MvRandom( *X );
    B = rcp( new Epetra_MultiVector( Map, numrhs ) );
    OPT::Apply( *A, *X, *B );
    MVT::MvInit( *X, 0.0 );
  }
  //
  //  Construct an unpreconditioned linear problem instance.
  //
  Belos::LinearProblem<ST,MV,OP> My_LP( A, X, B );
  My_LP.SetBlockSize( blockSize );
  //
  // *******************************************************************
  // *************Start the block CG iteration*************************
  // *******************************************************************
  //
  Belos::OutputManager<ST> My_OM( MyPID );
  if (verbose)
    My_OM.SetVerbosity( Belos::Errors + Belos::Warnings 
			+ Belos::TimingDetails + Belos::FinalSummary );
  //
  Belos::StatusTestMaxIters<ST,MV,OP> test1( maxits );
  Belos::StatusTestResNorm<ST,MV,OP>  test2( tol );
  Belos::StatusTestOutputter<ST,MV,OP> test3( frequency, false );
  test3.set_resNormStatusTest( rcp(&test2,false) );
  test3.set_outputManager( rcp(&My_OM,false) );
  Belos::StatusTestCombo<ST,MV,OP> My_Test( Belos::StatusTestCombo<ST,MV,OP>::OR, test1, test3 );
  //
  Belos::BlockCG<ST,MV,OP>
    MyBlockCG(rcp(&My_LP,false), rcp(&My_Test,false), rcp(&My_OM,false) );
  //
  // **********Print out information about problem*******************
  //
  if (verbose) {
    cout << endl << endl;
    cout << "Dimension of matrix: " << NumGlobalElements << endl;
    cout << "Number of right-hand sides: " << numrhs << endl;
    cout << "Block size used by solver: " << blockSize << endl;
    cout << "Max number of CG iterations: " << maxits << endl; 
    cout << "Relative residual tolerance: " << tol << endl;
    cout << endl;
  }
  //
  //
  if (verbose) {
    cout << endl << endl;
    cout << "Running Block CG -- please wait" << endl;
    cout << (numrhs+blockSize-1)/blockSize 
	 << " pass(es) through the solver required to solve for " << endl; 
    cout << numrhs << " right-hand side(s) -- using a block size of " << blockSize
	 << endl << endl;
  }

  MyBlockCG.Solve();
  
  if (My_Test.GetStatus()!=Belos::Converged) {
    if (verbose)
      cout << "End Result: TEST FAILED" << endl;	
    return -1;
  }
  //
  // Default return value
  //
  if (verbose)
    cout << "End Result: TEST PASSED" << endl;
  return 0;
  //
} // end test_bl_cg_hb.cpp




