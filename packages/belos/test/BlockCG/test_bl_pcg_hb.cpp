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
// As currently set up, this driver tests the case when the number of right-hand
// sides (numrhs = 15) is greater than the blocksize (block = 10) used by 
// the solver. Here, 2 passes through the solver are required to solve 
// for all right-hand sides. This information can be edited (see below - other
// information used by block solver - can be user specified) to solve for
// other sizes of systems. For example, one could set numrhs = 1 and block = 1,
// to solve a single right-hand side system in the traditional way, or, set
// numrhs = 1 and block > 1 to solve a single rhs-system with a block implementation. 
//
// 
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosOutputManager.hpp"
#include "BelosStatusTestMaxIters.hpp"
#include "BelosStatusTestResNorm.hpp"
#include "BelosStatusTestOutputter.hpp"
#include "BelosStatusTestCombo.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosBlockCG.hpp"
#include "createEpetraProblem.hpp"
#include "Trilinos_Util.h"
#include "Ifpack_CrsIct.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Teuchos_CommandLineProcessor.hpp"
//
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
  int frequency = -1; // how often residuals are printed by solver 
  int numrhs = 15;  // total number of right-hand sides to solve for
  int blockSize = 10;  // blocksize used by solver
  std::string filename("bcsstk14.hb");
  double tol = 1.0e-5;  // relative residual tolerance

  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("frequency",&frequency,"Solvers frequency for printing residuals (#iters).");
  cmdp.setOption("filename",&filename,"Filename for Harwell-Boeing test matrix.");
  cmdp.setOption("tol",&tol,"Relative residual tolerance used by CG solver.");
  cmdp.setOption("num-rhs",&numrhs,"Number of right-hand sides to be solved for.");
  cmdp.setOption("block-size",&blockSize,"Block size to be used by CG solver.");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }
  if (!verbose)
    frequency = -1;  // Reset frequency if verbosity is off
  //
  // Get the problem
  //
  int MyPID;
  RefCountPtr<Epetra_CrsMatrix> A;
  int return_val =Belos::createEpetraProblem(filename,NULL,&A,NULL,NULL,&MyPID);
  if(return_val != 0) return return_val;
  verbose &= (MyPID==0); /* Only print on the first processor */
  //
  // *****Select the Preconditioner*****
  //
  if (verbose) cout << endl << endl;
  if (verbose) cout << "Constructing ICT preconditioner" << endl;
  int Lfill = 0;
  // if (argc > 2) Lfill = atoi(argv[2]);
  if (verbose) cout << "Using Lfill = " << Lfill << endl;
  int Overlap = 0;
  // if (argc > 3) Overlap = atoi(argv[3]);
  if (verbose) cout << "Using Level Overlap = " << Overlap << endl;
  double Athresh = 0.0;
  // if (argc > 4) Athresh = atof(argv[4]);
  if (verbose) cout << "Using Absolute Threshold Value of " << Athresh << endl;
  double Rthresh = 1.0;
  // if (argc >5) Rthresh = atof(argv[5]);
  if (verbose) cout << "Using Relative Threshold Value of " << Rthresh << endl;
  double dropTol = 1.0e-6;
  //
  Teuchos::RefCountPtr<Ifpack_CrsIct> ICT;
  //
  if (Lfill > -1) {
    ICT = Teuchos::rcp( new Ifpack_CrsIct(*A, dropTol, Lfill) );
    ICT->SetAbsoluteThreshold(Athresh);
    ICT->SetRelativeThreshold(Rthresh);
    int initerr = ICT->InitValues(*A);
    if (initerr != 0) cout << "InitValues error = " << initerr;
    assert(ICT->Factor() == 0);
  }
  //
  bool transA = false;
  double Cond_Est;
  ICT->Condest(transA, Cond_Est);
  if (verbose) {
    cout << "Condition number estimate for this preconditoner = " << Cond_Est << endl;
    cout << endl;
  }
  //
  // Solve using Belos
  //
  typedef double                           ST;
  typedef Belos::Operator<ST>              OP;
  typedef Belos::MultiVec<ST>              MV;
  typedef Belos::OperatorTraits<ST,MV,OP> OPT;
  typedef Belos::MultiVecTraits<ST,MV>    MVT;
  //
  // Construct a Belos::Operator instance through the Epetra interface.
  //
  Belos::EpetraOp Amat( A );
  //
  // call the ctor for the preconditioning object
  //
  Belos::EpetraPrecOp Prec( ICT );
  //
  // ********Other information used by block solver***********
  // *****************(can be user specified)******************
  //
  const Epetra_Map &Map = A->RowMap();
  const int NumGlobalElements = Map.NumGlobalElements();
  int maxits = NumGlobalElements/blockSize - 1; // maximum number of iterations to run
  //
  // *****Construct initial guess and random right-hand-sides *****
  //
  Belos::EpetraMultiVec soln(Map, numrhs);
  Belos::EpetraMultiVec rhs(Map, numrhs);
  rhs.MvRandom();
  //
  // *****Create Linear Problem for Belos Solver
  //
  Belos::LinearProblem<ST,MV,OP> My_LP( rcp(&Amat, false), rcp(&soln, false), rcp(&rhs,false) );
  My_LP.SetLeftPrec( rcp(&Prec,false) );
  My_LP.SetBlockSize( blockSize );
  //
  // *****Create Status Test Class for the Belos Solver
  //
  Belos::OutputManager<ST> My_OM( MyPID );
  if (verbose)
    My_OM.SetVerbosity( Belos::Errors + Belos::Warnings 
			+ Belos::TimingDetails + Belos::FinalSummary );

  Belos::StatusTestMaxIters<ST,MV,OP> test1( maxits );
  Belos::StatusTestResNorm<ST,MV,OP> test2( tol );
  Belos::StatusTestOutputter<ST,MV,OP> test3( frequency, false );
  test3.set_resNormStatusTest( rcp(&test2,false) );
  test3.set_outputManager( rcp(&My_OM,false) );
  Belos::StatusTestCombo<ST,MV,OP> My_Test( Belos::StatusTestCombo<ST,MV,OP>::OR, test1, test3 );
  
  //
  // *******************************************************************
  // *************Start the block CG iteration*************************
  // *******************************************************************
  //
  Belos::BlockCG<ST,MV,OP> MyBlockCG( rcp(&My_LP, false), rcp(&My_Test,false), rcp(&My_OM,false));
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
  
  if (verbose) {
    cout << endl << endl;
    cout << "Running Block CG -- please wait" << endl;
    cout << (numrhs+blockSize-1)/blockSize 
	 << " pass(es) through the solver required to solve for " << endl; 
    cout << numrhs << " right-hand side(s) -- using a block size of " << blockSize
	 << endl << endl;
  }

  MyBlockCG.Solve();	

  //
  // Compute actual residuals.
  //
  std::vector<ST> actual_resids( numrhs );
  std::vector<ST> rhs_norm( numrhs );
  Belos::EpetraMultiVec resid( Map, numrhs );
  OPT::Apply( Amat, soln, resid );
  MVT::MvAddMv( -1.0, resid, 1.0, rhs, resid ); 
  MVT::MvNorm( resid, &actual_resids );
  MVT::MvNorm( rhs, &rhs_norm );
  if (verbose) {
    cout<< "---------- Actual Residuals (normalized) ----------"<<endl<<endl;
    for (int i=0; i<numrhs; i++) {
      cout<<"Problem "<<i<<" : \t"<< actual_resids[i]/rhs_norm[i] <<endl;
    }
  }

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
} // end test_bl_pcg_hb.cpp

