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
#include "BelosLinearProblemManager.hpp"
#include "BelosOutputManager.hpp"
#include "BelosStatusTestMaxIters.hpp"
#include "BelosStatusTestResNorm.hpp"
#include "BelosStatusTestCombo.hpp"
#include "BelosPetraInterface.hpp"
#include "BelosBlockCG.hpp"
#include "createEpetraProblem.hpp"
#include "TSFCoreEpetraLinearOp.hpp"
#include "TSFCoreEpetraMultiVector.hpp"
#include "TSFCoreEpetraVectorSpace.hpp"
#include "Trilinos_Util.h"
#include "Ifpack_CrsIct.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
//
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
  Teuchos::Time timer("Belos Preconditioned CG");
  bool success = true;
  bool verbose = true;
  try {
    //
    // Get the problem
    //
    RefCountPtr<Epetra_Map> rowMap;
    RefCountPtr<Epetra_CrsMatrix> A;
    int MyPID;
    bool verbose;
    int return_val =Belos::createEpetraProblem(argc,argv,&rowMap,&A,NULL,NULL,&MyPID,&verbose);
    if(return_val != 0) return return_val;
    //
    // *****Construct the Preconditioner*****
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
    Ifpack_CrsIct* ICT = 0;
    //
    if (Lfill > -1) {
      ICT = new Ifpack_CrsIct(*A, dropTol, Lfill);
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
    Epetra_Operator& prec = *ICT;
    //
    // Solve using Belos
    //
    typedef TSFCore::LinearOp<double> OP;
    typedef TSFCore::MultiVector<double> MV;
    typedef Belos::OperatorTraits<double, MV, OP> OPT;
    typedef Belos::MultiVecTraits<double, MV> MVT;
    //
    // Construct a Belos::Operator instance through the Epetra interface.
    //
    TSFCore::EpetraLinearOp Amat(A);
    RefCountPtr<TSFCore::EpetraVectorSpace> vs = rcp(new TSFCore::EpetraVectorSpace(rowMap));
    //
    // call the ctor for the preconditioning object
    //
    TSFCore::EpetraLinearOp Prec( rcp(&prec,false), TSFCore::NOTRANS, TSFCore::EPETRA_OP_APPLY_APPLY_INVERSE );
    //
    // ********Other information used by block solver***********
    //*****************(can be user specified)******************
    //
    const Epetra_Map &Map = A->RowMap();
    const int NumGlobalElements = Map.NumGlobalElements();
    //int numrhs = 15;  // total number of right-hand sides to solve for
    int numrhs = 10;
    int block = 10;  // blocksize used by solver
    //int block = 5;
    //int block = 15;
    int maxits = NumGlobalElements/block - 1; // maximum number of iterations to run
    double tol = 1.0e-6;  // relative residual tolerance
    //
    //*****Construct solution vector and random right-hand-sides *****
    //
    RefCountPtr<Epetra_MultiVector> B = rcp( new Epetra_MultiVector(*rowMap, numrhs) );
    B->SetSeed(0);
    B->Random();
    RefCountPtr<Epetra_MultiVector> X = rcp( new Epetra_MultiVector(*rowMap, numrhs) );
    TSFCore::EpetraMultiVector rhs( B, vs );
    TSFCore::EpetraMultiVector soln( X, vs );
    //
    //*****Create Linear Problem for Belos Solver
    //
    Belos::LinearProblemManager<double,OP,MV> My_LP( rcp(&Amat, false), rcp(&soln, false), rcp(&rhs,false) );
    My_LP.SetLeftPrec( rcp(&Prec, false) );
    My_LP.SetBlockSize( block );
    //
    //*****Create Status Test Class for the Belos Solver
    //
    Belos::StatusTestMaxIters<double,OP,MV> test1( maxits );
    Belos::StatusTestResNorm<double,OP,MV> test2( tol );
    Belos::StatusTestCombo<double,OP,MV> My_Test( Belos::StatusTestCombo<double,OP,MV>::OR, test1, test2 );

    Belos::OutputManager<double> My_OM( MyPID );
    if (verbose)
      My_OM.SetVerbosity( 2 );
    //
    // *******************************************************************
    // *************Start the block CG iteration*************************
    // *******************************************************************
    //
    Belos::BlockCG<double,OP,MV> MyBlockCG( rcp(&My_LP, false), rcp(&My_Test,false), rcp(&My_OM,false));
    //
    // **********Print out information about problem*******************
    //
    if (verbose) {
      cout << endl << endl;
      cout << "Dimension of matrix: " << NumGlobalElements << endl;
      cout << "Number of right-hand sides: " << numrhs << endl;
      cout << "Block size used by solver: " << block << endl;
      cout << "Max number of CG iterations: " << maxits << endl; 
      cout << "Relative residual tolerance: " << tol << endl;
      cout << endl;
    }

    if (verbose) {
      cout << endl << endl;
      cout << "Running Block CG -- please wait" << endl;
      cout << (numrhs+block-1)/block 
           << " pass(es) through the solver required to solve for " << endl; 
      cout << numrhs << " right-hand side(s) -- using a block size of " << block
           << endl << endl;
    }
    timer.start(true);
    MyBlockCG.Solve();	
    timer.stop();
    //
    if (My_Test.GetStatus() != Belos::Converged ) success = false; 
    //
    if(verbose)
      Teuchos::print_memory_usage_stats(Teuchos::get_default_workspace_store().get(),std::cout);
    //
    // Compute actual residuals.
    //
    double* actual_resids = new double[numrhs];
    double* rhs_norm = new double[numrhs];
    RefCountPtr<Epetra_MultiVector> R = rcp( new Epetra_MultiVector(*rowMap, numrhs) );
    TSFCore::EpetraMultiVector resid( R, vs );
    OPT::Apply( Amat, soln, resid );
    MVT::MvAddMv( -1.0, resid, 1.0, rhs, resid ); 
    MVT::MvNorm( resid, actual_resids );
    MVT::MvNorm( rhs, rhs_norm );
    if (verbose) {
      cout<< "---------- Actual Residuals (normalized) ----------"<<endl<<endl;
      for (int i=0; i<numrhs; i++) {
         cout<<"Problem "<<i<<" : \t"<< actual_resids[i]/rhs_norm[i] <<endl;
      }
    }
    // Release all objects  

    if (ICT) { delete ICT; ICT = 0; }
    delete [] actual_resids;
    delete [] rhs_norm;
  }
  catch ( const std::exception &excpt ) {
    if (verbose)
      std::cerr << "*** Caught a standard exception : "<< excpt.what() << std::endl;
    success = false;
  }
  catch ( ... ) {
    if (verbose)
      std::cerr << "*** Caught an unknown exception!\n";
    success = false;
  }
	
  if (verbose)
    cout << "Solution time: "<< timer.totalElapsedTime()<<endl;
  	
  if ( success ) {
    cout<< "***************** The test PASSED !!!********************"<<endl;
    return 0;
  }
  cout<< "********************The test FAILED!!! ********************"<<endl;
  return 1; 
  //
} // end test_bl_pcg_hb.cpp

