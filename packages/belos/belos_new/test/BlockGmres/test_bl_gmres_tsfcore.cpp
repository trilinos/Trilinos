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
// The right-hand-side from the problem is being used instead of multiple
// random right-hand-sides.  The initial guesses are all set to zero. 
//
// NOTE: No preconditioner is used in this case. 
//
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblemManager.hpp"
#include "BelosOutputManager.hpp"
#include "BelosStatusTestMaxIters.hpp"
#include "BelosStatusTestMaxRestarts.hpp"
#include "BelosStatusTestResNorm.hpp"
#include "BelosStatusTestCombo.hpp"
#include "BelosPetraInterface.hpp"
#include "BelosBlockGmres.hpp"
#include "createEpetraProblem.hpp"
#include "TSFCoreEpetraLinearOp.hpp"
#include "TSFCoreEpetraMultiVector.hpp"
#include "TSFCoreEpetraVectorSpace.hpp"
#include "Epetra_CrsMatrix.h"
#include "Teuchos_Time.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_Map.h"

int main(int argc, char *argv[]) {
#ifdef EPETRA_MPI	
  MPI_Init(&argc,&argv);
  Belos::MPIFinalize mpiFinalize; // Will call finalize with *any* return
#endif
  //
  using Teuchos::ParameterList;
  using Teuchos::RefCountPtr;
  using Teuchos::rcp;
  Teuchos::Time timer("Belos Gmres");
  bool success = true;
  bool verbose = true;
  try {
    //
    // Get the problem
    //
    RefCountPtr<Epetra_Map> rowMap;
    RefCountPtr<Epetra_CrsMatrix> A;
    RefCountPtr<Epetra_MultiVector> B, X;
    int MyPID;
    int return_val =Belos::createEpetraProblem(argc,argv,&rowMap,&A,&B,&X,&MyPID,&verbose);
    if(return_val != 0) return return_val;
    //
    // Solve using Belos
    //
    typedef TSFCore::LinearOp<double> OP;
    typedef TSFCore::MultiVector<double> MV;
    //
    // Construct a Belos::Operator instance through the Epetra interface.
    //
    TSFCore::EpetraLinearOp Amat(A);
    RefCountPtr<TSFCore::EpetraVectorSpace> vs = rcp(new TSFCore::EpetraVectorSpace(rowMap));
    //
    // ********Other information used by block solver***********
    // *****************(can be user specified)******************
    //
    const int NumGlobalElements = B->GlobalLength();
    int numrhs = 1;  // total number of right-hand sides to solve for
    int block = 1;  // blocksize used by solver
    int numrestarts = 3; // number of restarts allowed 
    int maxits = NumGlobalElements/block - 1; // maximum number of iterations to run
    double tol = 1.0e-6;  // relative residual tolerance
    //
    ParameterList My_PL;
    My_PL.set( "Length", maxits );    // Maximum number of blocks in Krylov factorization
    //
    // Construct the right-hand side and solution multivectors.
    //
    TSFCore::EpetraMultiVector rhs(B,vs);
    TSFCore::EpetraMultiVector soln(X,vs);
    //
    // Construct an unpreconditioned linear problem instance.
    //
    Belos::LinearProblemManager<double,OP,MV>
      My_LP( rcp(&Amat, false), rcp(&soln,false), rcp(&rhs,false));
    My_LP.SetBlockSize( block );
    //
    // *******************************************************************
    // *************Start the block Gmres iteration*************************
    // *******************************************************************
    //
    typedef Belos::StatusTestCombo<double,OP,MV>  StatusTestCombo_t;
    typedef Belos::StatusTestResNorm<double,OP,MV>  StatusTestResNorm_t;
    Belos::StatusTestMaxIters<double,OP,MV> test1( maxits );
    Belos::StatusTestMaxRestarts<double,OP,MV> test2( numrestarts );
    StatusTestCombo_t test3( StatusTestCombo_t::OR, test1, test2 );
    StatusTestResNorm_t test4( tol );
    StatusTestCombo_t My_Test( StatusTestCombo_t::OR, test3, test4 );
  
    Belos::OutputManager<double> My_OM( MyPID );
    if (verbose)
      My_OM.SetVerbosity( 2 );
  
    Belos::BlockGmres<double,OP,MV>
      MyBlockGmres( rcp(&My_LP,false), rcp(&My_Test,false), rcp(&My_OM,false), rcp(&My_PL,false) );
  
    //
    // **********Print out information about problem*******************
    //
    if (verbose) {
      cout << endl << endl;
      cout << "Dimension of matrix: " << NumGlobalElements << endl;
      cout << "Number of right-hand sides: " << numrhs << endl;
      cout << "Block size used by solver: " << block << endl;
      cout << "Number of restarts allowed: " << numrestarts << endl;
      cout << "Max number of Gmres iterations per restart cycle: " << maxits << endl; 
      cout << "Relative residual tolerance: " << tol << endl;
      cout << endl;
    }
    //
    //
    if (verbose) {
      cout << endl << endl;
      cout << "Running Block Gmres -- please wait" << endl;
      cout << (numrhs+block-1)/block 
	   << " pass(es) through the solver required to solve for " << endl; 
      cout << numrhs << " right-hand side(s) -- using a block size of " << block
	   << endl << endl;
    }
  
    timer.start(true);
    MyBlockGmres.Solve();
    timer.stop();

    if( My_Test.GetStatus() != Belos::Converged ) success = false;

    if(verbose)
      Teuchos::print_memory_usage_stats(Teuchos::get_default_workspace_store().get(),std::cout);

  }
  catch( const std::exception &excpt ) {
    if(verbose)
      std::cerr << "*** Caught a standard exception : " << excpt.what() << std::endl;
    success = false;
  }
  catch( ... ) {
    if(verbose)
      std::cerr << "*** Caught an unknown exception!\n";
    success = false;
  }
 
  if (verbose) {
    cout << "Solution time: "<<timer.totalElapsedTime()<<endl;
  }
  
  if ( success ) {
    cout<< "***************** The test PASSED !!!********************"<<endl;
    return 0;
  }
  cout<< "********************The test FAILED!!! ********************"<<endl;
  return 1;

} // end test_bl_gmres_tsfcore.cpp
