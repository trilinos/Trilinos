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
#include "Trilinos_Util.h"
#include "Epetra_CrsMatrix.h"
#include "Teuchos_Time.hpp"
//
//
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"

int main(int argc, char *argv[]) {
	//
	int i;
	int n_nonzeros, N_update;
	int *bindx=0, *update=0, *col_inds=0;
	double *val=0, *row_vals=0;
	double *xguess=0, *b=0, *xexact=0;
	Teuchos::Time timer("Belos");	

#ifdef EPETRA_MPI	
	// Initialize MPI	
	MPI_Init(&argc,&argv); 	
	Epetra_MpiComm Comm( MPI_COMM_WORLD );	
#else	
	Epetra_SerialComm Comm;	
#endif
	
	int MyPID = Comm.MyPID();
	
	bool verbose = 0;
	//
        if((argc < 2 || argc > 4)&& MyPID==0) {
        cerr << "Usage: " << argv[0]
         << " [ -v ] [ HB_filename ]" << endl
         << "where:" << endl
         << "-v                 - run test in verbose mode" << endl
         << "HB_filename        - filename and path of a Harwell-Boeing data set" << endl
         << endl;
        return(1);
        }
        //
        // Find verbosity flag
        //
        int file_arg = 1;
        for(i = 1; i < argc; i++)
        {
          if(argv[i][0] == '-' && argv[i][1] == 'v') {
            verbose = (MyPID == 0);
            if(i==1) file_arg = 2;
          }
        }
	//
	//**********************************************************************
	//******************Set up the problem to be solved*********************
	//**********************************************************************
    	//
    	int NumGlobalElements;  // total # of rows in matrix
	//
	// *****Read in matrix from HB file******
	//
        Trilinos_Util_read_hb(argv[file_arg], MyPID, &NumGlobalElements, &n_nonzeros,
                              &val, &bindx, &xguess, &b, &xexact);
        // 
        // *****Distribute data among processors*****
        //
        Trilinos_Util_distrib_msr_matrix(Comm, &NumGlobalElements, &n_nonzeros, &N_update, 
                                         &update, &val, &bindx, &xguess, &b, &xexact);
	//
	// *****Construct the matrix*****
	//
	int NumMyElements = N_update; // # local rows of matrix on processor
	//
    	// Create an integer vector NumNz that is used to build the Petra Matrix.
	// NumNz[i] is the Number of OFF-DIAGONAL term for the ith global equation 
	// on this processor
	//
	int * NumNz = new int[NumMyElements];
	for (i=0; i<NumMyElements; i++) {
		NumNz[i] = bindx[i+1] - bindx[i] + 1;
	}
	//
	Epetra_Map Map(NumGlobalElements, NumMyElements, update, 0, Comm);
	//
	// Create a Epetra_Matrix
	//
	Epetra_CrsMatrix A(Copy, Map, NumNz);
	//
	// Add rows one-at-a-time
	//
	int NumEntries;
	for (i=0; i<NumMyElements; i++) {
		row_vals = val + bindx[i];
		col_inds = bindx + bindx[i];
		NumEntries = bindx[i+1] - bindx[i];
		assert(A.InsertGlobalValues(update[i], NumEntries, row_vals, col_inds)==0);
		assert(A.InsertGlobalValues(update[i], 1, val+i, update+i)==0);
	}
	//
	// Finish up
	//
	assert(A.TransformToLocal()==0);
	assert(A.OptimizeStorage()==0);
	A.SetTracebackMode(1); // Shutdown Epetra Warning tracebacks
	//
	// Construct a Belos::Operator instance through the Epetra interface.
	//
	Belos::PetraMat<double> Amat( &A );
	//
    	// ********Other information used by block solver***********
	//*****************(can be user specified)******************
	//
	int numrhs = 1;  // total number of right-hand sides to solve for
    	int block = 1;  // blocksize used by solver
	int numrestarts = 3; // number of restarts allowed 
	int maxits = NumGlobalElements/block - 1; // maximum number of iterations to run
    	double tol = 1.0e-6;  // relative residual tolerance
	//
	// Construct the right-hand side and solution multivectors.
	//
	Belos::PetraVec<double> rhs(Map, b, numrhs, NumMyElements);
	Belos::PetraVec<double> soln( Map, numrhs );
	Belos::PetraVec<double> xx(Map, xexact, numrhs, NumMyElements);
	//
	// Construct an unpreconditioned linear problem instance.
	//
	Belos::LinearProblemManager<double> My_LP( rcp(&Amat, false), rcp(&soln,false), rcp(&rhs,false));
	My_LP.SetBlockSize( block );
	//
	//*******************************************************************
	// *************Start the block Gmres iteration*************************
	//*******************************************************************
	//
	Belos::StatusTestMaxIters<double> test1( maxits );
	Belos::StatusTestMaxRestarts<double> test2( numrestarts );
	Belos::StatusTestCombo<double> test3( Belos::StatusTestCombo<double>::OR, test1, test2 );
	Belos::StatusTestResNorm<double> test4( tol );
	Belos::StatusTestCombo<double> My_Test( Belos::StatusTestCombo<double>::OR, test3, test4 );

	Belos::OutputManager<double> My_OM( MyPID );
	if (verbose)
	  My_OM.SetVerbosity( 2 );

	Belos::BlockGmres<double> MyBlockGmres( rcp(&My_LP,false), rcp(&My_Test,false), rcp(&My_OM,false), maxits);

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
	
	timer.start();
	MyBlockGmres.Solve();
	timer.stop();

	if (verbose) {
		cout << "Solution time: "<<timer.totalElapsedTime()<<endl;
	}
	
// Release all objects  

  delete [] NumNz;
  delete [] update;
  delete [] val;
  delete [] bindx;
  if (xexact) delete [] xexact;
  if (xguess) delete [] xguess;
  if (b) delete [] b;
	
  if (My_Test.GetStatus() == Belos::Converged) {
    cout<< "***************** The test PASSED !!!********************"<<endl;
    return 0;
  }
  cout<< "********************The test FAILED!!! ********************"<<endl;
  return 1;
  //
} // end test_bl_gmres_hb.cpp
