//
// test_bl_cg_hb.cpp
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
// numrhs = 1 and block > 1 to sove a single rhs-system with a block implementation. 
//
// NOTE: No preconditioner is used in this case. 
//
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblemManager.hpp"
#include "BelosOutputManager.hpp"
#include "BelosStatusTestMaxIters.hpp"
#include "BelosStatusTestResNorm.hpp"
#include "BelosStatusTestCombo.hpp"
#include "BelosPetraInterface.hpp"
#include "BelosBlockCG.hpp"
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
	int i, j;
	int n_nonzeros, N_update;
	int *bindx=0, *update=0, *col_inds=0;
	double *val=0, *row_vals=0;
	double *xguess=0, *b=0, *xexact=0;
	Teuchos::Time timer("BelosCG");
	
#ifdef EPETRA_MPI	
	// Initialize MPI	
	MPI_Init(&argc,&argv); 	
	Epetra_MpiComm Comm( MPI_COMM_WORLD );	
#else	
	Epetra_SerialComm Comm;	
#endif
	
	int MyPID = Comm.MyPID();
	int NumProc = Comm.NumProc();
	bool verbose = (MyPID == 0);
	//
    	if(argc < 2 && verbose) {
     	cerr << "Usage: " << argv[0] 
	 << " HB_filename [level_fill [level_overlap [absolute_threshold [ relative_threshold]]]]" << endl
	 << "where:" << endl
	 << "HB_filename        - filename and path of a Harwell-Boeing data set" << endl
	 << endl;
    	return(1);
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
        Trilinos_Util_read_hb(argv[1], MyPID, &NumGlobalElements, &n_nonzeros,
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
	// Construct a Belos::Operator instance through the Petra interface.
	//
	Belos::PetraMat<double> Amat( &A );
	//
    	// ********Other information used by block solver***********
	//*****************(can be user specified)******************
	//
	int numrhs = 1;  // total number of right-hand sides to solve for
    	int block = 1;  // blocksize used by solver
	int maxits = NumGlobalElements/block - 1; // maximum number of iterations to run
    	double tol = 1.0e-6;  // relative residual tolerance
	//
	// Construct the right-hand side and solution multivectors.
	//
	Belos::PetraVec<double> rhs(Map, b, numrhs, NumMyElements);
	//Belos::PetraVec<double> rhs( Map, numrhs );
	//rhs.MvRandom();
	Belos::PetraVec<double> soln( Map, numrhs );
	Belos::PetraVec<double> xx(Map, xexact, numrhs, NumMyElements);
	//
	//  Construct an unpreconditioned linear problem instance.
	//
	Belos::LinearProblemManager<double> My_LP( &Amat, &soln, &rhs );
	My_LP.SetBlockSize( block );
	//
	//*******************************************************************
	// *************Start the block CG iteration*************************
	//*******************************************************************
	//
	Belos::StatusTestMaxIters<double> test1( maxits );
 	Belos::StatusTestResNorm<double> test2( tol );
	Belos::StatusTestCombo<double> My_Test( Belos::StatusTestCombo<double>::OR, test1, test2 );
	//
	Belos::OutputManager<double> My_OM( MyPID );
	//My_OM.SetVerbosity( 1 );
	//
	Belos::BlockCG<double> MyBlockCG(My_LP, My_Test, My_OM);
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
	//
    	//
	if (verbose) {
	   cout << endl << endl;
	   cout << "Running Block CG -- please wait" << endl;
	   cout << (numrhs+block-1)/block 
		    << " pass(es) through the solver required to solve for " << endl; 
	   cout << numrhs << " right-hand side(s) -- using a block size of " << block
			<< endl << endl;
	}
	timer.start();
	MyBlockCG.Solve();
	timer.stop();
	My_Test.Print(cout);

	if (verbose) {
		cout << "Solution time : "<< timer.totalElapsedTime()<<endl;
	}
	
// Release all objects  

  delete [] NumNz;
  delete [] val;
  delete [] xexact; 
  delete [] b; 
  delete [] xguess;
  delete [] update; 
  delete [] bindx;
  return 0;
  //
} // end test_bl_pgmrs_hb.cpp




