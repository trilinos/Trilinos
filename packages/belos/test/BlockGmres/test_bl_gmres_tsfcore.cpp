#include "TSFCoreEpetraMultiVector.hpp"
#include "TSFCoreEpetraLinearOp.hpp"
#include "TSFCoreEpetraVectorSpace.hpp"
#include "Trilinos_Util.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_MultiVector.h"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_Time.hpp"
#include "BelosTSFCoreInterface.hpp"
#include "BelosLinearProblemManager.hpp"
#include "BelosOutputManager.hpp"
#include "BelosStatusTestMaxIters.hpp"
#include "BelosStatusTestMaxRestarts.hpp"
#include "BelosStatusTestResNorm.hpp"
#include "BelosStatusTestCombo.hpp"
#include "BelosBlockGmres.hpp"
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
	double *xguess=0, *b=0, *xexact=0, *xsolve=0;
	Teuchos::Time timer("Gmres");

#ifdef EPETRA_MPI	
	// Initialize MPI	
	MPI_Init(&argc,&argv); 	
	Epetra_MpiComm Comm( MPI_COMM_WORLD );	
#else	
	Epetra_SerialComm Comm;	
#endif
	
        int MyPID = Comm.MyPID();

        bool verbose = (MyPID==0);
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
	Teuchos::RefCountPtr<Epetra_CrsMatrix> A =
		Teuchos::rcp( new Epetra_CrsMatrix(Copy, Map, NumNz) );
	// 
	// Add rows one-at-a-time
	//
	int NumEntries;
	for (i=0; i<NumMyElements; i++) {
		row_vals = val + bindx[i];
		col_inds = bindx + bindx[i];
		NumEntries = bindx[i+1] - bindx[i];
		assert(A->InsertGlobalValues(update[i], NumEntries, row_vals, col_inds)==0);
		assert(A->InsertGlobalValues(update[i], 1, val+i, update+i)==0);
	}
	//
	// Finish up
	//
	assert(A.get()->TransformToLocal()==0);
	assert(A.get()->OptimizeStorage()==0);
	//
	A->SetTracebackMode(1); // Shutdown Epetra Warning tracebacks
	//
	// Create the TSFCore/Belos Linear Operator
	//
	Teuchos::RefCountPtr<TSFCore::EpetraLinearOp> ELOp =
	  Teuchos::rcp( new TSFCore::EpetraLinearOp( A ) ); 
	Belos::TSFCoreMat<double> Amat( ELOp );
	//
        // ********Other information used by the GMRES block solver***********
        //*****************(can be user specified)******************
        //
        int numrhs = 1;  // total number of right-hand sides to solve for
        int block = 1;  // blocksize used by solver
        int numrestarts = 0; // number of restarts allowed
        int maxits = NumGlobalElements/block; // maximum number of iterations to run
        double tol = 1.0e-6;  // relative residual tolerance
	//
	// Create some Epetra_MultiVectors
	//
	Teuchos::RefCountPtr<Epetra_MultiVector> bb = 
	  Teuchos::rcp( new Epetra_MultiVector(Copy, Map, b, NumMyElements, numrhs) );
	Teuchos::RefCountPtr<Epetra_MultiVector> x = 
	  Teuchos::rcp( new Epetra_MultiVector(Map, numrhs) ); 
	Teuchos::RefCountPtr<Epetra_Map> rcp_map = Teuchos::rcp( &Map, false );
	Teuchos::RefCountPtr<TSFCore::EpetraVectorSpace> vs = 
	  Teuchos::rcp( new TSFCore::EpetraVectorSpace( rcp_map ) );
	// 
	// Create the TSFCore/Belos Multivectors
	//
	TSFCore::EpetraMultiVector RHS( bb, vs );
	TSFCore::EpetraMultiVector SOLN( x, vs );
	Belos::TSFCoreVec<double> rhs( RHS );
	Belos::TSFCoreVec<double> soln( SOLN );
	//
	// Set up linear problem instance.
	//
        Belos::LinearProblemManager<double> My_LP(&Amat, &soln, &rhs);
	My_LP.SetBlockSize( block );
	//
	//*******************************************************************
	// *************Start the GMRES iteration*************************
	//*******************************************************************
	//
        Belos::StatusTestMaxIters<double> test1( maxits );
        Belos::StatusTestMaxRestarts<double> test2( numrestarts );
        Belos::StatusTestCombo<double> test3( Belos::StatusTestCombo<double>::OR, test1, test2 );
        Belos::StatusTestResNorm<double> test4( tol );
        Belos::StatusTestCombo<double> My_Test( Belos::StatusTestCombo<double>::OR, test3, test4 );

        Belos::OutputManager<double> My_OM( MyPID );
	//My_OM.SetVerbosity( 3 );

        Belos::BlockGmres<double> MyBlockGmres(My_LP, My_Test, My_OM, maxits);

	timer.start();
	MyBlockGmres.Solve();
	timer.stop();
	
	My_Test.Print(cout);

	// Release all objects  	
	delete [] NumNz;
	delete [] bindx;
	delete [] val;
	if (xexact) delete [] xexact;
	if (xguess) delete [] xguess;
	if (b) delete [] b;

	return 0;
	//return success ? 0 : -1;

} // end TSFCoreSolversGmresTest.cpp
