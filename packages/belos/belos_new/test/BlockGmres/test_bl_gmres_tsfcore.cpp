#include "TSFCoreEpetraMultiVector.hpp"
#include "TSFCoreEpetraLinearOp.hpp"
#include "TSFCoreMultiVectorStdOps.hpp"
#include "Trilinos_Util.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_MultiVector.h"
#include "Teuchos_RefCountPtrDecl.hpp"
#include "Teuchos_Time.hpp"
#include "BelosTSFCoreInterface.hpp"
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
	int NumProc = Comm.NumProc();
	// Set verbosity of output
	bool verbose = false;
	for( i=1; i<argc; i++ ) {
		if (argv[i][0]=='-' && argv[i][1]=='v' && MyPID==0 ) { verbose = true; };
	}
	if ( verbose ) { argc--; } // Decrement argument counter if one of the arguments is the verbosity flag.
	//
	if( argc < 4 ) {
		if ( verbose ) {
			cerr
				<< "Usage: " << argv[0] 
				<< " HB_filename max_iter tol [-v] " << endl
				<< "where:" << endl
				<< "HB_filename        - filename and path of a Harwell-Boeing data set" << endl
				<< "max_iter           - maximum number of iterations allowed in GMRES solve" <<endl
				<< "tol                - relative residual tolerance for GMRES solve" << endl
				<< "[-v]               - verbosity flag for debugging" << endl
				<< endl;
		}
		return(1);
	}

	if ( verbose ) {
		std::cout
			<< std::endl << std::endl
			<< "***\n"
			<< "*** Testing simple GMRES solver on test problem \'" << argv[1] << "\'\n"
			<< "***\n\n";
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
	TSFCore::EpetraLinearOp ELOp( A ); 
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
	Teuchos::RefCountPtr<Epetra_MultiVector> bb = Teuchos::rcp( new Epetra_MultiVector(Copy, Map, b, NumMyElements, numrhs) );
	Teuchos::RefCountPtr<Epetra_MultiVector> x = Teuchos::rcp( new Epetra_MultiVector(Map, numrhs) ); 
	// 
	// Create the TSFCore/Belos Multivectors
	//
	TSFCore::EpetraMultiVector RHS( bb );
	TSFCore::EpetraMultiVector Soln( x );
	Belos::TSFCoreVec<double> rhs( RHS );
	Belos::TSFCoreVec<double> iguess( Soln );
	//
	// Set solver characteristics
	//
	//MyBlockGmres.SetInitGuess( iguess );
        //MyBlockGmres.SetRestart(numrestarts);
        //MyBlockGmres.SetDebugLevel(0);
	//
	//*******************************************************************
	// *************Start the GMRES iteration*************************
	//*******************************************************************
	//
	//timer.start();
	//MySolver.solve( ELOp, RHS, &Soln, TSFCore::NOTRANS, max_iter, tol );
	//timer.stop();
	//
	// Compute actual residual norm.
	//
	//double bnorm = norm_1( RHS );
	//ELOp.apply( TSFCore::NOTRANS, Soln, &RHS, 1.0, -1.0 );
	//double final_rel_err = norm_1( RHS ) / bnorm;
	
	if (verbose) {
	  //cout << "******************* Results ************************"<<endl;
	  //cout << "Iteration "<< MySolver.currIteration()<<" of "<<max_iter<< endl;
	  //cout << "Final Computed GMRES Relative Residual Norm : " << 
	  //  MySolver.currEstRelResidualNorm() << endl;
	  //cout << "Actual Computed GMRES Relative Residual Norm : " <<
	  //  final_rel_err << endl;
	  //cout << "Solution time:  "<< timer.totalElapsedTime() <<endl;
	  //cout << "****************************************************"<<endl;
	}
      	
	// Release all objects  	
	delete [] NumNz;
	delete [] bindx;

	//bool success = (final_rel_err <= tol);

	/*if (verbose) {
		if(success)
			std::cout << "\nCongradualtions! the system was solved to the specified tolerance!\n";
		else
			std::cout << "\nOh no! the system was not solved to the specified tolerance!\n";
	}*/
	
	return 0;
	//return success ? 0 : -1;

} // end TSFCoreSolversGmresTest.cpp
