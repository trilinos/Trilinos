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
// numrhs = 1 and block > 1 to sove a single rhs-system with a block implementation. 
//
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
#include "Ifpack_IlukGraph.h"
#include "Ifpack_CrsRiluk.h"
#include "Epetra_CrsMatrix.h"
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
	 << "level_fill         - The amount of fill to use for ILU preconditioner (default 0)" << endl
	 << "level_overlap      - The amount of overlap used for overlapping Schwarz subdomains (default 0)" << endl
	 << "absolute_threshold - The minimum value to place on the diagonal prior to factorization (default 0.0)" << endl
	 << "relative_threshold - The relative amount to perturb the diagonal prior to factorization (default 1.0)" << endl << endl
	 << "To specify a non-default value for one of these parameters, you must specify all" << endl
	 << " preceding values but not any subsequent parameters. Example:" << endl
	 << "bl_pgmres_hb_mpi.exe mymatrix.hb 1  - loads mymatrix.hb, uses level fill of one, all other values are defaults" << endl
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
	Trilinos_Util_read_hb(argv[1], MyPID, &NumGlobalElements, &n_nonzeros, &val, 
		                    &bindx);
	//
	// *****Distribute data among processors*****
	//
	Trilinos_Util_distrib_msr_matrix(Comm, &NumGlobalElements, &n_nonzeros, &N_update,
		                             &update, &val, &bindx);
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
	//
	// call the ctor that calls the petra ctor for a matrix
	//
	Belos::PetraMat<double> Amat( &A );
	//
	A.SetTracebackMode(1); // Shutdown Epetra Warning tracebacks
	//
	//
	//*****Construct the Preconditioner*****
	//
	if (verbose) cout << endl << endl;
	if (verbose) cout << "Constructing ILU preconditioner" << endl;
	int Lfill = 0;
	if (argc > 2) Lfill = atoi(argv[2]);
	if (verbose) cout << "Using Lfill = " << Lfill << endl;
	int Overlap = 0;
	if (argc > 3) Overlap = atoi(argv[3]);
	if (verbose) cout << "Using Level Overlap = " << Overlap << endl;
	double Athresh = 0.0;
	if (argc > 4) Athresh = atof(argv[4]);
	if (verbose) cout << "Using Absolute Threshold Value of " << Athresh << endl;
	double Rthresh = 1.0;
	if (argc >5) Rthresh = atof(argv[5]);
	if (verbose) cout << "Using Relative Threshold Value of " << Rthresh << endl;
	//
	Ifpack_IlukGraph * ilukGraph=0;
	Ifpack_CrsRiluk * ilukFactors=0;
	//
	if (Lfill > -1) {
		ilukGraph = new Ifpack_IlukGraph(A.Graph(), Lfill, Overlap);
		assert(ilukGraph->ConstructFilledGraph()==0);
		ilukFactors = new Ifpack_CrsRiluk(*ilukGraph);
		int initerr = ilukFactors->InitValues(A);
		if (initerr != 0) cout << "InitValues error = " << initerr;
		assert(ilukFactors->Factor() == 0);
	}
	//
	bool transA = false;
	double Cond_Est;
	ilukFactors->Condest(transA, Cond_Est);
	if (verbose) {
		cout << "Condition number estimate for this preconditoner = " << Cond_Est << endl;
		cout << endl;
	}
	Epetra_Operator& prec = dynamic_cast<Epetra_Operator&>(*ilukFactors);
	//
	// call the ctor for the preconditioning object
	//
	Belos::PetraPrec<double> EpetraOpPrec( &prec );
	//
    	// ********Other information used by block solver***********
	//*****************(can be user specified)******************
	//
	int numrhs = 15;  // total number of right-hand sides to solve for
    	int block = 10;  // blocksize used by solver
	int numrestarts = 20; // number of restarts allowed 
    	int maxits = NumGlobalElements/block-1; // maximum number of iterations to run
	int length = 15;
    	double tol = 1.0e-6;  // relative residual tolerance
	//
	//*****Construct solution vector and random right-hand-sides *****
	//
	Belos::PetraVec<double> soln(Map, numrhs);
	Belos::PetraVec<double> rhs(Map, numrhs);
	rhs.MvRandom();
	Belos::LinearProblemManager<double> My_LP( &Amat, &soln, &rhs);
	My_LP.SetLeftPrec( &EpetraOpPrec );
	My_LP.SetBlockSize( block );

	Belos::StatusTestMaxIters<double> test1( maxits );
	Belos::StatusTestMaxRestarts<double> test2( numrestarts );
	Belos::StatusTestCombo<double> test3( Belos::StatusTestCombo<double>::OR, test1, test2 );
	Belos::StatusTestResNorm<double> test4( tol );
	Belos::StatusTestCombo<double> My_Test( Belos::StatusTestCombo<double>::OR, test3, test4 );

	Belos::OutputManager<double> My_OM( MyPID );
	//My_OM.SetVerbosity( 1 );
	//
	//*******************************************************************
	// *************Start the block Gmres iteration*************************
	//*******************************************************************
	//
	Belos::BlockGmres<double> MyBlockGmres(My_LP, My_Test, My_OM, length);
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
	MyBlockGmres.Solve();
	My_Test.Print(cout);

	// Release all objects  
	if (ilukGraph) delete ilukGraph;
	if (ilukFactors) delete ilukFactors;
	delete [] NumNz; 
	delete [] bindx;
	delete [] update;
	delete [] val;
	
	
  return 0;
  //
} // end test_bl_pgmrs_hb.cpp
