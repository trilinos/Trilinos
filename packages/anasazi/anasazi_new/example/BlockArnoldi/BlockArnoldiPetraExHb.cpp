// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
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
//  This example compute the eigenvalues of a Harwell-Boeing matrix using the block Arnoldi
//  method.  The matrix is passed to the example routine through the command line, and 
//  converted to an Epetra matrix through some utilty routines.  This matrix is passed to the
//  eigensolver and then used to construct the Krylov decomposition.  The specifics of the 
//  block Arnoldi method can be set by the user.

#include "AnasaziPetraInterface.hpp"
#include "AnasaziBlockArnoldi.hpp"
#include "AnasaziBasicSort.hpp"
#include "AnasaziConfigDefs.hpp"
#include "Epetra_CrsMatrix.h"

#include "Trilinos_Util.h"

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
	 << " HB_filename " << endl
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
	A.SetTracebackMode(1); // Shutdown Epetra Warning tracebacks
	//
        //************************************
        // Start the block Arnoldi iteration
        //***********************************         
	//
        //  Variables used for the Block Arnoldi Method
        // 
        int block = 5;
        int length = 10;
        int nev = 5;
        double tol = 1.0e-8;
        string which="LM";
        int step = 5;
        int restarts = 10;
	//
        // create a PetraAnasaziVec. Note that the decision to make a view or
        // or copy is determined by the petra constructor called by Anasazi::PetraVec.
        // This is possible because I pass in arguements needed by petra.

        Anasazi::PetraVec ivec(Map, block);
        ivec.MvRandom();

        // call the ctor that calls the petra ctor for a matrix

        Anasazi::PetraOp Amat(A);
        Anasazi::Eigenproblem<double> MyProblem(&Amat, &ivec);

	// Inform the eigenproblem that the matrix A is symmetric
	//MyProblem.SetSymmetric(true);

	// Set the number of eigenvalues requested and the blocksize the solver should use
	MyProblem.SetNEV( nev );
	MyProblem.SetBlockSize( block );

        // Create a sorting manager to handle the sorting of eigenvalues in the solver
        Anasazi::BasicSort<double> MySort( which );

        // Create an output manager to handle the I/O from the solver
        Anasazi::OutputManager<double> MyOM( MyPID );
        //MyOM.SetVerbosity( 2 );
	//
	//  Initialize the Block Arnoldi solver
	//
        Anasazi::BlockArnoldi<double> MyBlockArnoldi(MyProblem, MySort, MyOM, tol, length,
                                         step, restarts);

#ifdef UNIX
        Epetra_Time & timer = *new Epetra_Time(Comm);
#endif

        // iterate a few steps (if you wish)
        //MyBlockArnoldi.iterate(5);

        // solve the problem to the specified tolerances or length
        MyBlockArnoldi.solve();
#ifdef UNIX
        double elapsed_time = timer.ElapsedTime();
        double total_flops = A.Flops(); 
        double MFLOPs = total_flops/elapsed_time/1000000.0;
#endif

        // obtain results directly
        double* resids = MyBlockArnoldi.getResiduals();
        double* evalr = MyBlockArnoldi.getEvals();
        double* evali = MyBlockArnoldi.getiEvals();

        // retrieve eigenvectors
        Anasazi::PetraVec evecr(Map, nev);
        MyBlockArnoldi.getEvecs( evecr );
        Anasazi::PetraVec eveci(Map, nev);
        MyBlockArnoldi.getiEvecs( eveci );

        // output results to screen
        MyBlockArnoldi.currentStatus();

#ifdef UNIX
        if (verbose)
                cout << "\n\nTotal MFLOPs for Arnoldi = " << MFLOPs << " Elapsed Time = "<<  elapsed_time <<endl;
#endif

        // Release all objects
        delete [] resids;
	delete [] evalr;
	delete [] evali;
        delete [] NumNz;

  	return 0;

} // end BlockArnoldiPetraExHb.cpp
