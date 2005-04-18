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

#include "AnasaziBlockKrylovSchur.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziBasicSort.hpp"
#include "AnasaziConfigDefs.hpp"
#include "AnasaziEpetraAdapter.hpp"

#include "Epetra_CrsMatrix.h"
#include "Teuchos_LAPACK.hpp"

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
	double zero = 0.0;

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
    	// Create an integer vector NumNz that is used to build the Epetra Matrix.
	// NumNz[i] is the Number of OFF-DIAGONAL term for the ith global equation 
	// on this processor
	//
	std::vector<int> NumNz(NumMyElements);
	for (i=0; i<NumMyElements; i++) {
		NumNz[i] = bindx[i+1] - bindx[i] + 1;
	}
	//
	Epetra_Map Map(NumGlobalElements, NumMyElements, update, 0, Comm);
	//
	// Create a Epetra_Matrix
	//
	Teuchos::RefCountPtr<Epetra_CrsMatrix> A = Teuchos::rcp( new Epetra_CrsMatrix(Copy, Map, &NumNz[0]) );
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
	assert(A->FillComplete()==0);
	A->SetTracebackMode(1); // Shutdown Epetra Warning tracebacks
	//
        //************************************
        // Start the block Arnoldi iteration
        //***********************************         
	//
        //  Variables used for the Block Arnoldi Method
        // 
        int nev = 5;
        int blockSize = 5;
        int maxBlocks = 10;
        int maxRestarts = 10;
        int step = 5;
        double tol = 1.0e-8;
        string which="LM";
	//
	// Create parameter list to pass into solver
	//
	Teuchos::ParameterList MyPL;
	MyPL.set( "Block Size", blockSize );
	MyPL.set( "Max Blocks", maxBlocks );
	MyPL.set( "Max Restarts", maxRestarts );
	MyPL.set( "Tol", tol );
	MyPL.set( "Step Size", step );

	typedef Epetra_MultiVector MV;
	typedef Epetra_Operator OP;
	typedef Anasazi::MultiVecTraits<double, MV> MVT;
	//
        // Create the eigenproblem to be solved.
        //
        Teuchos::RefCountPtr<Epetra_MultiVector> ivec = Teuchos::rcp( new Epetra_MultiVector(Map, blockSize) );
        ivec->Random();
	
	Teuchos::RefCountPtr<Anasazi::BasicEigenproblem<double, MV, OP> > MyProblem =
	  Teuchos::rcp( new Anasazi::BasicEigenproblem<double, MV, OP>(A, ivec) );
	
	// Inform the eigenproblem that the matrix A is symmetric
	//MyProblem->SetSymmetric(true);

	// Set the number of eigenvalues requested 
	MyProblem->SetNEV( nev );

       // Inform the eigenproblem that you are finishing passing it information
        assert( MyProblem->SetProblem() == 0 );

        // Create a sorting manager to handle the sorting of eigenvalues in the solver
	Teuchos::RefCountPtr<Anasazi::BasicSort<double, MV, OP> > MySort = 
	  Teuchos::rcp( new Anasazi::BasicSort<double, MV, OP>(which) );
	
        // Create an output manager to handle the I/O from the solver
	Teuchos::RefCountPtr<Anasazi::OutputManager<double> > MyOM =
	  Teuchos::rcp( new Anasazi::OutputManager<double>( MyPID ) );
        MyOM->SetVerbosity( Anasazi::FinalSummary );
	//
	//  Initialize the Block Arnoldi solver
	//
        Anasazi::BlockKrylovSchur<double, MV, OP> MySolver(MyProblem, MySort, MyOM, MyPL);
	
        // solve the problem to the specified tolerances or length
        MySolver.solve();
	
        // obtain results directly
        Teuchos::RefCountPtr<std::vector<double> > evals = MyProblem->GetEvals();
	
	// retrieve eigenvectors
	// The size of the eigenvector storage is nev + block, but the eigenvectors are stored in the first nev vectors.
        Teuchos::RefCountPtr<Epetra_MultiVector> evecs = MyProblem->GetEvecs();
        Teuchos::RefCountPtr<Epetra_MultiVector> evecr, eveci;
	if (MyProblem->IsSymmetric())
          evecr = evecs;
        else {
          evecr = Teuchos::rcp( new Epetra_MultiVector( View, *evecs, 0, nev ) );
          eveci = Teuchos::rcp( new Epetra_MultiVector( View, *evecs, nev, nev ) );
        }                 	

	// Compute residuals.
	Teuchos::LAPACK<int,double> lapack;
	Epetra_MultiVector tempevecr(Map,nev), tempAevec(Map,nev);
	Epetra_MultiVector tempeveci(Map,nev);
	Teuchos::SerialDenseMatrix<int,double> Breal(nev,nev), Breal2(nev,nev);
	Teuchos::SerialDenseMatrix<int,double> Bimag(nev,nev), Bimag2(nev,nev);
	std::vector<double> normA(nev), tempnrm(nev);
	cout<<endl<< "Actual Residuals"<<endl;
	cout<<"------------------------------------------------------"<<endl;
	Breal.putScalar(0.0); 
	if (!MyProblem->IsSymmetric())
	  Bimag.putScalar(0.0);
	for (i=0; i<nev; i++) { 
	  normA[i] = 0.0;
	  Breal(i,i) = (*evals)[i]; 
	  if (!MyProblem->IsSymmetric())
	    Bimag(i,i) = (*evals)[nev + i]; 
	}
	A->Apply( *evecr, tempAevec );
	MVT::MvTimesMatAddMv( -1.0, *evecr , Breal, 1.0, tempAevec );
	if (!MyProblem->IsSymmetric()) {
	  MVT::MvTimesMatAddMv( 1.0, *eveci, Bimag, 1.0, tempAevec );
	  MVT::MvNorm( tempAevec, &normA );
	  A->Apply( *eveci, tempAevec );
	  MVT::MvTimesMatAddMv( -1.0, *evecr, Bimag, 1.0, tempAevec );
	  MVT::MvTimesMatAddMv( -1.0, *eveci, Breal, 1.0, tempAevec );
	}
	MVT::MvNorm( tempAevec, &tempnrm );
	i = 0;
	while (i < nev) {
	  normA[i] = lapack.LAPY2( normA[i], tempnrm[i] );
	  if (MyProblem->IsSymmetric()) {
	    normA[i] /= Teuchos::ScalarTraits<double>::magnitude((*evals)[i]);
	    i++;
	  } else {
	    normA[i] /= lapack.LAPY2( (*evals)[i], (*evals)[nev + i] );
	    if ((*evals)[nev + i] != zero) {
	      normA[i+1] = normA[i];
	      i = i+2;
	    } else {
	      i++;
	    }
	  }
	}
	if (MyProblem->IsSymmetric()) {
	  cout<<"Real Part"<<"\t"<<"Residual"<<endl;
	  cout<<"------------------------------------------------------"<<endl;
	  for (i=0; i<nev; i++) {
	    cout<< (*evals)[i] << "\t\t"<< normA[i] << endl;
	  }  
	  cout<<"------------------------------------------------------"<<endl;
	} else {
	  cout<<"Real Part"<<"\t"<<"Imag Part"<<"\t"<<"Residual"<<endl;
	  cout<<"------------------------------------------------------"<<endl;
	  for (i=0; i<nev; i++) {
	    cout<< (*evals)[i] << "\t\t" << (*evals)[nev + i] << "\t\t"<< normA[i] << endl;
	  }  
	  cout<<"------------------------------------------------------"<<endl;
	}	

	if (bindx) delete [] bindx;
	if (update) delete [] update;
	if (col_inds) delete [] col_inds;
	if (val) delete [] val;
	if (row_vals) delete [] row_vals;

  	return 0;

} // end BlockKrylovSchurEpetraExHb.cpp
