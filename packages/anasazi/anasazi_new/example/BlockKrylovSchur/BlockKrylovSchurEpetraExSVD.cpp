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
//  This example shows how to use the block Krylov-Schur method to compute a few
//  of the largest singular values (sigma) and corresponding right singular 
//  vectors (v) for the matrix A by solving the symmetric problem:
//
//                             (A'*A)*v = sigma*v
//
//  where A is an m by n real matrix that is derived from the simplest finite
//  difference discretization of the 2-dimensional kernel K(s,t)dt where
//
//                    K(s,t) = s(t-1)   if 0 <= s <= t <= 1
//                             t(s-1)   if 0 <= t <= s <= 1
//
//  NOTE:  This example came from the ARPACK SVD driver dsvd.f
//
#include "AnasaziBlockKrylovSchur.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziBasicSort.hpp"
#include "AnasaziConfigDefs.hpp"
#include "AnasaziEpetraAdapter.hpp"

#include "Epetra_CrsMatrix.h"
#include "Teuchos_LAPACK.hpp"

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"

int main(int argc, char *argv[]) {
	int i, j;
	const double one = 1.0;
	const double zero = 0.0;
	Teuchos::LAPACK<int,double> lapack;

#ifdef EPETRA_MPI

	// Initialize MPI
	MPI_Init(&argc,&argv);

#endif


#ifdef EPETRA_MPI
	Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
	Epetra_SerialComm Comm;
#endif

	int MyPID = Comm.MyPID();
	int NumProc = Comm.NumProc();
	cout << "Processor "<<MyPID<<" of "<< NumProc << " is alive."<<endl;

	bool verbose = (MyPID==0);

	//  Dimension of the matrix
	int m = 500;
	int n = 100;

	// Construct a Map that puts approximately the same number of
	// equations on each processor.

	Epetra_Map RowMap(m, 0, Comm);
	Epetra_Map ColMap(n, 0, Comm);

	// Get update list and number of local equations from newly created Map.

	int NumMyRowElements = RowMap.NumMyElements();
	
	std::vector<int> MyGlobalRowElements(NumMyRowElements);
    	RowMap.MyGlobalElements(&MyGlobalRowElements[0]);

	/* We are building an m by n matrix with entries
	  
	            A(i,j) = k*(si)*(tj - 1) if i <= j
		           = k*(tj)*(si - 1) if i  > j
	
	where si = i/(m+1) and tj = j/(n+1) and k = 1/(n+1).
	*/

	// Create an Epetra_Matrix
	Epetra_CrsMatrix A(Copy, RowMap, n);

	// Compute coefficients for discrete integral operator
	std::vector<double> Values(n);
	std::vector<int> Indices(n);
	double inv_mp1 = one/(m+1);
	double inv_np1 = one/(n+1);
	for (i=0; i<n; i++) { Indices[i] = i; }
	
	for (i=0; i<NumMyRowElements; i++)	{
	  //
	  for (j=0; j<n; j++) {
	    //
	    if ( MyGlobalRowElements[i] <= j )
	      Values[j] = inv_np1 * ( (MyGlobalRowElements[i]+one)*inv_mp1 ) * ( (j+one)*inv_np1 - one );  // k*(si)*(tj-1)
	    else
	      Values[j] = inv_np1 * ( (j+one)*inv_np1 ) * ( (MyGlobalRowElements[i]+one)*inv_mp1 - one );  // k*(tj)*(si-1)
	  }
	  assert(A.InsertGlobalValues(MyGlobalRowElements[i], n, &Values[0], &Indices[0])==0);
	}

	// Finish up
	assert(A.TransformToLocal(&ColMap, &RowMap)==0);
	assert(A.OptimizeStorage()==0);
	A.SetTracebackMode(1); // Shutdown Epetra Warning tracebacks

	//************************************
	// Start the block Arnoldi iteration
	//***********************************
	//
	//  Variables used for the Block Arnoldi Method
	//
	int block = 1;
	int length = 10;
	int nev = 4;
	double tol = lapack.LAMCH('E');
	string which="LM";
	int restarts = 300;
	int step = restarts*length*block;

	typedef Anasazi::MultiVec<double> MV;
	typedef Anasazi::Operator<double> OP;

	// Create an Anasazi::EpetraVec for an initial vector to start the solver. 
	// Note:  This needs to have the same number of columns as the blocksize.
	Teuchos::RefCountPtr<Anasazi::EpetraVec> ivec = Teuchos::rcp( new Anasazi::EpetraVec(ColMap, block) );
	ivec->MvRandom();

	// Call the constructor for the (A^T*A) operator
	Teuchos::RefCountPtr<Anasazi::EpetraSymOp>  Amat = Teuchos::rcp( new Anasazi::EpetraSymOp(A) );	
	Teuchos::RefCountPtr<Anasazi::BasicEigenproblem<double, MV, OP> > MyProblem =
	  Teuchos::rcp( new Anasazi::BasicEigenproblem<double, MV, OP>(Amat, ivec) );

	// Inform the eigenproblem that the matrix A is symmetric
	MyProblem->SetSymmetric(true);

	// Set the number of eigenvalues requested and the blocksize the solver should use
	MyProblem->SetNEV( nev );
	MyProblem->SetBlockSize( block );

	// Inform the eigenproblem that you are finishing passing it information
	assert( MyProblem->SetProblem() == 0 );

	// Create a sorting manager to handle the sorting of eigenvalues in the solver
	Teuchos::RefCountPtr<Anasazi::BasicSort<double, MV, OP> > MySort = 
	  Teuchos::rcp( new Anasazi::BasicSort<double, MV, OP>(which) );

	// Create an output manager to handle the I/O from the solver
	Teuchos::RefCountPtr<Anasazi::OutputManager<double> > MyOM =
	  Teuchos::rcp( new Anasazi::OutputManager<double>( MyPID ) );
	//MyOM->SetVerbosity( 2 );	

	// Initialize the Block Arnoldi solver
	Anasazi::BlockKrylovSchur<double, MV, OP> MySolver(MyProblem, MySort, MyOM, tol, 
							   length, step, restarts);	
	
	// Iterate a few steps (if you wish)
	//MySolver.iterate(5);

	// Solve the problem to the specified tolerances or length
	MySolver.solve();

	// Obtain results directly
	double* evals = MyProblem->GetEvals();

	// Retrieve eigenvectors
	//	Anasazi::EpetraVec* evecr = dynamic_cast<Anasazi::EpetraVec*>(&(MyProblem->GetEvecs()));

	// Output results to screen
	MySolver.currentStatus();
	
	// Compute singular values/vectors and direct residuals.
	//
	// Compute singular values which are the square root of the eigenvalues
	if (MyOM->doOutput(-1)) {
	  cout<<"------------------------------------------------------"<<endl;
	  cout<<"Computed Singular Values: "<<endl;
	  cout<<"------------------------------------------------------"<<endl;
	}
	for (i=0; i<nev; i++) { evals[i] = Teuchos::ScalarTraits<double>::squareroot( evals[i] ); }
	//
	// Compute left singular vectors :  u = Av/sigma
	//
	std::vector<double> tempnrm(nev), directnrm(nev);
	std::vector<int> index(nev);
	for (i=0; i<nev; i++) { index[i] = i; }
	Anasazi::EpetraVec Av(RowMap,nev), u(RowMap,nev);
	Anasazi::EpetraVec* evecs = dynamic_cast<Anasazi::EpetraVec* >(MyProblem->GetEvecs()->CloneView( &index[0], nev ));
	Teuchos::SerialDenseMatrix<int,double> S(nev,nev);
        A.Apply( *evecs, Av );
	Av.MvNorm( &tempnrm[0] );
	for (i=0; i<nev; i++) { S(i,i) = one/tempnrm[i]; };
	u.MvTimesMatAddMv( one, Av, S, zero );
	//
	// Compute direct residuals : || Av - sigma*u ||
	//
	for (i=0; i<nev; i++) { S(i,i) = evals[i]; }
	Av.MvTimesMatAddMv( -one, u, S, one );
	Av.MvNorm( &directnrm[0] );
	if (MyOM->doOutput(-1)) {
	  cout<<"Singular Value"<<"\t\t"<<"Direct Residual"<<endl;
	  cout<<"------------------------------------------------------"<<endl;
	  for (i=0; i<nev; i++) {
	    cout<< evals[i] << "\t\t" << directnrm[i] << endl;
	  }  
	  cout<<"------------------------------------------------------"<<endl;
	}
	
	// Clean up.
	if (evecs) delete evecs;
	
	return 0;
}
