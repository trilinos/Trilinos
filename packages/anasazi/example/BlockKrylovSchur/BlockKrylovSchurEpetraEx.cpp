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
//  This example computes the specified eigenvalues of the discretized 2D Convection-Diffusion
//  equation using the block Krylov-Schur method.  This discretized operator is constructed as an
//  Epetra matrix, then passed into the Anasazi::EpetraOp to be used in the construction of the
//  Krylov decomposition.  The specifics of the block Krylov-Schur method can be set by the user.

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
	int i;
	double zero = 0.0;

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
        int nx = 10;  			// Discretization points in any one direction.
	int NumGlobalElements = nx*nx;	// Size of matrix nx*nx

	// Construct a Map that puts approximately the same number of
	// equations on each processor.

	Epetra_Map Map(NumGlobalElements, 0, Comm);

	// Get update list and number of local equations from newly created Map.

	int NumMyElements = Map.NumMyElements();

	std::vector<int> MyGlobalElements(NumMyElements);
    	Map.MyGlobalElements(&MyGlobalElements[0]);

	// Create an integer vector NumNz that is used to build the Petra Matrix.
	// NumNz[i] is the Number of OFF-DIAGONAL term for the ith global equation
	// on this processor
	std::vector<int> NumNz(NumMyElements);

	/* We are building a matrix of block structure:
	
			| T -I          |
			|-I  T -I       |
			|   -I  T       |
			|        ...  -I|
			|           -I T|

	 where each block is dimension nx by nx and the matrix is on the order of
	 nx*nx.  The block T is a tridiagonal matrix. 
	*/

	for (i=0; i<NumMyElements; i++) {
		if (MyGlobalElements[i] == 0 || MyGlobalElements[i] == NumGlobalElements-1 || 
		    MyGlobalElements[i] == nx-1 || MyGlobalElements[i] == nx*(nx-1) )
			NumNz[i] = 3;
		else if (MyGlobalElements[i] < nx || MyGlobalElements[i] > nx*(nx-1) || 
                         MyGlobalElements[i]%nx == 0 || (MyGlobalElements[i]+1)%nx == 0)
			NumNz[i] = 4;
		else
			NumNz[i] = 5;
	}

	// Create an Epetra_Matrix

	Teuchos::RefCountPtr<Epetra_CrsMatrix> A = Teuchos::rcp( new Epetra_CrsMatrix(Copy, Map, &NumNz[0]) );

	// Diffusion coefficient, can be set by user.
	// When rho*h/2 <= 1, the discrete convection-diffusion operator has real eigenvalues.
	// When rho*h/2 > 1, the operator has complex eigenvalues.
	double rho = 0.0;  

	// Compute coefficients for discrete convection-diffution operator
	const double one = 1.0;
	std::vector<double> Values(4);
	std::vector<int> Indices(4);
	double h = one /(nx+1);
	double h2 = h*h;
	double c = 5.0e-01*rho/ h;
	Values[0] = -one/h2 - c; Values[1] = -one/h2 + c; Values[2] = -one/h2; Values[3]= -one/h2;
	double diag = 4.0 / h2;
	int NumEntries;
	
	for (i=0; i<NumMyElements; i++)
	{
		if (MyGlobalElements[i]==0)
		{
			Indices[0] = 1;
			Indices[1] = nx;
			NumEntries = 2;
			assert(A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[1], &Indices[0])==0);
		}
		else if (MyGlobalElements[i] == nx*(nx-1))
		{
			Indices[0] = nx*(nx-1)+1;
			Indices[1] = nx*(nx-2);
			NumEntries = 2;
			assert(A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[1], &Indices[0])==0);
		}
		else if (MyGlobalElements[i] == nx-1)
		{
			Indices[0] = nx-2;
			NumEntries = 1;
			assert(A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[0], &Indices[0])==0);
			Indices[0] = 2*nx-1;
			assert(A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[2], &Indices[0])==0);
		}
		else if (MyGlobalElements[i] == NumGlobalElements-1)
		{
			Indices[0] = NumGlobalElements-2;
			NumEntries = 1;
			assert(A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[0], &Indices[0])==0);
			Indices[0] = nx*(nx-1)-1;
			assert(A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[2], &Indices[0])==0);
		}
		else if (MyGlobalElements[i] < nx)
		{
                        Indices[0] = MyGlobalElements[i]-1;
                        Indices[1] = MyGlobalElements[i]+1;
			Indices[2] = MyGlobalElements[i]+nx;
                        NumEntries = 3;
                        assert(A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[0], &Indices[0])==0);
		}
		else if (MyGlobalElements[i] > nx*(nx-1))
		{
                        Indices[0] = MyGlobalElements[i]-1;
                        Indices[1] = MyGlobalElements[i]+1;
			Indices[2] = MyGlobalElements[i]-nx;
                        NumEntries = 3;
                        assert(A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[0], &Indices[0])==0);
		}
                else if (MyGlobalElements[i]%nx == 0)
		{
			Indices[0] = MyGlobalElements[i]+1;
			Indices[1] = MyGlobalElements[i]-nx;
			Indices[2] = MyGlobalElements[i]+nx;
			NumEntries = 3;
			assert(A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[1], &Indices[0])==0);
		}
		else if ((MyGlobalElements[i]+1)%nx == 0)
		{
			Indices[0] = MyGlobalElements[i]-nx;
                        Indices[1] = MyGlobalElements[i]+nx;
                        NumEntries = 2;
                        assert(A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[2], &Indices[0])==0);
                        Indices[0] = MyGlobalElements[i]-1;
                        NumEntries = 1;
                        assert(A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[0], &Indices[0])==0);
		}
		else
		{
			Indices[0] = MyGlobalElements[i]-1;
			Indices[1] = MyGlobalElements[i]+1;
			Indices[2] = MyGlobalElements[i]-nx;
			Indices[3] = MyGlobalElements[i]+nx;
			NumEntries = 4;
			assert(A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[0], &Indices[0])==0);
		}
		// Put in the diagonal entry
		assert(A->InsertGlobalValues(MyGlobalElements[i], 1, &diag, &MyGlobalElements[i])==0);
	}

	// Finish up
	assert(A->TransformToLocal()==0);
	A->SetTracebackMode(1); // Shutdown Epetra Warning tracebacks

	//************************************
	// Start the block Arnoldi iteration
	//***********************************
	//
	//  Variables used for the Block Krylov Schur Method
	//	  
	int nev = 4;
	int blockSize = 2;
	int maxBlocks = 20;
	int maxRestarts = 300;
	int step = 1;
	double tol = 1e-8;
	string which="SM";	
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
	typedef Anasazi::MultiVecTraits<double, Epetra_MultiVector> MVT;

	// Create an Epetra_MultiVector for an initial vector to start the solver.
	// Note:  This needs to have the same number of columns as the blocksize.
	Teuchos::RefCountPtr<Epetra_MultiVector> ivec = Teuchos::rcp( new Epetra_MultiVector(Map, blockSize) );
	ivec->Random();

	// Create the eigenproblem.
	Teuchos::RefCountPtr<Anasazi::BasicEigenproblem<double, MV, OP> > MyProblem =
	  Teuchos::rcp( new Anasazi::BasicEigenproblem<double, MV, OP>(A, ivec) );
	
	// Inform the eigenproblem that the operator A is symmetric
	MyProblem->SetSymmetric(rho==0.0); 
	
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

	// Initialize the Block Arnoldi solver
	Anasazi::BlockKrylovSchur<double, MV, OP> MySolver(MyProblem, MySort, MyOM, MyPL);
	
	// Solve the problem to the specified tolerances or length
	MySolver.solve();

	// Retrieve eigenvalues
	Teuchos::RefCountPtr<std::vector<double> > evals = MyProblem->GetEvals();

	// Retrieve eigenvectors
	// The size of the eigenvector storage is 2 x nev.  
	// The real part of the eigenvectors is stored in the first nev vectors.
	// The imaginary part of the eigenvectors is stored in the second nev vectors.
	Teuchos::RefCountPtr<Epetra_MultiVector> evecr, eveci;
	std::vector<int> index(nev);
	
	// Get real part.
	Teuchos::RefCountPtr<Epetra_MultiVector> evecs = MyProblem->GetEvecs();

	// Get imaginary part, if needed.
	if (MyProblem->IsSymmetric()) {
          evecr = evecs;
        }
        else {
          evecr = Teuchos::rcp( new Epetra_MultiVector( View, *evecs, 0, nev ) );
          eveci = Teuchos::rcp( new Epetra_MultiVector( View, *evecs, nev, nev ) );	  
	}	  
	
	// Compute residuals.
	Teuchos::LAPACK<int,double> lapack;
	Epetra_MultiVector tempAevec(Map,nev);
	Teuchos::SerialDenseMatrix<int,double> Breal(nev,nev), Breal2(nev,nev);
	Teuchos::SerialDenseMatrix<int,double> Bimag(nev,nev), Bimag2(nev,nev);
	std::vector<double> normA(nev);
	std::vector<double> tempnrm(nev);
	cout<<endl<< "Actual Residuals"<<endl;
	cout<<"------------------------------------------------------"<<endl;
	Breal.putScalar(0.0); 
	if (!MyProblem->IsSymmetric())
	  Bimag.putScalar(0.0);
	for (i=0; i<nev; i++) { 
	  normA[i] = 0.0;
	  Breal(i,i) = (*evals)[i]; 
	  if (!MyProblem->IsSymmetric())
	    Bimag(i,i) = (*evals)[nev+i]; 
	}
	A->Apply( *evecr, tempAevec );
	MVT::MvTimesMatAddMv( -1.0, *evecr, Breal, 1.0, tempAevec );
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
	    normA[i] /= lapack.LAPY2( (*evals)[i], (*evals)[nev+i] );
	    if ((*evals)[nev + i] != zero) {
	      normA[i+1] = normA[i];
	      i = i+2;
	    } else {
	      i++;
	    }
	  }
	}
	if (MyProblem->IsSymmetric()) {
	  cout<<"Real Part"<<"\t"<<"Direct Residual"<<endl;
	  cout<<"------------------------------------------------------"<<endl;
	  for (i=0; i<nev; i++) {
	    cout<< (*evals)[i] << "\t\t"<< normA[i] << endl;
	  }  
	  cout<<"------------------------------------------------------"<<endl;
	} else {
	  cout<<"Real Part"<<"\t"<<"Imag Part"<<"\t"<<"Direct Residual"<<endl;
	  cout<<"------------------------------------------------------"<<endl;
	  for (i=0; i<nev; i++) {
	    cout<< (*evals)[i] << "\t\t" << (*evals)[nev + i] << "\t\t"<< normA[i] << endl;
	  }  
	  cout<<"------------------------------------------------------"<<endl;
	}	

	return 0;
}
