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
//  equation using the block Arnoldi method.  This discretized operator is constructed as an
//  Epetra matrix, then passed into the AnasaziPetraMat to be used in the construction of the
//  Krylov decomposition.  The specifics of the block Arnoldi method can be set by the user.

#include "AnasaziPetraInterface.hpp"
#include "AnasaziBlockArnoldi.hpp"
#include "AnasaziConfigDefs.hpp"
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
	Epetra_MpiComm & Comm = *new Epetra_MpiComm(MPI_COMM_WORLD);
#else
	Epetra_SerialComm & Comm = *new Epetra_SerialComm();
#endif

	int MyPID = Comm.MyPID();
	int NumProc = Comm.NumProc();
	cout << "Processor "<<MyPID<<" of "<< NumProc << " is alive."<<endl;

	//bool verbose = (MyPID==0);

	//  Dimension of the matrix
        int nx = 10;  			// Discretization points in any one direction.
	int NumGlobalElements = nx*nx;	// Size of matrix nx*nx

	// Construct a Map that puts approximately the same number of
	// equations on each processor.

	Epetra_Map& Map = *new Epetra_Map(NumGlobalElements, 0, Comm);

	// Get update list and number of local equations from newly created Map.

	int NumMyElements = Map.NumMyElements();

	int * MyGlobalElements = new int[NumMyElements];
    	Map.MyGlobalElements(MyGlobalElements);

	// Create an integer vector NumNz that is used to build the Petra Matrix.
	// NumNz[i] is the Number of OFF-DIAGONAL term for the ith global equation
	// on this processor
	int * NumNz = new int[NumMyElements];

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

	Epetra_CrsMatrix& A = *new Epetra_CrsMatrix(Copy, Map, NumNz);

	// Diffusion coefficient, can be set by user.
	// When rho*h/2 <= 1, the discrete convection-diffusion operator has real eigenvalues.
	// When rho*h/2 > 1, the operator has complex eigenvalues.
	double rho = 0.0;  

	// Compute coefficients for discrete convection-diffution operator
	const double one = 1.0;
	double *Values = new double[4];
	double h = one /(nx+1);
	double h2 = h*h;
	double c = 5.0e-01*rho/ h;
	Values[0] = -one/h2 - c; Values[1] = -one/h2 + c; Values[2] = -one/h2; Values[3]= -one/h2;
	int *Indices = new int[4];
	double diag = 4.0 / h2;
	int NumEntries;
	
	for (i=0; i<NumMyElements; i++)
	{
		if (MyGlobalElements[i]==0)
		{
			Indices[0] = 1;
			Indices[1] = nx;
			NumEntries = 2;
			assert(A.InsertGlobalValues(MyGlobalElements[i], NumEntries, Values+1, Indices)==0);
		}
		else if (MyGlobalElements[i] == nx*(nx-1))
		{
			Indices[0] = nx*(nx-1)+1;
			Indices[1] = nx*(nx-2);
			NumEntries = 2;
			assert(A.InsertGlobalValues(MyGlobalElements[i], NumEntries, Values+1, Indices)==0);
		}
		else if (MyGlobalElements[i] == nx-1)
		{
			Indices[0] = nx-2;
			NumEntries = 1;
			assert(A.InsertGlobalValues(MyGlobalElements[i], NumEntries, Values, Indices)==0);
			Indices[0] = 2*nx-1;
			assert(A.InsertGlobalValues(MyGlobalElements[i], NumEntries, Values+2, Indices)==0);
		}
		else if (MyGlobalElements[i] == NumGlobalElements-1)
		{
			Indices[0] = NumGlobalElements-2;
			NumEntries = 1;
			assert(A.InsertGlobalValues(MyGlobalElements[i], NumEntries, Values, Indices)==0);
			Indices[0] = nx*(nx-1)-1;
			assert(A.InsertGlobalValues(MyGlobalElements[i], NumEntries, Values+2, Indices)==0);
		}
		else if (MyGlobalElements[i] < nx)
		{
                        Indices[0] = MyGlobalElements[i]-1;
                        Indices[1] = MyGlobalElements[i]+1;
			Indices[2] = MyGlobalElements[i]+nx;
                        NumEntries = 3;
                        assert(A.InsertGlobalValues(MyGlobalElements[i], NumEntries, Values, Indices)==0);
		}
		else if (MyGlobalElements[i] > nx*(nx-1))
		{
                        Indices[0] = MyGlobalElements[i]-1;
                        Indices[1] = MyGlobalElements[i]+1;
			Indices[2] = MyGlobalElements[i]-nx;
                        NumEntries = 3;
                        assert(A.InsertGlobalValues(MyGlobalElements[i], NumEntries, Values, Indices)==0);
		}
                else if (MyGlobalElements[i]%nx == 0)
		{
			Indices[0] = MyGlobalElements[i]+1;
			Indices[1] = MyGlobalElements[i]-nx;
			Indices[2] = MyGlobalElements[i]+nx;
			NumEntries = 3;
			assert(A.InsertGlobalValues(MyGlobalElements[i], NumEntries, Values+1, Indices)==0);
		}
		else if ((MyGlobalElements[i]+1)%nx == 0)
		{
			Indices[0] = MyGlobalElements[i]-nx;
                        Indices[1] = MyGlobalElements[i]+nx;
                        NumEntries = 2;
                        assert(A.InsertGlobalValues(MyGlobalElements[i], NumEntries, Values+2, Indices)==0);
                        Indices[0] = MyGlobalElements[i]-1;
                        NumEntries = 1;
                        assert(A.InsertGlobalValues(MyGlobalElements[i], NumEntries, Values, Indices)==0);
		}
		else
		{
			Indices[0] = MyGlobalElements[i]-1;
			Indices[1] = MyGlobalElements[i]+1;
			Indices[2] = MyGlobalElements[i]-nx;
			Indices[3] = MyGlobalElements[i]+nx;
			NumEntries = 4;
			assert(A.InsertGlobalValues(MyGlobalElements[i], NumEntries, Values, Indices)==0);
		}
		// Put in the diagonal entry
		assert(A.InsertGlobalValues(MyGlobalElements[i], 1, &diag, MyGlobalElements+i)==0);
	}

	// Finish up
	assert(A.TransformToLocal()==0);
	A.SetTracebackMode(1); // Shutdown Epetra Warning tracebacks

	//************************************
	// Start the block Arnoldi iteration
	//***********************************
	//
	//  Variables used for the Block Arnoldi Method
	//
	int block = 1;
	int length = 20;
	int nev = 4;
	double tol = 1.0e-8;
	string which="SM";
	int restarts = 300;
	int step = 1;
	//int step = restarts*length*block;

	// Create a PetraAnasaziVec. Note that the decision to make a view or
	// or copy is determined by the petra constructor called by Anasazi::PetraVec.
	// This is possible because I pass in arguements needed by petra.
	Anasazi::PetraVec ivec(Map, block);
	ivec.MvRandom();

	// Call the ctor that calls the petra ctor for a matrix
	Anasazi::PetraMat Amat(A);	
	Anasazi::Eigenproblem<double> MyProblem(&Amat, &ivec);

	// Inform the eigenproblem that the matrix A is symmetric
	MyProblem.SetSymmetric(rho==0.0); 

	// Initialize the Block Arnoldi solver
	Anasazi::BlockArnoldi<double> MyBlockArnoldi(MyProblem, tol, nev, length, block, 
						which, step, restarts);	
	MyBlockArnoldi.setDebugLevel(0);

#ifdef UNIX
	Epetra_Time & timer = *new Epetra_Time(Comm);
#endif

	// Iterate a few steps (if you wish)
	//MyBlockArnoldi.iterate(5);

	// Solve the problem to the specified tolerances or length
	MyBlockArnoldi.solve();

#ifdef UNIX
	double elapsed_time = timer.ElapsedTime();
	double total_flops = A.Flops();
	double MFLOPs = total_flops/elapsed_time/1000000.0;
#endif

	// Obtain results directly
	double* resids = MyBlockArnoldi.getResiduals();
	double* evalr = MyBlockArnoldi.getEvals(); 
	double* evali = MyBlockArnoldi.getiEvals();

	// Retrieve eigenvectors
	Anasazi::PetraVec evecr(Map, nev);
	MyBlockArnoldi.getEvecs( evecr );
	Anasazi::PetraVec eveci(Map, nev);
	MyBlockArnoldi.getiEvecs( eveci );

	// Output results to screen
	MyBlockArnoldi.currentStatus();
	
	// Compute residuals.
	Teuchos::LAPACK<int,double> lapack;
	Anasazi::PetraVec tempevecr(Map,nev), tempeveci(Map,nev);
	Anasazi::PetraVec tempAevec(Map,nev);
	Teuchos::SerialDenseMatrix<int,double> Breal(nev,nev), Bimag(nev,nev);
	Teuchos::SerialDenseMatrix<int,double> Breal2(nev,nev), Bimag2(nev,nev);
	double* normA = new double[nev];
	double* tempnrm = new double[nev];
	cout<<endl<< "Actual Residuals"<<endl;
	cout<<"------------------------------------------------------"<<endl;
	Breal.putScalar(0.0); Bimag.putScalar(0.0);
	for (i=0; i<nev; i++) { Breal(i,i) = evalr[i]; Bimag(i,i) = evali[i]; }
	Amat.ApplyMatrix( evecr, tempAevec );
	tempAevec.MvTimesMatAddMv( -1.0, evecr, Breal, 1.0 );
	tempAevec.MvTimesMatAddMv( 1.0, eveci, Bimag, 1.0 );
	tempAevec.MvNorm( normA );
	Amat.ApplyMatrix( eveci, tempAevec );
	tempAevec.MvTimesMatAddMv( -1.0, evecr, Bimag, 1.0 );
	tempAevec.MvTimesMatAddMv( -1.0, eveci, Breal, 1.0 );
	tempAevec.MvNorm( tempnrm );
	i = 0;
	while (i < nev) {
	  normA[i] = lapack.LAPY2( normA[i], tempnrm[i] );
	  normA[i] /= lapack.LAPY2( evalr[i], evali[i] );
	  if (evali[i] != zero) {
	    normA[i+1] = normA[i];
	    i = i+2;
	  } else {
	    i++;
	  }
	}
	cout<<"Real Part"<<"\t"<<"Imag Part"<<"\t"<<"Residual"<<endl;
	cout<<"------------------------------------------------------"<<endl;
	for (i=0; i<nev; i++) {
	  cout<< evalr[i] << "\t\t" << evali[i] << "\t\t"<< normA[i] << endl;
	}  
	cout<<"------------------------------------------------------"<<endl;
	
#ifdef UNIX
	if (verbose)
		cout << "\n\nTotal MFLOPs for Arnoldi = " << MFLOPs << " Elapsed Time = "<<  elapsed_time << endl;
#endif


	// Release all objects
	delete [] resids;
	delete [] evalr;
	delete [] evali;
	delete [] NumNz;
	delete [] Values;
	delete [] Indices;
	delete [] MyGlobalElements;

	delete &A;
	delete &Map;
	delete &Comm;

	return 0;
}
