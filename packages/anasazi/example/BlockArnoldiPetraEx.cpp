//
//  File: BlockArnoldiPetraEx.cpp
//
//  This example computes the specified eigenvalues of the discretized 2D Convection-Diffusion
//  equation using the block Arnoldi method.  This discretized operator is constructed as an
//  Epetra matrix, then passed into the AnasaziPetraMat to be used in the construction of the
//  Krylov decomposition.  The specifics of the block Arnoldi method can be set by the user.

#include "AnasaziPetraInterface.hpp"
#include "AnasaziBlockArnoldi.hpp"
#include "AnasaziConfigDefs.hpp"
#include "Epetra_CrsMatrix.h"

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"

int main(int argc, char *argv[]) {
	int ierr = 0, i, j;

#ifdef EPETRA_MPI

	// Initialize MPI

	MPI_Init(&argc,&argv);

	int size, rank; // Number of MPI processes, My process ID

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

#else

	int size = 1; // Serial case (not using MPI)
	int rank = 0;

#endif


#ifdef EPETRA_MPI
	Epetra_MpiComm & Comm = *new Epetra_MpiComm(MPI_COMM_WORLD);
#else
	Epetra_SerialComm & Comm = *new Epetra_SerialComm();
#endif

	int MyPID = Comm.MyPID();
	int NumProc = Comm.NumProc();
	cout << "Processor "<<MyPID<<" of "<< NumProc << " is alive."<<endl;

	bool verbose = (MyPID==0);

	//  Dimension of the matrix
        int nx = 10;  			// Discretization points in any one direction.
	int NumGlobalElements = nx*nx;	// Size of matrix nx*nx

	// We will use zero based indices
	int IndexBase = 0;

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
	int block = 2;
	int length = 20;
	int nev = 4;
	double tol = 1.0e-14;
	string which="SM";
	int step = 1;
	int restarts = 300;

	// create a PetraAnasaziVec. Note that the decision to make a view or
	// or copy is determined by the petra constructor called by Anasazi::PetraVec.
	// This is possible because I pass in arguements needed by petra.
	Anasazi::PetraVec<double> ivec(Map, block);
	ivec.MvRandom();

	// call the ctor that calls the petra ctor for a matrix
	Anasazi::PetraMat<double> Amat(A);	
	Anasazi::Eigenproblem<double> MyProblem(&Amat, &ivec);

	// initialize the Block Arnoldi solver
	Anasazi::BlockArnoldi<double> MyBlockArnoldi(MyProblem, tol, nev, length, block, 
						which, step, restarts);
	
	// inform the solver that the problem is symmetric
	//MyBlockArnoldi.setSymmetric(true);
	MyBlockArnoldi.setDebugLevel(0);

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
	Anasazi::PetraVec<double> evecr(Map, nev);
	MyBlockArnoldi.getEvecs( evecr );
	Anasazi::PetraVec<double> eveci(Map, nev);
	MyBlockArnoldi.getiEvecs( eveci );

	// output results to screen
	MyBlockArnoldi.currentStatus();



#ifdef UNIX
	if (verbose)
		cout << "\n\nTotal MFLOPs for Arnoldi = " << MFLOPs << " Elapsed Time = "<<  elapsed_time << endl;
#endif


	// Release all objects
	delete [] resids, evalr, evali;
	delete [] NumNz;
	delete [] Values;
	delete [] Indices;
	delete [] MyGlobalElements;

	delete &A;
	delete &Map;
	delete &Comm;

	return 0;
}
