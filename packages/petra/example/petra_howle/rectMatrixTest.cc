#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#ifdef PETRA_MPI
#include "mpi.h"
#endif
#ifndef __cplusplus
#define __cplusplus
#endif
#include "Petra_Comm.h"
#include "Petra_Map.h"
#include "Petra_Time.h"
#include "Petra_RDP_MultiVector.h"
#include "Petra_RDP_Vector.h"
#include "Petra_RDP_CRS_Matrix.h"

// Test code to make a rectangular petra matrix from data in a file
// -- vh

// prototypes
Petra_RDP_CRS_Matrix* readMatrixIn(FILE *dataFile, Petra_Comm& Comm);
Petra_RDP_CRS_Matrix* readRectMatrixIn(FILE *dataFile, Petra_Comm& Comm);
Petra_RDP_Vector* readVectorIn(FILE *dataFile, Petra_Comm& Comm);
void matVecTest(Petra_RDP_CRS_Matrix* Aptr,
								Petra_RDP_Vector* xptr,
								Petra_RDP_Vector* bptr);


int main(int argc, char *argv[])
{
	int ierr = 0, i, j;
	
#ifdef PETRA_MPI
  // Initialize MPI
  MPI_Init(&argc,&argv);
  int size, rank; // Number of MPI processes, My process ID
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
  int size = 1; // Serial case (not using MPI)
  int rank = 0;
#endif

#ifdef PETRA_MPI
  Petra_Comm & Comm = *new Petra_Comm( MPI_COMM_WORLD );
#else
  Petra_Comm & Comm = *new Petra_Comm();
#endif

  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();
  // cout << "Processor "<<MyPID<<" of "<< NumProc << " is alive."<<endl;

  bool verbose = (MyPID==0);


	// read in matrices:
  FILE *dataFile;

	// This is a square matrix
	cout << " reading in F matrix " << endl;
	dataFile = fopen("F.data","r");
  Petra_RDP_CRS_Matrix* Fptr = readMatrixIn(dataFile, Comm);
	fclose(dataFile);

	// Read in my rectangular matrix
	cout << " reading in B matrix " << endl;
	dataFile = fopen("B.data","r");
  Petra_RDP_CRS_Matrix* Bptr = readRectMatrixIn(dataFile, Comm);
	fclose(dataFile);
	Petra_RDP_CRS_Matrix& B = *Bptr;
	cout << "global rows " << B.NumGlobalRows() << endl;
	cout << "global cols " << B.NumGlobalCols() << endl;
	cout << "my local rows " << B.NumMyRows() << endl;
	cout << "my local cols  " << B.NumMyCols() << endl;

	// read in vectors (b and soln)
  cout << "reading in vector b " << endl;
  dataFile = fopen("rhs.data","r");
  Petra_RDP_Vector* bptr = readVectorIn(dataFile, Comm);
  fclose(dataFile);

}


Petra_RDP_CRS_Matrix* readMatrixIn(FILE *dataFile, Petra_Comm& Comm)
{
	/// Read in the matrix B (16x24 in little example)
	/// Assuming following file format:
  ///   line 1 Number of rows in F
	///   line 2: max of elements in vector NumNz
  ///   line 3:numnz+2  %d  (number of nonzeros per row)
	///   next line is total number of nonzeros
	///   rest is %d %d %g (i,j,s) matlab sparse format

	int NumGlobalElements = 0;
	int maxNumNz = 0;
	int i, j;
	// dataFile = fopen("B.data","r");
	fscanf(dataFile, "%d", &NumGlobalElements);
	// cout << "NumGlobalElements = " << NumGlobalElements << "\n" << endl;
	fscanf(dataFile, "%d", &maxNumNz);

  // Construct a Map that puts approximately the same number of 
  // equations on each processor.
  Petra_Map& Map = *new Petra_Map(NumGlobalElements, 0, Comm);

  // Get update list and number of local equations from newly created Map.
  int NumMyElements = Map.NumMyElements();

	int * MyGlobalElements = new int[NumMyElements];
	Map.MyGlobalElements(MyGlobalElements);

	int * NumNz = new int[NumMyElements];
  for (i=0; i<NumMyElements; i++)
		{
			fscanf(dataFile, "%d", &NumNz[i]);
			// NumNz[i] = NumNz[i] - 1; // subtracting off one for the diagonals
			// figure out what to do for this in the non-square matrices (B)
			// cout << NumNz[i] << endl;
		}

	Petra_RDP_CRS_Matrix* Mptr = new Petra_RDP_CRS_Matrix(Copy, Map, NumNz);
  Petra_RDP_CRS_Matrix& M = *Mptr;
  double *Values = new double[maxNumNz];
  int *Indices = new int[maxNumNz];
  int NumEntries;
	/// read in i,j,s values for B matrix -- vh
	/// next line is number of nonzeros
	/// rest is %d %d %g
	int nnzM = 0; /// number of nonzeros in matrix and length of i,j,s
	
	fscanf(dataFile, "%d", &nnzM);
	// cout << "\nnumber of nonzeros in B: " << nnzB << "\n" << endl;


	int * iM = new int[nnzM];
	int * jM = new int[nnzM];
	double * sM = new double[nnzM];
	for (i=0; i<nnzM; i++)
		{
			fscanf(dataFile, "%d %d %lg", &iM[i], &jM[i], &sM[i]);
			// matlab used indexing from 1, C++ starts at 0
			// so subtract 1:
			iM[i]--;
			jM[i]--;
			// cout << "iM: " << iM[i]
			// << "\tjM: " << jM[i]
			// << "\tsM: " << sM[i]
			// << endl;
		}

	/// now fill in the matrix values row by row
	int offset = 0;
	for (i=0; i<NumGlobalElements; i++) /// i.e., for each row...
		{
			// cout << "NumNz[" << i << "] = " << NumNz[i] << endl;
			for (j=0; j<NumNz[i]; j++)
				{
					/// set number of entries, values, column indices
					Indices[j] = jM[offset + j];
					Values[j] = sM[offset + j];
					// cout << "iM: " << iM[offset + j]
					//		 << "\tjM: " << jM[offset + j]
					//		 << "\tsM: " << sM[offset + j]
					//		 << endl;
				}
			NumEntries = NumNz[i];
			assert(M.InsertGlobalValues(MyGlobalElements[i], 
						    NumEntries, Values, Indices)==0);
			// cout << "offset = " << offset << endl;
			offset = offset + NumNz[i];
			// cout << "Got to here in B matrix." <<endl;
		}

  // Finish up
  assert(M.TransformToLocal()==0);

	cout << "nonzeros = " << M.NumGlobalNonzeros() << endl;
	cout << "rows = " << M.NumGlobalRows() << endl;
	cout << "cols = " << M.NumGlobalCols() << endl;
	return Mptr;
}


Petra_RDP_CRS_Matrix* readRectMatrixIn(FILE *dataFile, Petra_Comm& Comm)
{
	/// Read in the rectangular matrix 
	/// Assuming following file format:
  ///   line 1 Number of rows
	///   line 2: max of elements in vector NumNz
  ///   line 3:numnz+2  %d  (number of nonzeros per row)
	///   next line is total number of nonzeros
	///   rest is %d %d %g (i,j,s) matlab sparse format

	int NumGlobalElements = 0;
	int maxNumNz = 0;
	int i, j;
	// dataFile = fopen("B.data","r");
	fscanf(dataFile, "%d", &NumGlobalElements);
	// cout << "NumGlobalElements = " << NumGlobalElements << "\n" << endl;
	fscanf(dataFile, "%d", &maxNumNz);

  // Construct a Map that puts approximately the same number of 
  // equations on each processor.
  Petra_Map& Map = *new Petra_Map(NumGlobalElements, 0, Comm);

	// make another map for the columns'
	// Note -- figure out how to pass in the number of columns
	Petra_Map& ColMap = *new Petra_Map(24, 0, Comm);

  // Get update list and number of local equations from newly created Map.
  int NumMyElements = Map.NumMyElements();

	int * MyGlobalElements = new int[NumMyElements];
	Map.MyGlobalElements(MyGlobalElements);

	int * NumNz = new int[NumMyElements];
  for (i=0; i<NumMyElements; i++)
		{
			fscanf(dataFile, "%d", &NumNz[i]);
			// NumNz[i] = NumNz[i] - 1; // subtracting off one for the diagonals
			// figure out what to do for this in the non-square matrices (B)
			// cout << NumNz[i] << endl;
		}

	Petra_RDP_CRS_Matrix* Mptr = new Petra_RDP_CRS_Matrix(Copy, Map, NumNz);
  Petra_RDP_CRS_Matrix& M = *Mptr;
  double *Values = new double[maxNumNz];
  int *Indices = new int[maxNumNz];
  int NumEntries;
	/// read in i,j,s values for B matrix -- vh
	/// next line is number of nonzeros
	/// rest is %d %d %g
	int nnzM = 0; /// number of nonzeros in matrix and length of i,j,s
	
	fscanf(dataFile, "%d", &nnzM);
	// cout << "\nnumber of nonzeros in B: " << nnzB << "\n" << endl;


	int * iM = new int[nnzM];
	int * jM = new int[nnzM];
	double * sM = new double[nnzM];
	for (i=0; i<nnzM; i++)
		{
			fscanf(dataFile, "%d %d %lg", &iM[i], &jM[i], &sM[i]);
			// matlab used indexing from 1, C++ starts at 0
			// so subtract 1:
			iM[i]--;
			jM[i]--;
			// cout << "iM: " << iM[i]
			// << "\tjM: " << jM[i]
			// << "\tsM: " << sM[i]
			// << endl;
		}

	/// now fill in the matrix values row by row
	int offset = 0;
	for (i=0; i<NumGlobalElements; i++) /// i.e., for each row...
		{
			// cout << "NumNz[" << i << "] = " << NumNz[i] << endl;
			for (j=0; j<NumNz[i]; j++)
				{
					/// set number of entries, values, column indices
					Indices[j] = jM[offset + j];
					Values[j] = sM[offset + j];
					// cout << "iM: " << iM[offset + j]
					//		 << "\tjM: " << jM[offset + j]
					//		 << "\tsM: " << sM[offset + j]
					//  	 << endl;
				}
			NumEntries = NumNz[i];
			assert(M.InsertGlobalValues(MyGlobalElements[i], 
						    NumEntries, Values, Indices)==0);
			// cout << "offset = " << offset << endl;
			offset = offset + NumNz[i];
			// cout << "Got to here in B matrix." <<endl;
		}

  // Finish up
  assert(M.TransformToLocal(&ColMap, &Map)==0);

	cout << "nonzeros = " << M.NumGlobalNonzeros() << endl;
	cout << "rows = " << M.NumGlobalRows() << endl;
	cout << "cols = " << M.NumGlobalCols() << endl;
	return Mptr;
}





Petra_RDP_Vector* readVectorIn(FILE *dataFile, Petra_Comm& Comm)
{
	/// Read in a vector.
	/// Assuming following file format:
  ///   line 1 Length of (dense) vector
	///   line 2 is number of nonzero elements in vector
	///   rest is %d %g (i, s) indices/values of non-zeros

	int NumGlobalElements = 0;
	int NzElms = 0;
	int i;

	fscanf(dataFile, "%d", &NumGlobalElements);
	fscanf(dataFile, "%d", &NzElms);
  // Construct a Map that puts approximately the same number of 
  // equations on each processor.
  Petra_Map& Map = *new Petra_Map(NumGlobalElements, 0, Comm);

  // Get update list and number of local equations from newly created Map.
  int NumMyElements = Map.NumMyElements();
	int * MyGlobalElements = new int[NumMyElements];
	Map.MyGlobalElements(MyGlobalElements);

	// make a petra map filled with zeros
	Petra_RDP_Vector* vptr = new Petra_RDP_Vector(Map);
	Petra_RDP_Vector& v = *vptr;

	// now fill in the nonzero elements
	double * myArray = new double[NumMyElements];
	int tempInd;
	double tempVal;
	// cout << "Length v " << NumGlobalElements << endl;
	// cout << "NzElms " << NzElms << endl;
  for (i=0; i<NzElms; i++)
		{
			fscanf(dataFile, "%d %lg", &tempInd, &tempVal);
			v[tempInd] = tempVal;
			// cout << tempVal << endl;
		}
	//  Petra_RDP_CRS_Matrix& M = *new Petra_RDP_CRS_Matrix(Copy, Map, NumNz);
 
	return vptr;

}

// small matrix vector multiply test
void matVecTest(Petra_RDP_CRS_Matrix* Aptr,
								Petra_RDP_Vector* xptr,
								Petra_RDP_Vector* bptr)
{
	Petra_RDP_CRS_Matrix& A = *Aptr;
	Petra_RDP_Vector& x = *xptr;
	Petra_RDP_Vector& b = *bptr; // we're going to overwrite b anyway
	// will that overwrite x, too? Look at ExtractCopy for alternative.

  A.Multiply(1, x, b);
}



// matrix vector multiply for our block matrix
/************
Petra_RDP_Vector* = myMatVecMult(Petra_RDP_CRS_Matrix* Bptr, 
																 Petra_RDP_CRS_Matrix* Fptr,
																 Petra_RDP_CRS_Matrix* Cptr,
																 Petra_RDP_Vector* xptr)

{
	// A = [F B'; B C]
	// return Ax pointer
  // cout << "block matrix-vector multiply" << endl;

}
*************/
