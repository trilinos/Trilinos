
// @HEADER
// ***********************************************************************
// 
//            Trilinos: An Object-Oriented Solver Framework
//                 Copyright (2001) Sandia Corporation
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

// Trilinos Tutorial
// -----------------
// Output a CRS matrix in MATLAB format
//
// (output reported at the end of the file)
//
// Marzio Sala, SNL, 9214, 20-Nov-2003

#include "Epetra_config.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"

// ============================================================
// define a class, derived from Epetra_CrsMatrix, which
// initializes the matrix entires. User has to provide
// a valid Epetra_Map in the constructor, plus the diagonal
// value, and the sub- and super-diagonal values.
// ============================================================

class TridiagonalCrsMatrix : public Epetra_CrsMatrix { 
  
public:
  TridiagonalCrsMatrix(const Epetra_Map & Map,
			      double a,
			      double diag, double c) :
    Epetra_CrsMatrix(Copy,Map,3) 
  {

    // global number of rows
    int NumGlobalElements = Map.NumGlobalElements();
    // local number of rows
    int NumMyElements = Map.NumMyElements();
    // get update list
    int * MyGlobalElements = new int [NumMyElements];
    Map.MyGlobalElements( MyGlobalElements );

    // Add  rows one-at-a-time
    // Need some vectors to help
    // Off diagonal Values will always be -1
    
    double *Values = new double[2];
    Values[0] = a; Values[1] = c;
    int *Indices = new int[2];
    int NumEntries;
    
    for( int i=0 ; i<NumMyElements; ++i ) {
      if (MyGlobalElements[i]==0) {
	Indices[0] = 1;
	NumEntries = 1;
      } else if (MyGlobalElements[i] == NumGlobalElements-1) {
	Indices[0] = NumGlobalElements-2;
	NumEntries = 1;
      } else {
	Indices[0] = MyGlobalElements[i]-1;
	Indices[1] = MyGlobalElements[i]+1;
	NumEntries = 2;
      }
      assert(this->InsertGlobalValues(MyGlobalElements[i], NumEntries, Values, Indices)==0);
      // Put in the diagonal entry
      assert(this->InsertGlobalValues(MyGlobalElements[i], 1, &diag, MyGlobalElements+i)==0);
    }

    TransformToLocal();
  }
  
};

/* ======== ================ *
 * function CrsMatrix2MATLAB *
 * ======== ================ *
 *
 * Print out a CrsMatrix in a MATLAB format. Each processor prints out
 * its part, starting from proc 0 to proc NumProc-1. The first line of
 * each processor's output states the number of local rows and of
 * local nonzero elements. Output is finished by "End of Matrix Output".
 *
 *
 * Return code:        true if matrix has been printed out
 * -----------         false otherwise
 *
 * Parameters:
 * ----------
 *
 * - Epetra_CrsMatrix  reference to the ditributed CrsMatrix to 
 *                     print out
 */

bool CrsMatrix2MATLAB( const Epetra_CrsMatrix & A ) 

{

  int MyPID = A.Comm().MyPID(); 
  int NumProc = A.Comm().NumProc();

  // work only on transformed matrices;
  if( A.IndicesAreLocal() == false ) {
    if( MyPID == 0 ) { 
      cerr << "ERROR in "<< __FILE__ << ", line " << __LINE__ << endl;
      cerr << "Function CrsMatrix2MATLAB accepts\n";
      cerr << "transformed matrices ONLY. Please call A.TransformToLoca()\n";
      cerr << "on your matrix A to that purpose.\n";
      cerr << "Now returning...\n";
    }
    return false;
  }

  int NumMyRows = A.NumMyRows(); // number of rows on this process
  int NumNzRow;   // number of nonzero elements for each row
  int NumEntries; // number of extracted elements for each row
  int NumGlobalRows; // global dimensio of the problem
  int GlobalRow;  // row in global ordering
  int NumGlobalNonzeros; // global number of nonzero elements

  NumGlobalRows = A.NumGlobalRows();
  NumGlobalNonzeros = A.NumGlobalNonzeros();

  // print out on cout if no filename is provided

  int IndexBase = A.IndexBase(); // MATLAB start from 0
  if( IndexBase == 0 ) IndexBase = 1; 

  // write on file the dimension of the matrix

  if( MyPID==0 ) {
    cout << "A = spalloc(";
    cout << NumGlobalRows << ',' << NumGlobalRows;
    cout << ',' << NumGlobalNonzeros << ");\n";
  }

  for( int Proc=0 ; Proc<NumProc ; ++Proc ) {

    if( MyPID == Proc ) {

      cout << "% On proc " << Proc << ": ";
      cout << NumMyRows << " rows and ";
      cout << A.NumMyNonzeros() << " nonzeros\n";

      // cycle over all local rows to find out nonzero elements
      for( int MyRow=0 ; MyRow<NumMyRows ; ++MyRow ) {
	
	GlobalRow = A.GRID(MyRow);
	
	NumNzRow = A.NumMyEntries(MyRow);
	double *Values = new double[NumNzRow];
	int *Indices = new int[NumNzRow];
	
	A.ExtractMyRowCopy(MyRow, NumNzRow, 
			   NumEntries, Values, Indices);
	// print out the elements with MATLAB syntax
	for( int j=0 ; j<NumEntries ; ++j ) {
	  cout << "A(" << GlobalRow  + IndexBase 
	       << "," << A.GCID(Indices[j]) + IndexBase
	       << ") = " << Values[j] << ";\n";
	}
	
	delete Values;
	delete Indices;
      }
      
    }
    A.Comm().Barrier();
    if( MyPID == 0 ) {
      cout << " %End of Matrix Output\n";
    }
  }

  return true;

}

// =========== //
// main driver //
// =========== //

int main(int argc, char *argv[]) {

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  int ierr;
  
  // set global dimension to 5, could be any number
  int NumGlobalElements = 5;

  // define a linear map
  Epetra_Map Map(NumGlobalElements,0,Comm);
  
  // create the matrix
  TridiagonalCrsMatrix A( Map, -1.0, 2.0, -1.0);

  // output informationto stdout
  CrsMatrix2MATLAB( A );

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return( EXIT_SUCCESS );

}

/*

Output of this program

[msala:epetra]> mpirun -np 2 ./ex15.exe
A = spalloc(5,5,13);
% On proc 0: 3 rows and 8 nonzeros
A(1,1) = 2;
A(1,2) = -1;
A(2,1) = -1;
A(2,2) = 2;
A(2,3) = -1;
A(3,2) = -1;
A(3,3) = 2;
A(3,4) = -1;
% On proc 1: 2 rows and 5 nonzeros
A(4,4) = 2;
A(4,5) = -1;
A(4,3) = -1;
A(5,4) = -1;
A(5,5) = 2;
*/
