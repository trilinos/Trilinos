
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
// Print out some information about a CRS matrix
//
// (output reported at the end of the file)

#include "float.h"
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

    // Finish up
    assert(TransformToLocal()==0);

  }
  
};

// =============================================================
// This function print out some information about the input
// matrix. This includes number of rows and columns, some norms,
// few statistics about the nonzero structure of the matrix.
// Some details about the Trilinos storage of the matrix
// are also reported.
//
// Return code:        true if matrix has been printed out
// -----------         false otherwise
//
// Parameters:
// ----------
//
// - Epetra_CrsMatrix  reference to the ditributed CrsMatrix to 
//                     print out
// - os                output stream (can be cout)
//===============================================================

bool CrsMatrixInfo( const Epetra_CrsMatrix & A,
		    ostream & os ) 

{

  int MyPID = A.Comm().MyPID(); 
  int NumProc = A.Comm().NumProc();

  // take care that matrix is already trasformed
  bool IndicesAreGlobal = A.IndicesAreGlobal();
  if( IndicesAreGlobal == true ) {
    if( MyPID == 0 ) {
      os << "WARNING : matrix must be transformed to local\n";
      os << "          before calling CrsMatrixInfo\n";
      os << "          Now returning...\n";
    }
    return false;
  }

  int NumGlobalRows = A.NumGlobalRows();
  int NumGlobalNonzeros = A.NumGlobalNonzeros();
  int NumGlobalCols = A.NumGlobalCols();
  double NormInf = A.NormInf();
  double NormOne = A.NormOne();
  int NumGlobalDiagonals = A.NumGlobalDiagonals();
  int GlobalMaxNumEntries = A.GlobalMaxNumEntries();
  int IndexBase = A.IndexBase();
  bool StorageOptimized = A.StorageOptimized();
  bool LowerTriangular = A.LowerTriangular();
  bool UpperTriangular = A.UpperTriangular();
  bool NoDiagonal = A.NoDiagonal();
  bool Sorted = A.Sorted();

  // these variables identifies quantities I have to compute,
  // since not provided by Epetra_CrsMatrix

  double MyFrobeniusNorm( 0.0 ), FrobeniusNorm( 0.0 );
  double MyMinElement( DBL_MAX ), MinElement( DBL_MAX );
  double MyMaxElement( DBL_MIN ), MaxElement( DBL_MIN );
  double MyMinAbsElement( DBL_MAX ), MinAbsElement( DBL_MAX );
  double MyMaxAbsElement( 0.0 ), MaxAbsElement( 0.0 );
  double MinDiagonalElement( DBL_MAX );
  double MaxDiSagonalElement( DBL_MIN );
  int NumDiagonalDominantRows( 0 );
  int MaxBandwidth;

  int NumMyRows = A.NumMyRows();
  int * NzPerRow = new int[NumMyRows];
  int Row; // iterator on rows
  int Col; // iterator on cols
  int MaxNumEntries = A.MaxNumEntries();
  double * Values = new double[MaxNumEntries];
  int * Indices = new int[MaxNumEntries];
  double Element, AbsElement; // generic nonzero element and its abs value
  int NumEntries;
  double * Diagonal = new double [NumMyRows];
  // SumOffDiagonal is the sum of absolute values for off-diagonals
  double * SumOffDiagonal = new double [NumMyRows];  
  for( Row=0 ;  Row<NumMyRows ; ++Row ) {
    SumOffDiagonal[Row] = 0.0;
  }
  int * IsDiagonallyDominant = new int [NumMyRows];
  int GlobalRow;

  // cycle over all matrix elements
  for( Row=0 ; Row<NumMyRows ; ++Row ) {
    GlobalRow = A.GRID(Row);
    NzPerRow[Row] = A.NumMyEntries(Row);
    A.ExtractMyRowCopy(Row,NzPerRow[Row],NumEntries,Values,Indices);
    for( Col=0 ; Col<NumEntries ; ++Col ) {
      Element = Values[Col];
      AbsElement = abs(Element);
      if( Element<MyMinElement ) MyMinElement = Element;
      if( Element>MyMaxElement ) MyMaxElement = Element;
      if( AbsElement<MyMinAbsElement ) MyMinAbsElement = AbsElement;
      if( AbsElement>MyMaxAbsElement ) MyMaxAbsElement = AbsElement;
      if( Indices[Col] == Row ) Diagonal[Row] = Element;
      else
	SumOffDiagonal[Row] += abs(Element);
      MyFrobeniusNorm += pow(Element,2);
    }
  }   

  // analise storage per row 
  int MyMinNzPerRow( NumMyRows ), MinNzPerRow( NumMyRows );
  int MyMaxNzPerRow( 0 ), MaxNzPerRow( 0 );

  for( Row=0 ; Row<NumMyRows ; ++Row ) {
    if( NzPerRow[Row]<MyMinNzPerRow ) MyMinNzPerRow=NzPerRow[Row];
    if( NzPerRow[Row]>MyMaxNzPerRow ) MyMaxNzPerRow=NzPerRow[Row];
  }

  // a test to see if matrix is diagonally-dominant

  int MyDiagonalDominance( 0 ), DiagonalDominance( 0 );
  int MyWeakDiagonalDominance( 0 ), WeakDiagonalDominance( 0 );

  for( Row=0 ; Row<NumMyRows ; ++Row ) {
    if( abs(Diagonal[Row])>SumOffDiagonal[Row] ) 
      ++MyDiagonalDominance;
    else if( abs(Diagonal[Row])==SumOffDiagonal[Row] ) 
      ++MyWeakDiagonalDominance;
  }

  // reduction operations
  
  A.Comm().SumAll(&MyFrobeniusNorm, &FrobeniusNorm, 1);
  A.Comm().MinAll(&MyMinElement, &MinElement, 1);
  A.Comm().MaxAll(&MyMaxElement, &MaxElement, 1);
  A.Comm().MinAll(&MyMinAbsElement, &MinAbsElement, 1);
  A.Comm().MaxAll(&MyMaxAbsElement, &MaxAbsElement, 1);
  A.Comm().MinAll(&MyMinNzPerRow, &MinNzPerRow, 1);
  A.Comm().MaxAll(&MyMaxNzPerRow, &MaxNzPerRow, 1);
  A.Comm().SumAll(&MyDiagonalDominance, &DiagonalDominance, 1);
  A.Comm().SumAll(&MyWeakDiagonalDominance, &WeakDiagonalDominance, 1);

  // free memory

  delete Values;
  delete Indices;
  delete Diagonal;
  delete SumOffDiagonal;
  delete IsDiagonallyDominant;
  delete NzPerRow;

  // simply no output for MyPID>0, only proc 0 write on os
  if( MyPID != 0 ) return true;

  os << "*** general Information about the matrix\n";
  os << "Number of Global Rows = " << NumGlobalRows << endl;
  os << "Number of Global Cols = " << NumGlobalCols << endl;
  os << "is the matrix square  = " <<
    ((NumGlobalRows==NumGlobalCols)?"yes":"no") << endl;
  os << "||A||_\\infty          = " << NormInf << endl;
  os << "||A||_1               = " << NormOne << endl;
  os << "||A||_F               = " << sqrt(FrobeniusNorm) << endl;
  os << "Number of nonzero diagonal entries = "
     << NumGlobalDiagonals
     << "( " << 1.0* NumGlobalDiagonals/NumGlobalRows*100
     << " %)\n";
  os << "Nonzero per row : min = " << MinNzPerRow 
     << " average = " << 1.0*NumGlobalNonzeros/NumGlobalRows
     << " max = " << MaxNzPerRow << endl; 
  os << "Maximum number of nonzero elements/row = " 
     << GlobalMaxNumEntries << endl;
  os << "min( a_{i,j} )      = " << MinElement << endl;
  os << "max( a_{i,j} )      = " << MaxElement << endl;
  os << "min( abs(a_{i,j}) ) = " << MinAbsElement << endl;
  os << "max( abs(a_{i,j}) ) = " << MaxAbsElement << endl;
  os << "Number of diagonal dominant rows        = " << DiagonalDominance 
     << " (" << 100.0*DiagonalDominance/NumGlobalRows << " % of total)\n";
  os << "Number of weakly diagonal dominant rows = " 
     << WeakDiagonalDominance 
     << " (" << 100.0*WeakDiagonalDominance/NumGlobalRows << " % of total)\n";

  os << "*** Information about the Trilinos storage\n";
  os << "Base Index                 = " << IndexBase << endl;
  os << "is storage optimized       = " 
     << ((StorageOptimized==true)?"yes":"no") << endl;
  os << "are indices global         = "
     << ((IndicesAreGlobal==true)?"yes":"no") << endl;
  os << "is matrix lower triangular = " 
     << ((LowerTriangular==true)?"yes":"no") << endl;
  os << "is matrix upper triangular = " 
     << ((UpperTriangular==true)?"yes":"no") << endl;
  os << "are there diagonal entries = " 
     <<  ((NoDiagonal==false)?"yes":"no") << endl;
  os << "is matrix sorted           = " 
     <<  ((Sorted==true)?"yes":"no") << endl;

  return true;

}

// =========== //
// Main driver //
// =========== //

int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  // set global dimension to 5, could be any number
  int NumGlobalElements = 5;

  // create a map
  Epetra_Map Map(NumGlobalElements,0,Comm);

  TridiagonalCrsMatrix A( Map, -1.0, 2.0, -1.0);

  // query the matrix
  CrsMatrixInfo( A, cout );

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return( EXIT_SUCCESS );

}

/*

Output of this program (NOTE: the output produced by our code can be
slightly different)

[msala:epetra]> mpirun -np 2 ./ex14.exe
*** general Information about the matrix
Number of Global Rows = 5
Number of Global Cols = 5
is the matrix square  = yes
||A||_\infty          = 4
||A||_1               = 4
||A||_F               = 5.2915
Number of nonzero diagonal entries = 5( 100 %)
Nonzero per row : min = 2 average = 2.6 max = 3
Maximum number of nonzero elements/row = 3
min( a_{i,j} )      = -1
max( a_{i,j} )      = 2
min( abs(a_{i,j}) ) = 1
max( abs(a_{i,j}) ) = 2
Number of diagonal dominant rows        = 2 (40 % of total)
Number of weakly diagonal dominant rows = 3 (60 % of total)
*** Information about the Trilinos storage
Base Index                 = 0
is storage optimized       = no
are indices global         = no
is matrix lower triangular = no
is matrix upper triangular = no
are there diagonal entries = yes
is matrix sorted           = yes

*/
