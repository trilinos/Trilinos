
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
// A simple distributed 2D finite element code (Laplacian)
// using the grid
//
//     3      4      5
//     o------o------o
//     | (0) /| (2) /|
//     |    / |    / |
//     |   /  |   /  |
//     |  /   |  /   |
//     | /    | /    |
//     |/ (1) |/ (3) |
//     o------o------o
//     0      1      2
//
// processor 0 hosts the nodes 0, 3,and 4, while processor 1
// hosts 1, 2, and 5.
// Note: The grid information are hardwired, but they can easily
// be replaced by I/O functions.
// This code must be run with two processes.

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
#include "Epetra_Import.h"
#include "Epetra_SerialDenseMatrix.h"

// function declaration
void compute_loc_matrix( double *x_triangle, double *y_triangle,
			 Epetra_SerialDenseMatrix &Ke );
int find( const int list[], const int length, const int index);

// main driver
int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  int NumMyElements; // NODES assigned to this processor
  int NumMyExternalElements; // nodes used by this proc, but not hosted
  int NumMyTotalElements;
  int FE_NumMyElements; // TRIANGLES assigned to this processor
  int * MyGlobalElements; // nodes assigned to this processor
  Epetra_IntSerialDenseMatrix T; // store the grid connectivity

  int MyPID=Comm.MyPID();

  cout << MyPID << endl;

  switch( MyPID ) {

  case 0:
    NumMyElements = 3;
    NumMyExternalElements = 2;
    NumMyTotalElements = NumMyElements + NumMyExternalElements;
    FE_NumMyElements = 3;

    MyGlobalElements = new int[NumMyTotalElements];
    MyGlobalElements[0] = 0;
    MyGlobalElements[1] = 4;
    MyGlobalElements[2] = 3;
    MyGlobalElements[3] = 1;
    MyGlobalElements[4] = 5;

    break;
  case 1:
    NumMyElements = 3;
    NumMyExternalElements = 2;
    NumMyTotalElements = NumMyElements + NumMyExternalElements;
    FE_NumMyElements = 3;

    MyGlobalElements = new int[NumMyTotalElements];
    MyGlobalElements[0] = 1;
    MyGlobalElements[1] = 2;
    MyGlobalElements[2] = 5;
    MyGlobalElements[3] = 0;
    MyGlobalElements[4] = 4;
    break;

  }

  // build Map corresponding to update
  Epetra_Map Map(-1,NumMyElements,MyGlobalElements,0,Comm);

  // vector containing coordinates BEFORE exchanging external nodes
  Epetra_Vector CoordX_noExt(Map);
  Epetra_Vector CoordY_noExt(Map);

  switch( MyPID ) {

  case 0:
    T.Shape(3,FE_NumMyElements);

    // fill x-coordinates
    CoordX_noExt[0] = 0.0; 
    CoordX_noExt[1] = 1.0; 
    CoordX_noExt[2] = 0.0;
    // fill y-coordinates
    CoordY_noExt[0] = 0.0; 
    CoordY_noExt[1] = 1.0; 
    CoordY_noExt[2] = 1.0;
    // fill connectivity
    T(0,0) = 0; T(0,1) = 4; T(0,2) = 3;
    T(1,0) = 0; T(1,1) = 1; T(1,2) = 4;
    T(2,0) = 4; T(2,1) = 1; T(2,2) = 5;
    break;
    
  case 1:

    T.Shape(3,FE_NumMyElements);

    // fill x-coordinates
    CoordX_noExt[0] = 1.0; 
    CoordX_noExt[1] = 2.0; 
    CoordX_noExt[2] = 2.0;
    // fill y-coordinates
    CoordY_noExt[0] = 0.0; 
    CoordY_noExt[1] = 0.0; 
    CoordY_noExt[2] = 1.0;
    // fill connectivity
    T(0,0) = 0; T(0,1) = 1; T(0,2) = 4;
    T(1,0) = 1; T(1,1) = 5; T(1,2) = 4;
    T(2,0) = 1; T(2,1) = 2; T(2,2) = 5;
    break;
  }

  // - - - - - - - - - - - - - - - - - - - - //
  // E X T E R N A L   N O D E S   S E T U P //
  // - - - - - - - - - - - - - - - - - - - - //

  // build target map to exchange the valus of external nodes
  Epetra_Map TargetMap(-1,NumMyTotalElements,
			  MyGlobalElements, 0, Comm);
  // !@# rename Map -> SourceMap ?????
  Epetra_Import Importer(TargetMap,Map);
  Epetra_Vector CoordX(TargetMap);
  Epetra_Vector CoordY(TargetMap);

  CoordX.Import(CoordX_noExt,Importer,Insert);
  CoordY.Import(CoordY_noExt,Importer,Insert);
 
  // now CoordX_noExt and CoordY_noExt are no longer required
  // NOTE: better to construct CoordX and CoordY as MultiVector

  // - - - - - - - - - - - - //
  // M A T R I X   S E T U P //
  // - - - - - - - - - - - - //

  // build the CRS matrix corresponding to the grid
  // some vectors are allocated
  const int MaxNnzRow = 5;

  Epetra_CrsMatrix A(Copy,Map,MaxNnzRow);

  int Element, MyRow, GlobalRow, GlobalCol, i, j, k;
  double zero=0.0;
  Epetra_IntSerialDenseMatrix Struct; // temp to create the matrix connectivity
  Struct.Shape(NumMyElements,MaxNnzRow);
  for( i=0 ; i<NumMyElements ; ++i ) 
    for( j=0 ; j<MaxNnzRow ; ++j )
      Struct(i,j) = -1;

  // cycle over all the finite elements
  for( Element=0 ; Element<FE_NumMyElements ; ++Element ) {
    // cycle over each row
    for( i=0 ; i<3 ; ++i ) {
      // get the global and local number of this row
      GlobalRow = T(Element,i);
      MyRow = A.LRID(GlobalRow);
      if( MyRow != -1 ) { // only rows stored on this proc
	// cycle over the columns
	for( j=0 ; j<3 ; ++j ) {
	  // get the global number only of this column
	  GlobalCol = T(Element,j);
	  // look if GlobalCol was already put in Struct
	  for( k=0 ; k<MaxNnzRow ; ++k ) {
	    if( Struct(MyRow,k) == GlobalCol ||
		Struct(MyRow,k) == -1 ) break; 
	  }
	  if( Struct(MyRow,k) == -1 ) { // new entry
	    Struct(MyRow,k) = GlobalCol;
	  } else if( Struct(MyRow,k) != GlobalCol ) {
	    // maybe not enough space has beenn allocated
	    cerr << "ERROR: not enough space for element "
		 << GlobalRow << "," << GlobalCol << endl;
	    return( 0 );
	  }
	}
      }
    }
  }

  int * Indices = new int [MaxNnzRow];
  double * Values  = new double [MaxNnzRow];
  for( i=0 ; i<MaxNnzRow ; ++i ) Values[i] = 0.0;

  // now use Struct to fill build the matrix structure
  for( int Row=0 ; Row<NumMyElements ; ++Row ) {
    int Length = 0;
    for( int j=0 ; j<MaxNnzRow ; ++j ) {
      if( Struct(Row,j) == -1 ) break;
      Indices[Length] = Struct(Row,j);
      Length++;
    }
    GlobalRow = MyGlobalElements[Row];    
    A.InsertGlobalValues(GlobalRow, Length, Values, Indices);
  }

  // replace global numbering with local one in T
  for( int Element=0 ; Element<FE_NumMyElements ; ++Element ) {
    for( int i=0 ; i<3 ; ++i ) {
      int global = T(Element,i);
      int local = find(MyGlobalElements,NumMyTotalElements,
			global);
      if( global == -1 ) {
	cerr << "ERROR\n";
	return( EXIT_FAILURE );
      }
      T(Element,i) = local;
    }
  }

  // - - - - - - - - - - - - - - //
  // M A T R I X   F I L L - I N //
  // - - - - - - - - - - - - - - //

  // room for the local matrix
  Epetra_SerialDenseMatrix Ke;
  Ke.Shape(3,3);

  // now fill the matrix
  for(  int Element=0 ; Element<FE_NumMyElements ; ++Element ) {
    // variables used inside
    int GlobalRow;
    int MyRow;
    int GlobalCol;
    double x_triangle[3];
    double y_triangle[3];
    // get the spatial coordinate of each local node
    for( int i=0 ; i<3 ; ++i ) {
      MyRow = T(Element,i);
      y_triangle[i] = CoordX[MyRow]; 
      x_triangle[i] = CoordY[MyRow]; 
    }
    // compute the local matrix for Element

    compute_loc_matrix( x_triangle, y_triangle,Ke ); 

    // insert it in the global one
    // cycle over each row
    for( int i=0 ; i<3 ; ++i ) {
      // get the global and local number of this row
      MyRow = T(Element,i);
      if( MyRow < NumMyElements ) {
	for( int j=0 ; j<3 ; ++j ) {
	  // get global column number
	  GlobalRow = MyGlobalElements[MyRow];
	  GlobalCol = MyGlobalElements[T(Element,j)];
	  A.SumIntoGlobalValues(GlobalRow,1,&(Ke(i,j)),&GlobalCol);
	}
      }
    }
  }

  A.TransformToLocal();

  // - - - - - - - - - - - - - //
  // R H S  &  S O L U T I O N //
  // - - - - - - - - - - - - - //

  Epetra_Vector x(Map), b(Map);
  x.Random(); b.PutScalar(0.0);

  // Solution can be obtained using Aztecoo

  // free memory before leaving
  delete MyGlobalElements;
  delete Indices;
  delete Values;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return( EXIT_SUCCESS );

} /* main */

int find( const int list[], const int length, const int index) 
{

  int pos=-1;
  for( int i=0 ; i<length ; ++i )
    if( list[i] == index ) {
      pos = i;
      break;
    }

  return pos;

} /* find */

void compute_loc_matrix( double *x_triangle, double *y_triangle,
			 Epetra_SerialDenseMatrix & Ke ) 
{
  
  int    ii, jj;
  double det_J;
  double xa, ya, xb, yb, xc, yc;

  xa = x_triangle[0];
  xb = x_triangle[1];
  xc = x_triangle[2];
  ya = y_triangle[0];
  yb = y_triangle[1];
  yc = y_triangle[2];

  Ke(0,0) = (yc-yb)*(yc-yb) + (xc-xb)*(xc-xb);
  Ke(0,1) = (yc-yb)*(ya-yc) + (xc-xb)*(xa-xc);
  Ke(0,2) = (yb-ya)*(yc-yb) + (xb-xa)*(xc-xb);
  Ke(1,0) = (yc-yb)*(ya-yc) + (xc-xb)*(xa-xc);
  Ke(1,1) = (yc-ya)*(yc-ya) + (xc-xa)*(xc-xa);
  Ke(1,2) = (ya-yc)*(yb-ya) + (xa-xc)*(xb-xa);
  Ke(2,0) = (yb-ya)*(yc-yb) + (xb-xa)*(xc-xb);
  Ke(2,1) = (ya-yc)*(yb-ya) + (xa-xc)*(xb-xa);
  Ke(2,2) = (yb-ya)*(yb-ya) + (xb-xa)*(xb-xa);

  det_J = (xb-xa)*(yc-ya)-(xc-xa)*(yb-ya);
  det_J = 2*det_J;
  if( det_J<0 ) det_J *=-1;
  
  for (ii = 0; ii < 3; ii++) {
    for (jj = 0; jj < 3; jj++) {
      Ke(ii,jj) = Ke(ii,jj) / det_J;
    }
  }

  return;

} /* compute_loc_matrix */
