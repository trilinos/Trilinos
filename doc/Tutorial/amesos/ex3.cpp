
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
// Solve a linear system with SuperLU
//
// Marzio Sala, SNL, 9214, 20-Nov-2003

#include <iostream>
#include <fstream>

#include "Epetra_config.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_IntVector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Time.h"
#include "Amesos_Umfpack.h"
#include "Amesos_Parameter_List.h"
#include "Amesos_Superludist.h"

void  get_neighbours( const int i, const int nx, const int ny,
		      int & left, int & right, 
		      int & lower, int & upper);

/* for each of gallery functions, I must provide:
   - this one
   - one without exact solution
   - one without vectors
   - one with Map in input
   - one with default values
*/

#include <map>
#include <string>

typedef map<string,string> Triutils_MatrixGalleryList;

void Triutils_GetFromList( Triutils_MatrixGalleryList list,
			   string label, int & value,
			   bool & found ) 
{

  bool f = false;
  
  for( map<string,string>::const_iterator ci = list.begin() ;
       ci != list.end() ; ++ci ) {
    if( ci->first == label ) {
      value = atoi(list[label].c_str());
      found = true;
      return;
    }
  }

  return;

}

void Triutils_GetFromList( Triutils_MatrixGalleryList list,
			   string label, double & value,
			   bool & found ) 
{

  bool f = false;
  
  for( map<string,string>::const_iterator ci = list.begin() ;
       ci != list.end() ; ++ci ) {
    if( ci->first == label ) {
      value = atof(list[label].c_str());
      found = true;
      return;
    }
  }

  return;

}

  
    
int Triutils_MatrixGallery(const char *name,
			   Triutils_MatrixGalleryList list,
			   Epetra_CrsMatrix * & A,
			   Epetra_MultiVector * & Xex,
			   Epetra_MultiVector * & X,
			   Epetra_MultiVector * & B,
			   Epetra_Map * & Map,
			   const Epetra_Comm & Comm)
  
{

  // variables used by all matrices
  int NumVectors;

  // default values
  NumVectors = 0;

  int * MyGlobalElements = NULL;
  double * Values = NULL;
  int * Indices = NULL;
  int * NumNz = NULL;

  bool found;  
  
  Triutils_GetFromList( list, "NumVectors", NumVectors, found );
  if( found == false ) NumVectors = 0;
  
  if( strcmp(name,"Laplacian2D") == 0 ) {

    int nx, ny;

    Triutils_GetFromList( list, "nx", nx, found );
    if( found == false ) nx = 5;
    Triutils_GetFromList( list, "ny", ny, found );
    if( found == false ) ny = 6;
    
    int NumGlobalElements = nx * ny;
    
    // create a map
    Map = new Epetra_Map(NumGlobalElements,0,Comm);
    
    // local number of rows
    int NumMyElements = Map->NumMyElements();
    // get update list
    MyGlobalElements = new int [NumMyElements];
    Map->MyGlobalElements( MyGlobalElements );

    // Create an integer vector NumNz that is used to build the Petra Matrix.
    // NumNz[i] is the Number of OFF-DIAGONAL term for the ith global equation 
    // on this processor

    NumNz = new int[NumMyElements];
  
    double off_left  = -1.0;
    double off_right = -1.0;
    double off_lower = -1.0;
    double off_upper = -1.0;
    double diag      =  4.0;
    int left, right, lower, upper;
  
    for ( int i=0; i<NumMyElements; i++) {
      NumNz[i] = 1;
      get_neighbours( MyGlobalElements[i], nx, ny, 
		      left, right, lower, upper); 
      if( left  != -1 ) ++NumNz[i];
      if( right != -1 ) ++NumNz[i];
      if( lower != -1 ) ++NumNz[i];
      if( upper != -1 ) ++NumNz[i];
    }
  
    // Create a Epetra_Matrix
    // create a CRS matrix
    
    A = new Epetra_CrsMatrix(Copy,*Map,NumNz);
    
    // Add  rows one-at-a-time
    // Need some vectors to help
    // Off diagonal Values will always be -1
    
    Values = new double[4];
    Indices = new int[4];
    int NumEntries;
    
    for( int i=0 ; i<NumMyElements; ++i ) {
      int NumEntries=0;
      get_neighbours(  MyGlobalElements[i], nx, ny, 
		       left, right, lower, upper);
      if( left != -1 ) {
	Indices[NumEntries] = left;
	Values[NumEntries] = off_left;
	++NumEntries;
      }
      if( right != -1 ) {
	Indices[NumEntries] = right;
	Values[NumEntries] = off_right;
	++NumEntries;
      }
      if( lower != -1 ) {
	Indices[NumEntries] = lower;
	Values[NumEntries] = off_lower;
	++NumEntries;
      }
      if( upper != -1 ) {
	Indices[NumEntries] = upper;
	Values[NumEntries] = off_upper;
	++NumEntries;
      }
      // put the off-diagonal entries
      assert(A->InsertGlobalValues(MyGlobalElements[i], NumEntries, 
				   Values, Indices)==0);
      // Put in the diagonal entry
      assert(A->InsertGlobalValues(MyGlobalElements[i], 1, 
				   &diag, MyGlobalElements+i)==0);
    }
    
  } else if( strcmp(name,"Laplacian1D") == 0 ) {

    int nx, NumGlobalElements;
    
    Triutils_GetFromList( list, "nx", NumGlobalElements, found );
    if( found == false ) NumGlobalElements = 5;    

    // create a map
    Map = new Epetra_Map(NumGlobalElements,0,Comm);
    
    // local number of rows
    int NumMyElements = Map->NumMyElements();
    // get update list
    MyGlobalElements = new int [NumMyElements];
    Map->MyGlobalElements( MyGlobalElements );

    A = new Epetra_CrsMatrix(Copy,*Map,3);
    
    // Add  rows one-at-a-time
    // Need some vectors to help
    // Off diagonal Values will always be -1
    
    Values = new double[2];
    Indices = new int[2];
    int NumEntries;
    double diag = 2.0;
    
    for( int i=0 ; i<NumMyElements; ++i ) {
      int NumEntries=0;
      int GID = MyGlobalElements[i];

      if( GID != 0 ) {
	Indices[NumEntries] = GID-1;
	Values[NumEntries] = -1.0;
	NumEntries++;
      }

      if( GID != NumGlobalElements-1 ) {
	Indices[NumEntries] = GID+1;
	Values[NumEntries] = -1.0;
	NumEntries++;
      }
      
      // put the off-diagonal entries
      assert(A->InsertGlobalValues(MyGlobalElements[i], NumEntries, 
				   Values, Indices)==0);
      // Put in the diagonal entry
      assert(A->InsertGlobalValues(MyGlobalElements[i], 1, 
				   &diag, MyGlobalElements+i)==0);
    }
    
  } else if( strcmp(name,"diag") == 0 ) {

    int nx, NumGlobalElements;
    double diag_value;
    
    Triutils_GetFromList( list, "nx", NumGlobalElements, found );
    if( found == false ) NumGlobalElements = 5;    

    Triutils_GetFromList( list, "diag value", diag_value, found );
    if( found == false ) diag_value = 1.0;    

    // create a map
    Map = new Epetra_Map(NumGlobalElements,0,Comm);
    
    // local number of rows
    int NumMyElements = Map->NumMyElements();
    // get update list
    MyGlobalElements = new int [NumMyElements];
    Map->MyGlobalElements( MyGlobalElements );

    A = new Epetra_CrsMatrix(Copy,*Map,1);
    
    for( int i=0 ; i<NumMyElements; ++i ) {
      // Put in the diagonal entry
      assert(A->InsertGlobalValues(MyGlobalElements[i], 1, 
				   &diag_value, MyGlobalElements+i)==0);
    }
    
  }  
  
  if( MyGlobalElements != NULL ) delete MyGlobalElements;
  if( NumNz != NULL ) delete NumNz;
  if( Indices != NULL ) delete Indices;
  if( Values != NULL ) delete Values;
  
  // Finish up
  assert(A->TransformToLocal()==0);
  
  if( NumVectors == 0 ) return 0;
  
  // create x and b vectors
  X = new Epetra_MultiVector(A->DomainMap(), NumVectors);
  Xex = new Epetra_MultiVector(A->DomainMap(), NumVectors);
  B = new Epetra_MultiVector(A->RangeMap(), NumVectors ); 

  Xex->Random();
    
  A->Multiply(false,*Xex,*B);
    
  X->PutScalar(0.0);
    
  return 0;

};


	  
int Triutils_MatrixGallery(const char *name,
			   Triutils_MatrixGalleryList list,
			   Epetra_CrsMatrix * & A,
			   Epetra_Map * & Map,
			   const Epetra_Comm & Comm)
{

  int oldNumVectors = -1;
  
  for( map<string,string>::const_iterator ci = list.begin() ;
       ci != list.end() ; ++ci ) {
    oldNumVectors = atoi(list["NumVectors"].c_str());
  }      
  list["NumVectors"] = 1;

  Epetra_MultiVector *X, *B, *Xex;
  
  Triutils_MatrixGallery(name, list, A, Xex, X, B, Map, Comm);

  if( oldNumVectors != -1 ) {
    list["NumVectors"] = oldNumVectors;
  }

  return 0;
  
}

int main(int argc, char *argv[]) {

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  Triutils_MatrixGalleryList list;

  list["NumVectors"] = "1";
  list["nx"] = "10";
  list["ny"] = "10";
  
  Epetra_CrsMatrix * A;
  Epetra_MultiVector * X;
  Epetra_MultiVector * Xex;
  Epetra_MultiVector * B;
  Epetra_Map * Map;
  
  int NumVectors = 1;
    
  Triutils_MatrixGallery("Laplacian1D", list, A, Xex, X, B, Map, Comm );
  
  
  // create linear problem
  Epetra_LinearProblem Problem(A,X,B);

  AMESOS::Parameter::List params;
  
  Amesos_Superludist * SuperludistProblem =
    new Amesos_Superludist(Problem,params);

  SuperludistProblem->SymbolicFactorization();
  SuperludistProblem->NumericFactorization();
  SuperludistProblem->Solve();
  
  X->Update(-1.0,*Xex, 1.0);

  double error;
  
  X->Norm2(&error);

  if( Comm.MyPID() == 0 ) cout << "|| Xex - X ||_2 = " << error << endl;
  
  delete SuperludistProblem;
  
#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return( EXIT_SUCCESS );

}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void  get_neighbours( const int i, const int nx, const int ny,
		      int & left, int & right, 
		      int & lower, int & upper) 
{

  int ix, iy;
  ix = i%nx;
  iy = (i - ix)/nx;

  if( ix == 0 ) 
    left = -1;
  else 
    left = i-1;
  if( ix == nx-1 ) 
    right = -1;
  else
    right = i+1;
  if( iy == 0 ) 
    lower = -1;
  else
    lower = i-nx;
  if( iy == ny-1 ) 
    upper = -1;
  else
    upper = i+nx;

  return;

}
