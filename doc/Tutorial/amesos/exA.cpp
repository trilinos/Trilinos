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

/*
 * Solve a 2D Laplacian problem with AztecOO
 *
 * Marzio Sala, SNL 9214, 03-Oct-2003
 */

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
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Time.h"
#include "Amesos_Umfpack.h"
#include "Amesos_Parameter_List.h"
#include "Amesos_Superludist.h"


#define ML_SOLVE_WITH_AMESOS_UMFPACK      1
#define ML_SOLVE_WITH_AMESOS_SUPERLUDIST  2 
#define ML_SOLVE_WITH_AMESOS_MUMPS        3 
#define ML_SOLVE_WITH_AZTECOO             4

typedef struct ML_Amesos_Handle_Struct 
{
  void * AmesosProblem;
  void * Map;
  void * NewMap;
  void * Importer;
  void * Exporter;
} ML_Amesos_Handle;

int ML_Amesos_Gen(Epetra_CrsMatrix * Ematrix_orig,
		  ML_Amesos_Handle *Handle,
		  int choice, int DesiredProcs );

int ML_Amesos_Solve( ML_Amesos_Handle *Amesos_Handle, double ML_x[],
		     double ML_rhs[] );
void ML_Amesos_Destroy(ML_Amesos_Handle * Amesos_Handle);

void  get_neighbours( const int i, const int nx, const int ny,
		      int & left, int & right, 
		      int & lower, int & upper);

int main(int argc, char *argv[]) {

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  int nx = 5;
  int ny = 6;
  int NumGlobalElements = nx * ny;

  // create a map
  Epetra_Map Map(NumGlobalElements,0,Comm);
  // local number of rows
  int NumMyElements = Map.NumMyElements();
  // get update list
  int * MyGlobalElements = new int [NumMyElements];
  Map.MyGlobalElements( MyGlobalElements );

  // Create an integer vector NumNz that is used to build the Petra Matrix.
  // NumNz[i] is the Number of OFF-DIAGONAL term for the ith global equation 
  // on this processor

  int * NumNz = new int[NumMyElements];
  
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

  Epetra_CrsMatrix A(Copy,Map,NumNz);

  // Add  rows one-at-a-time
  // Need some vectors to help
  // Off diagonal Values will always be -1

  double *Values = new double[4];
  int *Indices = new int[4];
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
    assert(A.InsertGlobalValues(MyGlobalElements[i], NumEntries, 
				Values, Indices)==0);
    // Put in the diagonal entry
    assert(A.InsertGlobalValues(MyGlobalElements[i], 1, 
				&diag, MyGlobalElements+i)==0);
  }
  
  // Finish up
  assert(A.TransformToLocal()==0);

  // create x and b vectors
  Epetra_Vector x(Map);
  Epetra_Vector b(Map);
  x.PutScalar(1.0);

  A.Multiply(false,x,b);

  x.PutScalar(0.0);

  ML_Amesos_Handle Handle;

  ML_Amesos_Gen( &A, &Handle, ML_SOLVE_WITH_AMESOS_SUPERLUDIST, 3);

  ML_Amesos_Solve( &Handle, x.Values(), b.Values());

  b.PutScalar(1.0);
  x.Update(-1.0,b,1.0);

  double norm2;
  x.Norm2(&norm2);

  ML_Amesos_Destroy( &Handle );
  
  if( Comm.MyPID() == 0 ) 
    cout << "|| x - x_exact ||_2 = " << norm2 << endl;
  
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










/* ======================================================================== */
/*!
 \brief find out the number of local rows to be assigned to this process
 so that only the first DesiredProcs are used

*/
/* ------------------------------------------------------------------------ */

int ML_FindLocalRows_GivenNumberOfProcs( int NumLocalRows,
					 int DesiredProcs,
					 Epetra_RowMatrix * matrix )
{

  int NewNumLocalRows = 0;
  int MyPID = matrix->Comm().MyPID();
  
  // I suppose that no rows are replicated among the processes
  // (as usually done within ML, Epetra is more flexible)
  int NumGlobalRows;
  matrix->Comm().SumAll( &NumLocalRows, &NumGlobalRows, 1 );

  int base = NumGlobalRows/DesiredProcs;
  int mod = NumGlobalRows%DesiredProcs;
  int modmod = mod%DesiredProcs;
  int moddiv = mod/DesiredProcs;

  if( MyPID < DesiredProcs ) {
    NewNumLocalRows = base + moddiv;
  }

  if( MyPID == 0 ) NewNumLocalRows + modmod;

  return NewNumLocalRows;
  
} /* ML_FindLocalRows_GivenNumberOfProcs */


/* ======================================================================== */
/*!
 \brief Setup to use various Amesos solver on an ML_Operator

 This function performs the following actions:
 - convert the ML_Operator of current level to an Epetra_CrsMatrix
 - if required, redistribute this matrix among a subset of the number of
   processors, and creates the Epetra_Import and Epetra_Export objects
   required to redistribute solution and rhs
 - perform symbolic and numerical factorization of the Epetra_CrsMatrix

*/
/* ------------------------------------------------------------------------ */

int ML_Amesos_Gen(Epetra_CrsMatrix * Ematrix_orig,
		  ML_Amesos_Handle *Handle,
		  int choice, int DesiredProcs )
{

  Epetra_CrsMatrix * Ematrix;
  
  // I suppose the problem is square (DomainMap is equal to RangeMap)
  Handle->Map = (void*)&(Ematrix_orig->OperatorDomainMap());
  
  // check input
  int NumProcs = Ematrix_orig->Comm().NumProc();
  if( DesiredProcs > NumProcs) {
    cerr << "*ML*WRN* DesiredProcs > NumProcs (" << DesiredProcs;
    cerr << " - " << NumProcs << ")\n";
    cerr << "*ML*WRN* continuing with DesiredProcs = NumProcs\n";
    cerr << "*ML*WRN* (file " << __FILE__ << ",  line " << __LINE__ <<")\n";
    DesiredProcs = -1;
  }
  
  if( DesiredProcs > 0 && DesiredProcs != NumProcs ) {

    // redistribute the matrix

    int NumLocalRows = Ematrix_orig->NumMyRows();
    int NumGlobalRows = Ematrix_orig->NumGlobalRows();
    int NewNumLocalRows = ML_FindLocalRows_GivenNumberOfProcs(NumLocalRows,
							      DesiredProcs,
							      Ematrix_orig);

    // vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    //    if( Ematrix_orig->Comm().MyPID() == 0 ) NewNumLocalRows = 30;
    //    else NewNumLocalRows = 0;
    // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    
    Epetra_Map * NewMap = new Epetra_Map(NumGlobalRows, NewNumLocalRows, 0,
					 Ematrix_orig->Comm());
    // NOTE: both Importer and Exporter are Epetra_Import. Exporter means
    // to export what computed on NewMap to the "old" Map.
    Epetra_Import * Importer = new Epetra_Import(*NewMap, Ematrix_orig->Map());
    Epetra_Import * Exporter = new Epetra_Import(Ematrix_orig->Map(), *NewMap);
    
    Ematrix = new Epetra_CrsMatrix(Copy,*NewMap,0);

    Ematrix->Import(*Ematrix_orig,*Importer,Add);
    Ematrix->FillComplete();
    
    Handle->Importer = (void *)Importer;
    Handle->Exporter = (void *)Exporter;
    Handle->NewMap = NewMap;

  } else {

    // do nothing in this case
    Ematrix = Ematrix_orig;
    Handle->Importer = NULL;
    Handle->Exporter = NULL;
    Handle->NewMap = NULL;
    
  }
  
  // create the linear problem with the new distribution
  
  Epetra_LinearProblem * LinearProblem = new Epetra_LinearProblem;
  LinearProblem->SetOperator( Ematrix ); 

  // parameter vector for Amesos
  AMESOS::Parameter::List ParamList;

  // prepare the linear solver

  Amesos_BaseSolver * AmesosProblem;

  switch( choice ) {
  case ML_SOLVE_WITH_AMESOS_UMFPACK:
    AmesosProblem = new Amesos_Umfpack( *LinearProblem, ParamList );
    break;
  case ML_SOLVE_WITH_AMESOS_SUPERLUDIST:
    AmesosProblem = new Amesos_Superludist( *LinearProblem, ParamList );
    break;
  case ML_SOLVE_WITH_AMESOS_MUMPS:
    //    AmesosProblem = new Amesos_Mumps( *LinearProblem, ParamList );
    break;
  default:
    cerr << "*ML*ERR* wrong parameter specified (" << choice << ")\n";
    cerr << "*ML*ERR* (file " << __FILE__ << ", line " << __LINE__ << ")\n";
    exit( EXIT_FAILURE );
  }
  
  // solve and factorize  

  AmesosProblem->SymbolicFactorization();

  AmesosProblem->NumericFactorization();

  Handle->AmesosProblem = (void *) AmesosProblem ;
  
  return 0;

} /* ML_Amesos_Gen */

/* ======================================================================== */
/*!
 \brief Solve a linear problem using a handle created by `ML_Amesos_Gen'.

 \note this function makes use of the Amesos_BaseSolver class.
*/
/* ------------------------------------------------------------------------ */

int ML_Amesos_Solve( ML_Amesos_Handle *Amesos_Handle, double ML_x[],
		     double ML_rhs[] )
{
  
  Amesos_BaseSolver * AmesosProblem =
    (Amesos_BaseSolver *) (Amesos_Handle->AmesosProblem) ;
  Epetra_LinearProblem *Amesos_LinearProblem =
    (Epetra_LinearProblem *) AmesosProblem->GetProblem() ; 

  Epetra_Map * Map = (Epetra_Map*)Amesos_Handle->Map;
  Epetra_Map * NewMap = (Epetra_Map*)Amesos_Handle->NewMap;

  Epetra_Import * Importer =
     (Epetra_Import *) (Amesos_Handle->Importer) ;
  Epetra_Import * Exporter =
     (Epetra_Import *) (Amesos_Handle->Exporter) ;

  
  // here we should put some ReserView to avoid allocating vectors all times
  
  Epetra_Vector rhs( Copy, *Map, ML_rhs ) ;
  Epetra_Vector lhs( View, *Map, ML_x ) ;
  
  if( Importer == NULL ) {

    // matrix has not been redistributed, rhs and lhs have the same
    // map of matrix
    Amesos_LinearProblem->SetRHS( &rhs ) ; 
    Amesos_LinearProblem->SetLHS( &lhs ) ;

    AmesosProblem->Solve() ; 
    
  } else {

    // matrix has been redistributed. I need to Import the new rhs,
    // then to export the computer solution

    Epetra_Vector new_rhs( *NewMap ) ;
    Epetra_Vector new_lhs( *NewMap ) ;
    new_rhs.Import(rhs,*Importer,Insert);
    new_lhs.PutScalar(0.0);
    
    Amesos_LinearProblem->SetRHS( &new_rhs ) ; 
    Amesos_LinearProblem->SetLHS( &new_lhs ) ;

    // solve with Amesos on the new Map
    AmesosProblem->Solve() ; 

    // redistribute results according to the original map
    lhs.Import(new_lhs,*Exporter,Insert);

  }
  
  return 0;
  
} /* ML_Amesos_Solve */

/* ======================================================================== */
/*!
 \brief destroy the objects created by `ML_Amesos_Gen'.

*/
/* ------------------------------------------------------------------------ */

void ML_Amesos_Destroy(ML_Amesos_Handle * Amesos_Handle)
{

  if( Amesos_Handle->Importer != NULL ) {
    Epetra_Import * Importer =
      (Epetra_Import *) (Amesos_Handle->Importer) ;
    delete Importer;
  }

  if( Amesos_Handle->Exporter != NULL ) {
    Epetra_Import * Exporter =
      (Epetra_Import *) (Amesos_Handle->Exporter) ;
    delete Exporter;
  }
  
  Amesos_BaseSolver * AmesosProblem =
    (Amesos_BaseSolver *) (Amesos_Handle->AmesosProblem);

  Epetra_LinearProblem * Amesos_LinearProblem =
    (Epetra_LinearProblem *) AmesosProblem->GetProblem();

 // delete Amesos_LinearProblem->GetOperator(); 
  delete Amesos_LinearProblem;

  delete AmesosProblem;

//  Epetra_Map * Map = (Epetra_Map*)Amesos_Handle->Map;
//  delete Map;
  
  if( Amesos_Handle->NewMap != NULL ) {
    Epetra_Map * NewMap = (Epetra_Map*)Amesos_Handle->NewMap;
    delete NewMap;
  }
  
} /* ML_Amesos_Destroy */

