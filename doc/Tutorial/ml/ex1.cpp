
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
// Two-level domain decomposition preconditioner with AztecO and ML
//
// Marzio Sala, SNL, 9214, 19-Nov-2003

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
#include "AztecOO.h"
// includes required by ML
#include "ml_include.h"
#include "Epetra_LinearProblem.h"
#include "ml_epetra_operator.h"
#include "ml_epetra_utils.h"

void  get_neighbours( const int i, const int nx, const int ny,
		      int & left, int & right, 
		      int & lower, int & upper);

int main(int argc, char *argv[])
{

#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
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

  double off_left  = -1.0;
  double off_right = -1.0;
  double off_lower = -1.0;
  double off_upper = -1.0;
  double diag      =  4.0;
  int left, right, lower, upper;
  
  Epetra_CrsMatrix A(Copy,Map,5);

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

  // define vectors
  Epetra_Vector x(Map);
  Epetra_Vector b(Map);

  b.PutScalar(1.0);
  
  // define a linear problem
  Epetra_LinearProblem problem(&A, &x, &b);

  // Construct a solver object for this problem
  AztecOO solver(problem);

  // AZTEC settings
  solver.SetAztecOption(AZ_solver, AZ_cg);

  // Create and set an ML multilevel preconditioner
  ML *ml_handle;
  // Maximum number of levels
  int N_levels = 10;
  // output level
  ML_Set_PrintLevel(3);

  ML_Create(&ml_handle,N_levels);
  EpetraMatrix2MLMatrix(ml_handle, 0, &A);

  ML_Aggregate *agg_object;
  ML_Aggregate_Create(&agg_object);
  ML_Aggregate_Set_MaxCoarseSize(agg_object,1);
  N_levels = ML_Gen_MGHierarchy_UsingAggregation(ml_handle, 0,
                                                  ML_INCREASING, agg_object);
  // Set a symmetric Gauss-Seidel smoother for the MG method
  ML_Gen_Smoother_SymGaussSeidel(ml_handle, ML_ALL_LEVELS,
                                  ML_BOTH, 1, ML_DEFAULT);
  ML_Gen_Solver    (ml_handle, ML_MGV, 0, N_levels-1);

  Epetra_ML_Operator  MLop(ml_handle,Comm,Map,Map);

  solver.SetPrecOperator(&MLop);
 
  double rthresh = 1.4;
  solver.SetAztecParam(AZ_rthresh, rthresh);
  double athresh = 10.0;
  solver.SetAztecParam(AZ_athresh, athresh);
  solver.SetAztecParam(AZ_ill_cond_thresh, 1.0e200);

  int Niters = 500;
  solver.SetAztecOption(AZ_kspace, 160);
   
  solver.Iterate(Niters, 1e-12);

  // compute the real residual
  
  Epetra_Vector bcomp(Map);
  assert(A.Multiply(false, x, bcomp)==0);
 
  Epetra_Vector resid(Map);
 
  assert(resid.Update(1.0, b, -1.0, bcomp, 0.0)==0);

  double residual;
  assert(resid.Norm2(&residual)==0);
  if (Comm.MyPID()==0) cout << "Residual    = " << residual << endl;

  assert(resid.Update(1.0, x, -1.0, x, 0.0)==0);

  assert(resid.Norm2(&residual)==0);
  if (Comm.MyPID()==0)
    cout << "2-norm of difference between computed and exact solution  = " << residual << endl;

  if (residual>1.0e-5) {
    cout << "Difference between computed and exact solution is large..." << endl      << "Computing norm of A times this difference.  If this norm is small, then matrix is singular"
      << endl;
    assert(A.Multiply(false, resid, bcomp)==0);
    assert(bcomp.Norm2(&residual)==0);
  if (Comm.MyPID()==0)
    cout << "2-norm of A times difference between computed and exact solution  = " << residual << endl;
  
  }

#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

return 0 ;
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
