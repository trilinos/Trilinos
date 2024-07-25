// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "NOX_Epetra_Interface_Required.H"
#include "NOX_Epetra_Interface_Jacobian.H"

#include "Laplace2D.H"

// this is required to know the number of lower, upper, left and right
// node for each node of the Cartesian grid (composed by nx \timex ny
// elements)

void
Laplace2D::get_myNeighbours( const int i, const int nx, const int ny,
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

// This function creates a CrsMatrix, whose elements corresponds
// to the discretization of a Laplacian over a Cartesian grid,
// with nx grid point along the x-axis and and ny grid points
// along the y-axis. For the sake of simplicity, I suppose that
// all the nodes in the matrix are internal nodes (Dirichlet
// boundary nodes are supposed to have been already condensated)

Epetra_CrsMatrix *
Laplace2D::CreateLaplacian( const int nx, const int ny, const Epetra_Comm * Comm)
{

  int NumGlobalElements = nx * ny;

  // create a map
  Epetra_Map * Map = new Epetra_Map(NumGlobalElements,0,*Comm);
  // local number of rows
  int NumMyElements = Map->NumMyElements();
  // get update list
  int * MyGlobalElements = Map->MyGlobalElements();

  double hx = 1.0/(nx-1);
  double hy = 1.0/(ny-1);
  double off_left  = -1.0/(hx*hx);
  double off_right = -1.0/(hx*hx);
  double off_lower = -1.0/(hy*hy);
  double off_upper = -1.0/(hy*hy);
  double diag      =  2.0/(hx*hx) + 2.0/(hy*hy);

  int left, right, lower, upper;

  // a bit overestimated the nonzero per row

  Epetra_CrsMatrix * A = new Epetra_CrsMatrix(Copy,*Map,5);

  // Add  rows one-at-a-time

  double * Values = new double[4];
  int *   Indices = new int[4];

  for( int i = 0; i < NumMyElements; ++i )
  {
    int NumEntries=0;
    get_myNeighbours(  MyGlobalElements[i], nx, ny, left, right, lower, upper );
    if( left != -1 )
    {
      Indices[NumEntries] = left;
      Values[NumEntries] = off_left;
      ++NumEntries;
    }
    if( right != -1 )
    {
      Indices[NumEntries] = right;
      Values[NumEntries] = off_right;
      ++NumEntries;
    }
    if( lower != -1 )
    {
      Indices[NumEntries] = lower;
      Values[NumEntries] = off_lower;
      ++NumEntries;
    }
    if( upper != -1 )
    {
      Indices[NumEntries] = upper;
      Values[NumEntries] = off_upper;
      ++NumEntries;
    }
    // put the off-diagonal entries
    A->InsertGlobalValues(MyGlobalElements[i], NumEntries, Values, Indices);
    // Put in the diagonal entry
    A->InsertGlobalValues(MyGlobalElements[i], 1, &diag, MyGlobalElements+i);
  }

  // put matrix in local ordering
  A->FillComplete();

  delete [] Indices;
  delete [] Values;
  delete    Map;

  return A;

} /* createJacobian */

// ==========================================================================
// This class contians the main definition of the nonlinear problem at
// hand. A method is provided to compute F(x) for a given x, and another
// method to update the entries of the Jacobian matrix, for a given x.
// As the Jacobian matrix J can be written as
//    J = L + diag(lambda*exp(x[i])),
// where L corresponds to the discretization of a Laplacian, and diag
// is a diagonal matrix with lambda*exp(x[i]). Basically, to update
// the jacobian we simply update the diagonal entries. Similarly, to compute
// F(x), we reset J to be equal to L, then we multiply it by the
// (distributed) vector x, then we add the diagonal contribution
// ==========================================================================

// constructor. Requires the number of nodes along the x-axis
// and y-axis, the value of lambda, and the Epetra_Communicator
// (to define a Map, which is a linear map in this case)
PDEProblem::PDEProblem(const int nx, const int ny, const double lambda,
           const Epetra_Comm * Comm) :
  nx_(nx), ny_(ny), lambda_(lambda)
{
  hx_ = 1.0/(nx_-1);
  hy_ = 1.0/(ny_-1);
  Matrix_ = Laplace2D::CreateLaplacian(nx_,ny_,Comm);
}

// destructor
PDEProblem::~PDEProblem()
{
  delete Matrix_;
}

// compute F(x)
void PDEProblem::ComputeF( const Epetra_Vector & x, Epetra_Vector & f )
{
  // reset diagonal entries
  double diag =  2.0/(hx_*hx_) + 2.0/(hy_*hy_);

  int NumMyElements = Matrix_->Map().NumMyElements();
  // get update list
  int * MyGlobalElements = Matrix_->Map().MyGlobalElements();

  for( int i = 0; i < NumMyElements; ++i )
  {
    // Put in the diagonal entry
    Matrix_->ReplaceGlobalValues(MyGlobalElements[i], 1, &diag, MyGlobalElements+i);
  }
  // matrix-vector product (intra-processes communication occurs in this call)
  Matrix_->Multiply( false, x, f );

  // add diagonal contributions
  for( int i = 0; i < NumMyElements; ++i )
  {
    // Put in the diagonal entry
    f[i] += lambda_*exp(x[i]);
  }
}

// update the Jacobian matrix for a given x
void PDEProblem::UpdateJacobian( const Epetra_Vector & x )
{
  double diag =  2.0/(hx_*hx_) + 2.0/(hy_*hy_);

  int NumMyElements = Matrix_->Map().NumMyElements();
  // get update list
  int * MyGlobalElements = Matrix_->Map().MyGlobalElements();

  for( int i = 0; i < NumMyElements; ++i )
  {
    // Put in the diagonal entry
    double newdiag = diag + lambda_*exp(x[i]);
    Matrix_->ReplaceGlobalValues(MyGlobalElements[i], 1, &newdiag, MyGlobalElements+i);
  }

}


