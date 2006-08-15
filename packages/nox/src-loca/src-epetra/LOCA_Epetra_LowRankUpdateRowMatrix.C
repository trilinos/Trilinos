// $Id$
// $Source$

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#include "Epetra_config.h"
#include "LOCA_Epetra_LowRankUpdateRowMatrix.H"

#include "Epetra_Vector.h"
#include "Epetra_Map.h"
#include "Epetra_Comm.h"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"

LOCA::Epetra::LowRankUpdateRowMatrix::
LowRankUpdateRowMatrix(
        const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
	const Teuchos::RefCountPtr<Epetra_RowMatrix>& jacRowMatrix, 
	const Teuchos::RefCountPtr<Epetra_MultiVector>& U_multiVec, 
	const Teuchos::RefCountPtr<Epetra_MultiVector>& V_multiVec,
	bool include_UV_terms) :
  LOCA::Epetra::LowRankUpdateOp(global_data, jacRowMatrix, U_multiVec,
				V_multiVec),
  J_rowMatrix(jacRowMatrix),
  nonconst_U(U_multiVec),
  nonconst_V(V_multiVec),
  includeUV(include_UV_terms),
  m(U_multiVec->NumVectors()),
  U_map(U_multiVec->Map()),
  V_map(V_multiVec->Map()),
  row_map(jacRowMatrix->RowMatrixRowMap())
{
}

LOCA::Epetra::LowRankUpdateRowMatrix::
~LowRankUpdateRowMatrix()
{
}

const Epetra_BlockMap &
LOCA::Epetra::LowRankUpdateRowMatrix::
Map() const
{
  return J_rowMatrix->Map();
}

int
LOCA::Epetra::LowRankUpdateRowMatrix::
NumMyRowEntries(int MyRow, int & NumEntries) const
{
  return J_rowMatrix->NumMyRowEntries(MyRow, NumEntries);
}

int
LOCA::Epetra::LowRankUpdateRowMatrix::
MaxNumEntries() const
{
  return J_rowMatrix->MaxNumEntries();
}

int
LOCA::Epetra::LowRankUpdateRowMatrix::
ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, double *Values, 
		 int * Indices) const
{
  // Get row of J
  int res = J_rowMatrix->ExtractMyRowCopy(MyRow, Length, NumEntries, Values, 
					  Indices);

  if (!includeUV)
    return res;

  // Add U*V^T to each entry of row
  int jac_row_gid = row_map.GID(MyRow);            // Global row of J
  int u_row_lid = U_map.LID(jac_row_gid);          // Local row of U
  for (int col=0; col<NumEntries; col++) {
    int jac_col_gid = Indices[col];                // Global col of J
    if (V_map.MyGID(jac_col_gid)) {
      int v_row_lid = V_map.LID(jac_col_gid);      // Local row of V
      Values[col] += computeUV(u_row_lid, v_row_lid);
    }
  }

  return res;
}

int
LOCA::Epetra::LowRankUpdateRowMatrix::
ExtractDiagonalCopy(Epetra_Vector & Diagonal) const
{
  // Get diagonal of J
  int res = J_rowMatrix->ExtractDiagonalCopy(Diagonal);

  if (!includeUV)
    return res;

  // Add U*V^T to each entry
  int numMyRows = J_rowMatrix->NumMyRows();
  for (int row=0; row<numMyRows; row++) {
    int jac_row_gid = row_map.GID(row);
    int u_row_lid = U_map.LID(jac_row_gid);
    int v_row_lid = V_map.LID(jac_row_gid);
    Diagonal[row] += computeUV(u_row_lid, v_row_lid);
  }

  return res;
}

int
LOCA::Epetra::LowRankUpdateRowMatrix::
Multiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  // Compute op(J)*X
  int res = J_rowMatrix->Multiply(TransA, X, Y);

  // Number of input vectors
  int numVecs = X.NumVectors();

  // Create temporary matrix to store V^T*X or U^T*X
  if (tmpMat == Teuchos::null || tmpMat->NumVectors() != numVecs) {
    tmpMat = Teuchos::rcp(new Epetra_MultiVector(localMap, numVecs, false));
  }

  if (!TransA) {

    // Compute V^T*X
    tmpMat->Multiply('T', 'N', 1.0, *V, X, 0.0);

    // Compute J*X + U*(V^T*X)
    Y.Multiply('N', 'N', 1.0, *U, *tmpMat, 1.0);

  }
  else {

    // Compute U^T*X
    tmpMat->Multiply('T', 'N', 1.0, *U, X, 0.0);

    // Compute J^T*X + V*(U^T*X)
    Y.Multiply('N', 'N', 1.0, *V, *tmpMat, 1.0);

  }

  return res;
}

int
LOCA::Epetra::LowRankUpdateRowMatrix::
Solve(bool Upper, bool Trans, bool UnitDiagonal, const Epetra_MultiVector& X, 
      Epetra_MultiVector& Y) const
{
  // We just ignore the U*V^T terms here
  return J_rowMatrix->Solve(Upper, Trans, UnitDiagonal, X, Y);
}

int
LOCA::Epetra::LowRankUpdateRowMatrix::
InvRowSums(Epetra_Vector& x) const
{
//   // Compute the inverse row-sums of J
//   int res = J_rowMatrix->InvRowSums(x);

//   if (!includeUV)
//     return res;

//   int numMyRows = J_rowMatrix->NumMyRows();
//   for (int jac_row=0; jac_row<numMyRows; jac_row++) {
//     int jac_row_gid = row_map.GID(jac_row);
//     int u_row_lid = U_map.LID(jac_row_gid);

//     int num_v_rows = V->MyLength();
//     double val = 0.0;
//     for (int v_row=0; v_row<num_v_rows; v_row++)
//       val += fabs(computeUV(u_row_lid, v_row));

//     double total_val = 0.0;
//     J_rowMatrix->Comm().SumAll(&val, &total_val, 1);

//     x[jac_row] = 1.0 / (1.0/x[jac_row] + total_val);
//   }

  int res;

  if (!includeUV)
    res = J_rowMatrix->InvRowSums(x);
  else {
    Epetra_Vector ones(V_map,false);
    ones.PutScalar(1.0);
    res = Multiply(true, ones, x);
    x.Reciprocal(x);
  }

  return res;
}

int
LOCA::Epetra::LowRankUpdateRowMatrix::
LeftScale(const Epetra_Vector& x)
{

  
  int res;
  if (J_rowMatrix->UseTranspose())
    res = J_rowMatrix->RightScale(x); // Right scale J^T
  else
    res = J_rowMatrix->LeftScale(x); // Left scale J

  // Now scale U
  for (int j=0; j<m; j++)
    (*nonconst_U)(j)->Multiply(1.0, x, *((*nonconst_U)(j)), 0.0); 

  return res;
}

int
LOCA::Epetra::LowRankUpdateRowMatrix::
InvColSums(Epetra_Vector& x) const
{
//   // Compute the inverse column-sums of J
//   int res = J_rowMatrix->InvColSums(x);

//   if (!includeUV)
//     return res;

//   int numMyCols = x.MyLength();
//   for (int col=0; col<numMyCols; col++) {
//     int col_gid = x.Map().GID(col);
//     int v_row_lid = V_map.LID(col_gid);

//     int num_u_rows = U->MyLength();
//     double val = 0.0;
//     for (int u_row=0; u_row<num_u_rows; u_row++)
//       val += fabs(computeUV(u_row, v_row_lid));

//     double total_val = 0.0;
//     J_rowMatrix->Comm().SumAll(&val, &total_val, 1);

//     x[col] = 1.0 / (1.0/x[col] + total_val);
//   }

  int res;

  if (!includeUV)
    res = J_rowMatrix->InvColSums(x);
  else {
    Epetra_Vector ones(V_map,false);
    ones.PutScalar(1.0);
    res = Multiply(false, ones, x);
    x.Reciprocal(x);
  }

  return res;
}

int
LOCA::Epetra::LowRankUpdateRowMatrix::
RightScale(const Epetra_Vector& x)
{

  int res;
  if (J_rowMatrix->UseTranspose())
    res = J_rowMatrix->LeftScale(x);  // Left scale J^T
  else
    res = J_rowMatrix->RightScale(x); // Right scale J

  // Now scale V
  for (int j=0; j<m; j++)
    (*nonconst_V)(j)->Multiply(1.0, x, *((*nonconst_V)(j)), 0.0); 

  return res;
}

bool
LOCA::Epetra::LowRankUpdateRowMatrix::
Filled() const
{
  return J_rowMatrix->Filled();
}

double
LOCA::Epetra::LowRankUpdateRowMatrix::
NormInf() const
{
  if (!includeUV)
    return J_rowMatrix->NormInf();

  Epetra_Vector tmp(row_map);
  InvRowSums(tmp);
  tmp.Reciprocal(tmp);

  double val;
  tmp.MaxValue(&val);

  return val;
}

double
LOCA::Epetra::LowRankUpdateRowMatrix::
NormOne() const
{
  if (!includeUV)
    return J_rowMatrix->NormOne();

  Epetra_Vector tmp(row_map);
  InvColSums(tmp);
  tmp.Reciprocal(tmp);

  double val;
  tmp.MaxValue(&val);

  return val;
}

int
LOCA::Epetra::LowRankUpdateRowMatrix::
NumGlobalNonzeros() const
{
  return J_rowMatrix->NumGlobalNonzeros();
}

int
LOCA::Epetra::LowRankUpdateRowMatrix::
NumGlobalRows() const
{
  return J_rowMatrix->NumGlobalRows();
}

int
LOCA::Epetra::LowRankUpdateRowMatrix::
NumGlobalCols() const
{
  return J_rowMatrix->NumGlobalCols();
}

int
LOCA::Epetra::LowRankUpdateRowMatrix::
NumGlobalDiagonals() const
{
  return J_rowMatrix->NumGlobalDiagonals();
}

int
LOCA::Epetra::LowRankUpdateRowMatrix::
NumMyNonzeros() const
{
  return J_rowMatrix->NumMyNonzeros();
}

int
LOCA::Epetra::LowRankUpdateRowMatrix::
NumMyRows() const
{
  return J_rowMatrix->NumMyRows();
}

int
LOCA::Epetra::LowRankUpdateRowMatrix::
NumMyCols() const
{
  return J_rowMatrix->NumMyCols();
}

int
LOCA::Epetra::LowRankUpdateRowMatrix::
NumMyDiagonals() const
{
  return J_rowMatrix->NumMyDiagonals();
}

bool
LOCA::Epetra::LowRankUpdateRowMatrix::
LowerTriangular() const
{
  return J_rowMatrix->LowerTriangular();
}

bool
LOCA::Epetra::LowRankUpdateRowMatrix::
UpperTriangular() const
{
  return J_rowMatrix->UpperTriangular();
}

const Epetra_Map &
LOCA::Epetra::LowRankUpdateRowMatrix::
RowMatrixRowMap() const
{
  return J_rowMatrix->RowMatrixRowMap();
}

const Epetra_Map &
LOCA::Epetra::LowRankUpdateRowMatrix::
RowMatrixColMap() const
{
  return J_rowMatrix->RowMatrixColMap();
}

const Epetra_Import *
LOCA::Epetra::LowRankUpdateRowMatrix::
RowMatrixImporter() const
{
  return J_rowMatrix->RowMatrixImporter();
}

double
LOCA::Epetra::LowRankUpdateRowMatrix::
computeUV(int u_row_lid, int v_row_lid) const
{
  double val = 0.0;
      
  // val = sum_{k=0}^m U(i,k)*V(j,k)
  for (int k=0; k<m; k++)
    val += (*U)[k][u_row_lid] * (*V)[k][v_row_lid];

  return val;
}
