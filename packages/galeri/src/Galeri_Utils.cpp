// @HEADER
// ************************************************************************
//
//           Galeri: Finite Element and Matrix Generation Package
//                 Copyright (2006) ETHZ/Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions about Galeri? Contact Marzio Sala (marzio.sala _AT_ gmail.com)
//
// ************************************************************************
// @HEADER

#include "Galeri_Utils.h"
#include "Galeri_Exception.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Comm.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseSolver.h"
#include "Epetra_LAPACK.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_IntSerialDenseMatrix.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Teuchos_ParameterList.hpp"

namespace Galeri {

// ============================================================================
void
Solve(const Epetra_LinearProblem Problem)
{
  Solve(Problem.GetMatrix(), Problem.GetLHS(), Problem.GetRHS());
}

// ============================================================================
void
Solve(const Epetra_RowMatrix* Matrix, const Epetra_MultiVector* LHS,
      const Epetra_MultiVector* RHS)
{
  if (Matrix->Comm().NumProc() != 1)
    throw(Exception(__FILE__, __LINE__,
                    "Solve() works only in serial"));
  if (LHS->NumVectors() != RHS->NumVectors())
    throw(Exception(__FILE__, __LINE__,
                    "number of vectors in multivectors not consistent"));

  int n = Matrix->NumGlobalRows();
  int NumVectors = LHS->NumVectors();

  Epetra_SerialDenseMatrix DenseMatrix;
  DenseMatrix.Shape(n, n);

  for (int i = 0 ; i < n ; ++i)
    for (int j = 0 ; j < n ; ++j)
      DenseMatrix(i,j) = 0.0;

  // allocate storage to extract matrix rows.
  int Length = Matrix->MaxNumEntries();
  vector<double> Values(Length);
  vector<int>    Indices(Length);

  for (int j = 0 ; j < Matrix->NumMyRows() ; ++j)
  {
    int NumEntries;
    int ierr = Matrix->ExtractMyRowCopy(j, Length, NumEntries,
                                        &Values[0], &Indices[0]);

    for (int k = 0 ; k < NumEntries ; ++k)
      DenseMatrix(j,Indices[k]) = Values[k];
  }

  Epetra_SerialDenseMatrix DenseX(n, NumVectors);
  Epetra_SerialDenseMatrix DenseB(n, NumVectors);

  for (int i = 0 ; i < n ; ++i)
    for (int j = 0 ; j < NumVectors ; ++j)
      DenseB(i,j) = (*RHS)[j][i];

  Epetra_SerialDenseSolver DenseSolver;

  DenseSolver.SetMatrix(DenseMatrix);
  DenseSolver.SetVectors(DenseX,DenseB);

  DenseSolver.Factor();
  DenseSolver.Solve();

  for (int i = 0 ; i < n ; ++i)
    for (int j = 0 ; j < NumVectors ; ++j)
       (*LHS)[j][i] = DenseX(i,j);
}

// ============================================================================
double
ComputeNorm(const Epetra_MultiVector* LHS, const Epetra_MultiVector* RHS)
{
  double TotalNorm = 0.0;

  Epetra_MultiVector LHS2(*LHS);
  LHS2.Update(1.0, *RHS, -1.0);

  vector<double> norm(LHS->NumVectors());
  LHS2.Norm2(&norm[0]);

  for (int i = 0 ; i < LHS->NumVectors() ; ++i)
    TotalNorm += norm[i];

  return(TotalNorm);
}

// ============================================================================
double
ComputeNorm(const Epetra_RowMatrix* A, const Epetra_MultiVector* LHS,
            const Epetra_MultiVector* RHS)
{
  double TotalNorm = 0.0;

  Epetra_MultiVector Ax(*RHS);
  A->Multiply(false, *LHS, Ax);
  Ax.Update(1.0, *RHS, -1.0);

  vector<double> norm(LHS->NumVectors());
  Ax.Norm2(&norm[0]);

  for (int i = 0 ; i < LHS->NumVectors() ; ++i)
    TotalNorm += norm[i];

  return(TotalNorm);
}

// ============================================================================
Epetra_MultiVector*
CreateCartesianCoordinates(const string CoordType,
                           const Epetra_BlockMap* BlockMap,
                           Teuchos::ParameterList& List)
{
  // FIXME: pdes > 1
  double length = 1.0;
  double delta_x, delta_y, delta_z;

  double lx = List.get("lx", 1.0);
  double ly = List.get("ly", 1.0);
  double lz = List.get("lz", 1.0);

  int nx = List.get("nx", -1);
  int ny = List.get("ny", -1);
  int nz = List.get("nz", -1);

  int ix, iy, iz;

  int NumMyElements = BlockMap->NumMyElements();
  int * MyGlobalElements = BlockMap->MyGlobalElements();

  Epetra_MultiVector* Coord;

  if (CoordType == "1D")
  {
    Coord = new Epetra_MultiVector(*BlockMap, 1);

    delta_x = lx / (nx - 1);

    for (int i = 0 ; i < NumMyElements ; ++i)
    {
      ix = MyGlobalElements[i];
      (*Coord)[0][i] = delta_x * ix;
    }
  }
  else if (CoordType == "2D")
  {
    Coord = new Epetra_MultiVector(*BlockMap, 2);

    delta_x = lx / (nx - 1);
    delta_y = ly / (ny - 1);

    for (int i = 0 ; i < NumMyElements ; ++i)
    {
      ix = MyGlobalElements[i] % nx;
      iy = (MyGlobalElements[i] - ix) / nx;

      (*Coord)[0][i] = delta_x * ix;
      (*Coord)[1][i] = delta_y * iy;
    }
  }
  else if (CoordType == "3D")
  {
    Coord = new Epetra_MultiVector(*BlockMap, 3);

    delta_x = lx / (nx - 1);
    delta_y = ly / (ny - 1);
    delta_z = lz / (nz - 1);

    for (int i = 0 ; i < NumMyElements ; i++)
    {
      int ixy = MyGlobalElements[i] % (nx * ny);
      iz = (MyGlobalElements[i] - ixy) / (nx * ny);

      ix = ixy % nx;
      iy = (ixy - ix) / ny;

      (*Coord)[0][i] = delta_x * ix;
      (*Coord)[1][i] = delta_y * iy;
      (*Coord)[2][i] = delta_z * iz;
    }
  }
  else
  {
    throw(Exception(__FILE__, __LINE__,
                    "`CoordType' has incorrect value ("
                    + CoordType + ")",
                    "in input to function CreateCartesianCoordinates()",
                    "Check the documentation for a list of valid choices"));
  }

  return(Coord);
}

// ============================================================================
string toString(const int& x)
{
  char s[100];
  sprintf(s, "%d", x);
  return string(s);
}

// ============================================================================
string toString(const unsigned int& x)
{
  char s[100];
  sprintf(s, "%d", x);
  return string(s);
}

// ============================================================================
string toString(const long int& x)
{
  char s[100];
  sprintf(s, "%ld", x);
  return string(s);
}

// ============================================================================
string toString(const double& x)
{
  char s[100];
  sprintf(s, "%g", x);
  return string(s);
}

// ============================================================================
void GetNeighboursCartesian2d(const int i, const int nx, const int ny,
                              int & left, int & right,
                              int & lower, int & upper)
{
  int ix, iy;
  ix = i % nx;
  iy = (i - ix) / nx;

  if (ix == 0)      left = -1;
  else              left = i - 1;
  if (ix == nx - 1) right = -1;
  else              right = i + 1;
  if (iy == 0)      lower = -1;
  else              lower = i - nx;
  if (iy == ny - 1) upper = -1;
  else              upper = i + nx;
}

// ============================================================================
void GetNeighboursCartesian2d(const int i, const int nx, const int ny,
                              int& left, int& right, int& lower, int& upper,
                              int& left2, int& right2, int& lower2, int& upper2)
{
  int ix, iy;
  ix = i % nx;
  iy = (i - ix) / nx;

  if (ix == 0)      left = -1;
  else              left = i - 1;
  if (ix == nx - 1) right = -1;
  else              right = i + 1;
  if (iy == 0)      lower = -1;
  else              lower = i - nx;
  if (iy == ny - 1) upper = -1;
  else              upper = i + nx;

  if (ix <= 1)      left2 = -1;
  else              left2 = i - 2;
  if (ix >= nx - 2) right2 = -1;
  else              right2 = i + 2;
  if (iy <= 1)      lower2 = -1;
  else              lower2 = i - 2 * nx;
  if (iy >= ny - 2) upper2 = -1;
  else              upper2 = i + 2 * nx;
}

// ============================================================================
void GetNeighboursCartesian3d(const int i,
                              const int nx, const int ny, const int nz,
                              int& left, int& right, int& lower, int& upper,
                              int& below, int& above)
{
  int ixy, iz;
  ixy = i % (nx * ny);

  iz = (i - ixy) / (nx * ny);

  if (iz == 0)      below = -1;
  else              below = i - nx * ny;
  if (iz == nz - 1) above = -1;
  else              above = i + nx * ny;

  GetNeighboursCartesian2d(ixy, nx, ny, left, right, lower, upper);

  if (left  != -1) left  += iz * (nx * ny);
  if (right != -1) right += iz * (nx * ny);
  if (lower != -1) lower += iz * (nx * ny);
  if (upper != -1) upper += iz * (nx * ny);
}

// ============================================================================
void
PrintStencil2D(const Epetra_CrsMatrix* Matrix,
               const int nx, const int ny, int GID)
{
  if (nx <= 0 || ny <= 0)
      throw(Exception(__FILE__, __LINE__, "Input parameter not valid"));

  if (GID == -1)
  {
    if (ny == 1)
      GID = (int)(nx/2);
    else
      GID = (int)(nx*(ny/2) + nx/2);
  }

  int LID = Matrix->RowMatrixRowMap().LID(GID);

  // only processor having this node will go on
  if (LID == -1) return;

  int MaxPerRow = Matrix->MaxNumEntries();
  int NumEntriesRow;   // local entries on each row
  vector<double> Values(MaxPerRow);
  vector<int>    Indices(MaxPerRow);

  int ierr = Matrix->ExtractMyRowCopy(LID, MaxPerRow, NumEntriesRow,
                                      &Values[0], &Indices[0]);

  if (ierr)
    throw(Exception(__FILE__, __LINE__,
                    "Matrix->ExtractMyRowCopy() return an error"));

  // cycle over nonzero elements, look for elements in positions that we
  // can understand

  int size = 5;
  Epetra_IntSerialDenseMatrix SI(size, size);
  Epetra_SerialDenseMatrix    SV(size, size);

  for (int i = 0 ; i < size ; ++i)
    for (int j = 0 ; j < size ; ++j)
      SV(i, j) = 0.0;

  SI(0,0) = Matrix->RowMatrixColMap().LID(GID - 2 - 2 * nx);
  SI(1,0) = Matrix->RowMatrixColMap().LID(GID - 1 - 2 * nx);
  SI(2,0) = Matrix->RowMatrixColMap().LID(GID - 2 * nx);
  SI(3,0) = Matrix->RowMatrixColMap().LID(GID + 1 - 2 * nx);
  SI(4,0) = Matrix->RowMatrixColMap().LID(GID + 2 - 2 * nx);

  SI(0,1) = Matrix->RowMatrixColMap().LID(GID - 2 - nx);
  SI(1,1) = Matrix->RowMatrixColMap().LID(GID - 1 - nx);
  SI(2,1) = Matrix->RowMatrixColMap().LID(GID - nx);
  SI(3,1) = Matrix->RowMatrixColMap().LID(GID + 1 - nx);
  SI(4,1) = Matrix->RowMatrixColMap().LID(GID + 2 - nx);

  SI(0,2) = Matrix->RowMatrixColMap().LID(GID - 2);
  SI(1,2) = Matrix->RowMatrixColMap().LID(GID - 1);
  SI(2,2) = Matrix->RowMatrixColMap().LID(GID);
  SI(3,2) = Matrix->RowMatrixColMap().LID(GID + 1);
  SI(4,2) = Matrix->RowMatrixColMap().LID(GID + 2);

  SI(0,3) = Matrix->RowMatrixColMap().LID(GID - 2 + nx);
  SI(1,3) = Matrix->RowMatrixColMap().LID(GID - 1 + nx);
  SI(2,3) = Matrix->RowMatrixColMap().LID(GID - nx);
  SI(3,3) = Matrix->RowMatrixColMap().LID(GID + 1 + nx);
  SI(4,3) = Matrix->RowMatrixColMap().LID(GID + 2 + nx);

  SI(0,4) = Matrix->RowMatrixColMap().LID(GID - 2 + 2 * nx);
  SI(1,4) = Matrix->RowMatrixColMap().LID(GID - 1 + 2 * nx);
  SI(2,4) = Matrix->RowMatrixColMap().LID(GID - 2 * nx);
  SI(3,4) = Matrix->RowMatrixColMap().LID(GID + 1 + 2 * nx);
  SI(4,4) = Matrix->RowMatrixColMap().LID(GID + 2 + 2 * nx);

  for (int i = 0 ; i < NumEntriesRow ; ++i)
  {
    // convert into block row
    int LocalColID = Indices[i];
    // look for known positions
    for (int ix = 0 ; ix < size ; ++ix)
      for (int iy = 0 ; iy < size ; ++iy)
	if (SI(ix, iy) == LocalColID)
	  SV(ix,iy) = Values[i];
  }

  cout << "2D computational stencil at GID " << GID
       << " (grid is " << nx << " x " << ny << ")" << endl;
  cout << endl;
  for (int iy = 0 ; iy < size ; ++iy)
  {
    for (int ix = 0 ; ix < size ; ++ix)
    {
      cout << " " << std::setw(10) << SV(ix,iy);
    }
    cout << endl;
  }
  cout << endl;
}

} // namespace Galeri
