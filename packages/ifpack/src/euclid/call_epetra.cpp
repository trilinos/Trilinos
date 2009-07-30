/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
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
//@HEADER
*/

#include "call_epetra.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include <vector>


int ExtractIndicesView(void* A, int GlobalRow, int *NumEntries, int **Indices){
  int ierr = 0;
  Epetra_CrsMatrix *mat = (Epetra_CrsMatrix*) A;
  int MyRow = mat->LRID(GlobalRow);
  if(MyRow == -1) return 0;
  int & len = *NumEntries;
  int *& ind = *Indices;
  ierr += mat->Graph().ExtractMyRowView(MyRow, len, ind);
  return ierr;
} 

int ExtractValuesView(void *A, int GlobalRow, int *NumEntries, double** Values){
  int ierr = 0;
  Epetra_CrsMatrix *mat = (Epetra_CrsMatrix*) A;
  int MyRow = mat->LRID(GlobalRow);
  if(MyRow == -1) return 0;
  int &len = *NumEntries;
  double *& val = *Values;
  ierr += mat->ExtractMyRowView(MyRow, len, val);
  return ierr;
}

int MinMaxMyGID(void* A, bool Row, bool min){
  Epetra_CrsMatrix* mat = (Epetra_CrsMatrix*) A;
  if(Row){
    if(min){
      return mat->RowMap().MinMyGID();
    } else { 
      return mat->RowMap().MaxMyGID();
    }
  } else {
    if(min){
      return mat->ColMap().MinMyGID();
    } else {
      return mat->ColMap().MaxMyGID();
    }
  }
}

int NumGlobalRowCol(void* A, bool Row){
  Epetra_CrsMatrix *mat = (Epetra_CrsMatrix*) A;
  if(Row){
    return mat->NumGlobalRows();
  } else {
    return mat->NumGlobalCols();
  }
}

int NumMyRowEntries(void *A, int Row, int *numEntries){
  Epetra_CrsMatrix* mat = (Epetra_CrsMatrix*) A;
  return mat->NumMyRowEntries(mat->LRID(Row), *numEntries);
}
