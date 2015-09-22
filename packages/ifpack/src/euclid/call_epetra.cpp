/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
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
