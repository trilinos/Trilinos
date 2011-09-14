//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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

/*#############################################################################
# CVS File Information
#    Current revision: $Revision$
#    Last modified:    $Date$
#    Modified by:      $Author$
#############################################################################*/

#include "EpetraExt_ConfigDefs.h"
#ifdef HAVE_PETSC

#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "EpetraExt_PETScAIJMatrix.h"

//==============================================================================
Epetra_PETScAIJMatrix::Epetra_PETScAIJMatrix(Mat Amat)
  : Epetra_Object("Epetra::PETScAIJMatrix"),
    Amat_(Amat),
    Values_(0),
    Indices_(0),
    MaxNumEntries_(-1),
    ImportVector_(0),
    NormInf_(-1.0),
    NormOne_(-1.0)
{
#ifdef HAVE_MPI
  MPI_Comm comm;
  PetscObjectGetComm( (PetscObject)Amat, &comm);
  Comm_ = new Epetra_MpiComm(comm);
#else
  Comm_ = new Epetra_SerialComm();
#endif  
  int ierr;
  char errMsg[80];
  MatGetType(Amat, &MatType_);
  if ( strcmp(MatType_,MATSEQAIJ) != 0 && strcmp(MatType_,MATMPIAIJ) != 0 ) {
    sprintf(errMsg,"PETSc matrix must be either seqaij or mpiaij (but it is %s)",MatType_);
    throw Comm_->ReportError(errMsg,-1);
  }
  petscMatrixType mt;
  Mat_MPIAIJ* aij=0;
  if (strcmp(MatType_,MATMPIAIJ) == 0) {
    mt = PETSC_MPI_AIJ;
    aij = (Mat_MPIAIJ*)Amat->data;
  }
  else if (strcmp(MatType_,MATSEQAIJ) == 0) {
    mt = PETSC_SEQ_AIJ;
  }
  int numLocalRows, numLocalCols;
  ierr = MatGetLocalSize(Amat,&numLocalRows,&numLocalCols);
  if (ierr) {
    sprintf(errMsg,"EpetraExt_PETScAIJMatrix.cpp, line %d, MatGetLocalSize() returned error code %d",__LINE__,ierr);
    throw Comm_->ReportError(errMsg,-1);
  }
  NumMyRows_ = numLocalRows;
  NumMyCols_ = numLocalCols; //numLocalCols is the total # of unique columns in the local matrix (the diagonal block)
  //TODO what happens if some columns are empty?
  if (mt == PETSC_MPI_AIJ)
    NumMyCols_ += aij->B->cmap->n;
  MatInfo info;
  ierr = MatGetInfo(Amat,MAT_LOCAL,&info);
  if (ierr) {
    sprintf(errMsg,"EpetraExt_PETScAIJMatrix.cpp, line %d, MatGetInfo() returned error code %d",__LINE__,ierr);
    throw Comm_->ReportError(errMsg,-1);
  }
  NumMyNonzeros_ = (int) info.nz_used; //PETSc stores nnz as double
  Comm_->SumAll(&(info.nz_used), &NumGlobalNonzeros_, 1);

  //The PETSc documentation warns that this may not be robust.
  //In particular, this will break if the ordering is not contiguous!
  int rowStart, rowEnd;
  ierr = MatGetOwnershipRange(Amat,&rowStart,&rowEnd);
  if (ierr) {
    sprintf(errMsg,"EpetraExt_PETScAIJMatrix.cpp, line %d, MatGetOwnershipRange() returned error code %d",__LINE__,ierr);
    throw Comm_->ReportError(errMsg,-1);
  }

  PetscRowStart_ = rowStart;
  PetscRowEnd_   = rowEnd;

  int* MyGlobalElements = new int[rowEnd-rowStart];
  for (int i=0; i<rowEnd-rowStart; i++)
    MyGlobalElements[i] = rowStart+i;

  ierr = MatGetInfo(Amat,MAT_GLOBAL_SUM,&info);
  if (ierr) {
    sprintf(errMsg,"EpetraExt_PETScAIJMatrix.cpp, line %d, MatGetInfo() returned error code %d",__LINE__,ierr);
    throw Comm_->ReportError(errMsg,-1);
  }
  int tmp;
  ierr = MatGetSize(Amat,&NumGlobalRows_,&tmp);

  DomainMap_ = new Epetra_Map(NumGlobalRows_, NumMyRows_, MyGlobalElements, 0, *Comm_);

  // get the GIDs of the non-local columns
  //FIXME what if the matrix is sequential?

  int * ColGIDs = new int[NumMyCols_];
  for (int i=0; i<numLocalCols; i++) ColGIDs[i] = MyGlobalElements[i];
  for (int i=numLocalCols; i<NumMyCols_; i++) ColGIDs[i] = aij->garray[i-numLocalCols];

  ColMap_ = new Epetra_Map(-1, NumMyCols_, ColGIDs, 0, *Comm_);

  Importer_ = new Epetra_Import(*ColMap_, *DomainMap_);

  delete [] MyGlobalElements;
  delete [] ColGIDs;
} //Epetra_PETScAIJMatrix(Mat Amat)

//==============================================================================

Epetra_PETScAIJMatrix::~Epetra_PETScAIJMatrix(){
  if (ImportVector_!=0) delete ImportVector_;
  delete Importer_;
  delete ColMap_;
  delete DomainMap_;
  delete Comm_;

  if (Values_!=0) {delete [] Values_; Values_=0;}
  if (Indices_!=0) {delete [] Indices_; Indices_=0;}
} //Epetra_PETScAIJMatrix dtor

//==========================================================================

extern "C" {
  PetscErrorCode CreateColmap_MPIAIJ_Private(Mat);
}

int Epetra_PETScAIJMatrix::ExtractMyRowCopy(int Row, int Length, int & NumEntries, double * Values,
                     int * Indices) const 
{
  int nz;
  PetscInt *gcols, *lcols, ierr;
  PetscScalar *vals;
  bool alloc=false;

  // PETSc assumes the row number is global, whereas Trilinos assumes it's local.
  int globalRow = PetscRowStart_ + Row;
  assert(globalRow < PetscRowEnd_);
  ierr=MatGetRow(Amat_, globalRow, &nz, (const PetscInt**) &gcols, (const PetscScalar **) &vals);CHKERRQ(ierr);

  // I ripped this bit of code from PETSc's MatSetValues_MPIAIJ() in mpiaij.c.  The PETSc getrow returns everything in
  // global numbering, so we must convert to local numbering.
  if (strcmp(MatType_,MATMPIAIJ) == 0) {
    Mat_MPIAIJ  *aij = (Mat_MPIAIJ*)Amat_->data;
    lcols = (PetscInt *) malloc(nz * sizeof(int));
    alloc=true;
    if (!aij->colmap) {
      ierr = CreateColmap_MPIAIJ_Private(Amat_);CHKERRQ(ierr);
    }
    /*
      A PETSc parallel aij matrix uses two matrices to represent the local rows.
      The first matrix, A, is square and contains all local columns.
      The second matrix, B, is rectangular and contains all non-local columns.

      Matrix A:
      Local column ID's are mapped to global column id's by adding cmap.rstart.
      Given the global ID of a local column, the local ID is found by
      subtracting cmap.rstart.

      Matrix B:
      Non-local column ID's are mapped to global column id's by the local-to-
      global map garray.  Given the global ID of a local column, the local ID is
      found by the global-to-local map colmap.  colmap is either an array or
      hash table, the latter being the case when PETSC_USE_CTABLE is defined.
    */
    int offset = Amat_->cmap->n-1; //offset for non-local column indices

    for (int i=0; i<nz; i++) {
      if (gcols[i] >= Amat_->cmap->rstart && gcols[i] < Amat_->cmap->rend) {
        lcols[i] = gcols[i] - Amat_->cmap->rstart;
      } else {
#       ifdef PETSC_USE_CTABLE
        ierr = PetscTableFind(aij->colmap,gcols[i]+1,lcols+i);CHKERRQ(ierr);
        lcols[i] = lcols[i] + offset;
#       else
        lcols[i] = aij->colmap[gcols[i]] + offset;
#       endif
      }

    } //for i=0; i<nz; i++)
  }
  else lcols = gcols;

  NumEntries = nz;
  if (NumEntries > Length) return(-1);
  for (int i=0; i<NumEntries; i++) {
    Indices[i] = lcols[i];
    Values[i] = vals[i];
  }
  if (alloc) free(lcols);
  MatRestoreRow(Amat_,globalRow,&nz,(const PetscInt**) &gcols, (const PetscScalar **) &vals);
  return(0);
} //ExtractMyRowCopy()

//==========================================================================

int Epetra_PETScAIJMatrix::NumMyRowEntries(int Row, int & NumEntries) const 
{
  int globalRow = PetscRowStart_ + Row;
  MatGetRow(Amat_, globalRow, &NumEntries, PETSC_NULL, PETSC_NULL);
  MatRestoreRow(Amat_,globalRow,&NumEntries, PETSC_NULL, PETSC_NULL);
  return(0);
}

//==============================================================================

int Epetra_PETScAIJMatrix::ExtractDiagonalCopy(Epetra_Vector & Diagonal) const
{

  //TODO optimization: only get this diagonal once
  Vec petscDiag;
  double *vals=0;
  int length;

  int ierr=VecCreate(Comm_->Comm(),&petscDiag);CHKERRQ(ierr);
  VecSetSizes(petscDiag,NumMyRows_,NumGlobalRows_);
# ifdef HAVE_MPI
  ierr = VecSetType(petscDiag,VECMPI);CHKERRQ(ierr);
# else //TODO untested!!
  VecSetType(petscDiag,VECSEQ);
# endif

  MatGetDiagonal(Amat_, petscDiag);
  VecGetArray(petscDiag,&vals);
  VecGetLocalSize(petscDiag,&length);
  for (int i=0; i<length; i++) Diagonal[i] = vals[i];
  VecRestoreArray(petscDiag,&vals);
  VecDestroy(petscDiag);
  return(0);
}

//=============================================================================

int Epetra_PETScAIJMatrix::Multiply(bool TransA,
                               const Epetra_MultiVector& X,
                               Epetra_MultiVector& Y) const
{
  (void)TransA;
  int NumVectors = X.NumVectors();
  if (NumVectors!=Y.NumVectors()) EPETRA_CHK_ERR(-1);  // X and Y must have same number of vectors

  double ** xptrs;
  double ** yptrs;
  X.ExtractView(&xptrs);
  Y.ExtractView(&yptrs);
  if (RowMatrixImporter()!=0) {
    if (ImportVector_!=0) {
      if (ImportVector_->NumVectors()!=NumVectors) { delete ImportVector_; ImportVector_= 0;}
    }
    if (ImportVector_==0) ImportVector_ = new Epetra_MultiVector(RowMatrixColMap(),NumVectors);
    ImportVector_->Import(X, *RowMatrixImporter(), Insert);
    ImportVector_->ExtractView(&xptrs);
  }

  double *vals=0;
  int length;
  Vec petscX, petscY;
  int ierr;
  for (int i=0; i<NumVectors; i++) {
#   ifdef HAVE_MPI
    ierr=VecCreateMPIWithArray(Comm_->Comm(),X.MyLength(),X.GlobalLength(),xptrs[i],&petscX); CHKERRQ(ierr);
    ierr=VecCreateMPIWithArray(Comm_->Comm(),Y.MyLength(),Y.GlobalLength(),yptrs[i],&petscY); CHKERRQ(ierr);
#   else //FIXME  untested
    ierr=VecCreateSeqWithArray(Comm_->Comm(),X.MyLength(),X.GlobalLength(),xptrs[i],&petscX); CHKERRQ(ierr);
    ierr=VecCreateSeqWithArray(Comm_->Comm(),Y.MyLength(),Y.GlobalLength(),yptrs[i],&petscY); CHKERRQ(ierr);
#   endif

    ierr = MatMult(Amat_,petscX,petscY);CHKERRQ(ierr);

    ierr = VecGetArray(petscY,&vals);CHKERRQ(ierr);
    ierr = VecGetLocalSize(petscY,&length);CHKERRQ(ierr);
    for (int j=0; j<length; j++) yptrs[i][j] = vals[j];
    ierr = VecRestoreArray(petscY,&vals);CHKERRQ(ierr);
  }

  VecDestroy(petscX); VecDestroy(petscY);
  
  double flops = NumGlobalNonzeros();
  flops *= 2.0;
  flops *= (double) NumVectors;
  UpdateFlops(flops);
  return(0);
} //Multiply()

//=============================================================================

int Epetra_PETScAIJMatrix::Solve(bool Upper, bool Trans, bool UnitDiagonal, 
                            const Epetra_MultiVector& X,
                            Epetra_MultiVector& Y) const
{
  (void)Upper;
  (void)Trans;
  (void)UnitDiagonal;
  (void)X;
  (void)Y;
  return(-1); // Not implemented
}

//=============================================================================
// Utility routine to get the specified row and put it into Values_ and Indices_ arrays
int Epetra_PETScAIJMatrix::MaxNumEntries() const {
  int NumEntries;

  if (MaxNumEntries_==-1) {
    for (int i=0; i < NumMyRows_; i++) {
      NumMyRowEntries(i, NumEntries);
      if (NumEntries>MaxNumEntries_) MaxNumEntries_ = NumEntries;
    }
  }
  return(MaxNumEntries_);
}

//=============================================================================
// Utility routine to get the specified row and put it into Values_ and Indices_ arrays
int Epetra_PETScAIJMatrix::GetRow(int Row) const {

  int NumEntries;
  int MaxNumEntries = Epetra_PETScAIJMatrix::MaxNumEntries();

  if (MaxNumEntries>0) {
    if (Values_==0) Values_ = new double[MaxNumEntries];
    if (Indices_==0) Indices_ = new int[MaxNumEntries];
  }
  Epetra_PETScAIJMatrix::ExtractMyRowCopy(Row, MaxNumEntries, NumEntries, Values_, Indices_);
  
  return(NumEntries);
}

//=============================================================================
int Epetra_PETScAIJMatrix::InvRowSums(Epetra_Vector& x) const {
//
// Put inverse of the sum of absolute values of the ith row of A in x[i].
//

  if (!OperatorRangeMap().SameAs(x.Map())) EPETRA_CHK_ERR(-2); // x must have the same distribution as the range of A

  int ierr = 0;
  int i, j;
  for (i=0; i < NumMyRows_; i++) {
    int NumEntries = GetRow(i); // Copies ith row of matrix into Values_ and Indices_
    double scale = 0.0;
    for (j=0; j < NumEntries; j++) scale += fabs(Values_[j]);
    if (scale<Epetra_MinDouble) {
      if (scale==0.0) ierr = 1; // Set error to 1 to signal that zero rowsum found (supercedes ierr = 2)
      else if (ierr!=1) ierr = 2;
      x[i] = Epetra_MaxDouble;
    }
    else
      x[i] = 1.0/scale;
  }
  UpdateFlops(NumGlobalNonzeros());
  EPETRA_CHK_ERR(ierr);
  return(0);
}

//=============================================================================
int Epetra_PETScAIJMatrix::InvColSums(Epetra_Vector& x) const {
//
// Put inverse of the sum of absolute values of the jth column of A in x[j].
//

  if (!Filled()) EPETRA_CHK_ERR(-1); // Matrix must be filled.
  if (!OperatorDomainMap().SameAs(x.Map())) EPETRA_CHK_ERR(-2); // x must have the same distribution as the domain of A
  

  Epetra_Vector * xp = 0;
  Epetra_Vector * x_tmp = 0;
  

  // If we have a non-trivial importer, we must export elements that are permuted or belong to other processors
  if (RowMatrixImporter()!=0) {
    x_tmp = new Epetra_Vector(RowMatrixColMap()); // Create import vector if needed
    xp = x_tmp;
  }
  int ierr = 0;
  int i, j;

  for (i=0; i < NumMyCols_; i++) (*xp)[i] = 0.0;

  for (i=0; i < NumMyRows_; i++) {
    int NumEntries = GetRow(i);// Copies ith row of matrix into Values_ and Indices_
    for (j=0; j < NumEntries; j++) (*xp)[Indices_[j]] += fabs(Values_[j]);
  }

  if (RowMatrixImporter()!=0){
    x.Export(*x_tmp, *RowMatrixImporter(), Add); // Fill x with Values from import vector
    delete x_tmp;
    xp = &x;
  }
  // Invert values, don't allow them to get too large
  for (i=0; i < NumMyRows_; i++) {
    double scale = (*xp)[i];
    if (scale<Epetra_MinDouble) {
      if (scale==0.0) ierr = 1; // Set error to 1 to signal that zero rowsum found (supercedes ierr = 2)
      else if (ierr!=1) ierr = 2;
      (*xp)[i] = Epetra_MaxDouble;
    }
    else
      (*xp)[i] = 1.0/scale;
  }
  UpdateFlops(NumGlobalNonzeros());
  EPETRA_CHK_ERR(ierr);
  return(0);
}

//=============================================================================
int Epetra_PETScAIJMatrix::LeftScale(const Epetra_Vector& X) {
//
// This function scales the ith row of A by x[i].
//
  double *xptr;
  X.ExtractView(&xptr);
  Vec petscX;
# ifdef HAVE_MPI
  int ierr=VecCreateMPIWithArray(Comm_->Comm(),X.MyLength(),X.GlobalLength(),xptr,&petscX); CHKERRQ(ierr);
# else //FIXME  untested
  int ierr=VecCreateSeqWithArray(Comm_->Comm(),X.MyLength(),X.GlobalLength(),xptr,&petscX); CHKERRQ(ierr);
# endif

  MatDiagonalScale(Amat_, petscX, PETSC_NULL);

  ierr=VecDestroy(petscX); CHKERRQ(ierr);
  return(0);
} //LeftScale()

//=============================================================================
int Epetra_PETScAIJMatrix::RightScale(const Epetra_Vector& X) {
//
// This function scales the jth row of A by x[j].
//
  double *xptr;
  X.ExtractView(&xptr);
  Vec petscX;
# ifdef HAVE_MPI
  int ierr=VecCreateMPIWithArray(Comm_->Comm(),X.MyLength(),X.GlobalLength(),xptr,&petscX); CHKERRQ(ierr);
# else //FIXME  untested
  int ierr=VecCreateSeqWithArray(Comm_->Comm(),X.MyLength(),X.GlobalLength(),xptr,&petscX); CHKERRQ(ierr);
# endif

  MatDiagonalScale(Amat_, PETSC_NULL, petscX);

  ierr=VecDestroy(petscX); CHKERRQ(ierr);
  return(0);
} //RightScale()

//=============================================================================

double Epetra_PETScAIJMatrix::NormInf() const {

  if (NormInf_>-1.0) return(NormInf_);

  MatNorm(Amat_, NORM_INFINITY,&NormInf_);
  UpdateFlops(NumGlobalNonzeros());

  return(NormInf_);
}

//=============================================================================

double Epetra_PETScAIJMatrix::NormOne() const {

  if (NormOne_>-1.0) return(NormOne_);
  if (!Filled()) EPETRA_CHK_ERR(-1); // Matrix must be filled.

  MatNorm(Amat_, NORM_1,&NormOne_);
  UpdateFlops(NumGlobalNonzeros());

  return(NormOne_);
}

//=============================================================================

void Epetra_PETScAIJMatrix::Print(ostream& os) const {
  return;
}

#endif /*HAVE_PETSC*/
