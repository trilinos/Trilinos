
/* Copyright (2001) Sandia Corportation. Under the terms of Contract 
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2001, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 * 
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */


#include "Epetra_IntSerialDenseMatrix.h"
//=============================================================================
Epetra_IntSerialDenseMatrix::Epetra_IntSerialDenseMatrix(void)
  : Epetra_Object(),
    M_(0),
    N_(0),
    LDA_(0),
    A_Copied_(false),
    A_(0)
{
}

//=============================================================================
Epetra_IntSerialDenseMatrix::Epetra_IntSerialDenseMatrix(Epetra_DataAccess CV, int *A, int LDA, 
					     int NumRows, int NumCols)
  : Epetra_Object(),
    M_(NumRows),
    N_(NumCols),
    LDA_(LDA),
    A_Copied_(false),    
    A_(A)

{
  if (CV==Copy) {
    LDA_ = M_;
    A_ = new int[LDA_*N_];
    CopyMat(A, LDA, M_, N_, A_, LDA_);
    A_Copied_ = true;
  }

}
//=============================================================================
Epetra_IntSerialDenseMatrix::Epetra_IntSerialDenseMatrix(const Epetra_IntSerialDenseMatrix& Source)
  : Epetra_Object(Source),  
    M_(Source.M_),
    N_(Source.N_),
    LDA_(Source.LDA_),
    A_Copied_(true),
    A_(Source.A_)

{

  LDA_ = M_;
  A_ = new int[LDA_*N_];
  CopyMat(Source.A_, Source.LDA_, M_, N_, A_, LDA_);
}

//=============================================================================
int Epetra_IntSerialDenseMatrix::Reshape(int NumRows, int NumCols) {

  // Allocate space for new matrix
  int * A_tmp = new int[NumRows*NumCols];
  for (int k = 0; k < NumRows*NumCols; k++) A_tmp[k] = 0; // Zero out values
  int M_tmp = EPETRA_MIN(M_, NumRows);
  int N_tmp = EPETRA_MIN(N_, NumCols);
  if (A_ != 0) CopyMat(A_, LDA_, M_tmp, N_tmp, A_tmp, NumRows); // Copy principal submatrix of A to new A
  
  DeleteArrays(); // Get rid of anything that might be already allocated  
  M_ = NumRows;
  N_ = NumCols;
  LDA_ = M_;
  A_ = A_tmp; // Set pointer to new A
  A_Copied_ = true;

  return(0);
}
int Epetra_IntSerialDenseMatrix::Shape(int NumRows, int NumCols) {
  DeleteArrays(); // Get rid of anything that might be already allocated
  M_ = NumRows;
  N_ = NumCols;
  LDA_ = M_;
  A_ = new int[LDA_*N_];
  for (int k = 0; k < LDA_*N_; k++) A_[k] = 0; // Zero out values
  A_Copied_ = true;

  return(0);
}
//=============================================================================
Epetra_IntSerialDenseMatrix::~Epetra_IntSerialDenseMatrix()
{
  DeleteArrays();
}
//=============================================================================
void Epetra_IntSerialDenseMatrix::DeleteArrays(void)
{
  if (A_Copied_)   {delete [] A_; A_ = 0; A_Copied_ = false;}
}
//=============================================================================
Epetra_IntSerialDenseMatrix & Epetra_IntSerialDenseMatrix::operator = ( const Epetra_IntSerialDenseMatrix & Source) {
  if (this==&Source) return(*this); // Special case of source same as target

  if (M()!=Source.M()) throw ReportError("Row dimension of source = " + toString(Source.M()) +
					  " is different than  row dimension of target = " + toString(LDA()), -1);
  if (N()!=Source.N()) throw ReportError("Column dimension of source = " + toString(Source.N()) +
					  " is different than column dimension of target = " + toString(N()), -2);

  CopyMat(Source.A(), Source.LDA(), Source.M(), Source.N(), A(), LDA());
  return(*this);
}
//=============================================================================
void Epetra_IntSerialDenseMatrix::CopyMat(int * A, int LDA, int NumRows, int NumCols, 
					int * B, int LDB) {

  int i, j;
  int * ptr1 = B;
  int * ptr2;
  for (j=0; j<NumCols; j++) {
    ptr1 = B + j*LDB;
    ptr2 = A + j*LDA;
    for (i=0; i<NumRows; i++) *ptr1++ = *ptr2++;
  }
  return;
}
//=============================================================================
int Epetra_IntSerialDenseMatrix::OneNorm(void) {

  int i, j;

    int anorm = 0;
    int * ptr;
    for (j=0; j<N_; j++) {
      int sum=0;
      ptr = A_ + j*LDA_;
      for (i=0; i<M_; i++) sum += abs(*ptr++);
      anorm = EPETRA_MAX(anorm, sum);
    }
    return(anorm);
}
//=============================================================================
int Epetra_IntSerialDenseMatrix::InfNorm(void) {

  int i, j;

    int anorm = 0;
    int * ptr;

    // Loop across columns in inner loop.  Most expensive memory access, but 
    // requires no extra storage.
    for (i=0; i<M_; i++) {
      int sum=0;
      ptr = A_ + i;
      for (j=0; j<N_; j++) {
	sum += abs(*ptr);
	ptr += LDA_;
      }
      anorm = EPETRA_MAX(anorm, sum);
    }
    return(anorm);
}
//=========================================================================
int& Epetra_IntSerialDenseMatrix::operator () (int RowIndex, int ColIndex)  {

  if (RowIndex>=M_) throw ReportError("Row index = " +toString(RowIndex) + 
				      " Out of Range 0 - " + toString(M_-1),-1);
  if (ColIndex>=N_) throw ReportError("Column index = " +toString(ColIndex) + 
				      " Out of Range 0 - " + toString(N_-1),-2);

  return(A_[ColIndex*LDA_ + RowIndex]);
}

//=========================================================================
const int& Epetra_IntSerialDenseMatrix::operator () (int RowIndex, int ColIndex) const  {

  if (RowIndex>=M_) throw ReportError("Row index = " +toString(RowIndex) + 
				      " Out of Range 0 - " + toString(M_-1),-1);
  if (ColIndex>=N_) throw ReportError("Column index = " +toString(ColIndex) + 
				      " Out of Range 0 - " + toString(N_-1),-2);

   return(A_[ColIndex*LDA_ + RowIndex]);
}
//=========================================================================
const int* Epetra_IntSerialDenseMatrix::operator [] (int ColIndex) const  {

  if (ColIndex>=N_) throw ReportError("Column index = " +toString(ColIndex) + 
				      " Out of Range 0 - " + toString(N_-1),-2);
  return(A_ + ColIndex*LDA_);
}
//=========================================================================
int* Epetra_IntSerialDenseMatrix::operator [] (int ColIndex)  {

  if (ColIndex>=N_) throw ReportError("Column index = " +toString(ColIndex) + 
				      " Out of Range 0 - " + toString(N_-1),-2);
  return(A_+ ColIndex*LDA_);
}

void Epetra_IntSerialDenseMatrix::Print(ostream& os) const {

  os << endl;
  for (int i=0; i<M_; i++) {
    for (int j=0; j<N_; j++){
      os << (*this)(i,j) << " ";
    }
    os << endl;
  }
  return;
}
