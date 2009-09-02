
//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
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
// ************************************************************************
//@HEADER

#include "Epetra_SerialSymDenseMatrix.h"
//=============================================================================
Epetra_SerialSymDenseMatrix::Epetra_SerialSymDenseMatrix(void)
  : Epetra_SerialDenseMatrix(),
    Upper_(false),
    UPLO_('L')

{
}
//=============================================================================
Epetra_SerialSymDenseMatrix::Epetra_SerialSymDenseMatrix(Epetra_DataAccess CV, double *A, int LDA, int NumRowsCols)
  : Epetra_SerialDenseMatrix(CV, A, LDA, NumRowsCols, NumRowsCols),
    Upper_(false),
    UPLO_('L')

{
}
//=============================================================================
Epetra_SerialSymDenseMatrix::Epetra_SerialSymDenseMatrix(const Epetra_SerialSymDenseMatrix& Source)
  : Epetra_SerialDenseMatrix(Source),  
    Upper_(Source.Upper_),
    UPLO_(Source.UPLO_)
{
}
//=============================================================================
Epetra_SerialSymDenseMatrix::~Epetra_SerialSymDenseMatrix()
{
}
//=============================================================================
void Epetra_SerialSymDenseMatrix::CopyUPLOMat(bool Upper, double * A, int LDA, int NumRows) {

  int i, j;
  double * ptr1;
  double * ptr2;

  if (Upper) {
    for (j=1; j<NumRows; j++) {
      ptr1 = A + j;
      ptr2 = A + j*LDA;
      for (i=0; i<j; i++) {
	*ptr1 = *ptr2++;
	ptr1+=LDA;
      }
    }
  }
  else {
    for (i=1; i<NumRows; i++) {
      ptr1 = A + i;
      ptr2 = A + i*LDA;
      for (j=0; j<i; j++) {
	*ptr2++ = *ptr1;
	ptr1+=LDA;
      }
    }
  }
}
//=============================================================================
double Epetra_SerialSymDenseMatrix::NormOne(void) const{

    return(Epetra_SerialSymDenseMatrix::NormInf());
}
//=============================================================================
double Epetra_SerialSymDenseMatrix::NormInf(void) const {

  int i, j;

  double anorm = 0.0;
  double * ptr;

  if (!Upper()) {
    for (j=0; j<N_; j++) {
      double sum = 0.0;
      ptr = A_ + j + j*LDA_;
      for (i=j; i<N_; i++) sum += std::abs(*ptr++);
      ptr = A_ + j;
      for (i=0; i<j; i++) {
	sum += std::abs(*ptr);
	ptr += LDA_;
      }
      anorm = EPETRA_MAX(anorm, sum);
    }
  }
  else {
    for (j=0; j<N_; j++) {
      double sum = 0.0;
      ptr = A_ + j*LDA_;
      for (i=0; i<j; i++) sum += std::abs(*ptr++);
      ptr = A_ + j + j*LDA_;
      for (i=j; i<N_; i++) {
	sum += std::abs(*ptr);
	ptr += LDA_;
      }
      anorm = EPETRA_MAX(anorm, sum);
    }
  }
  UpdateFlops(N_*N_);
  return(anorm);
}

//=============================================================================
int Epetra_SerialSymDenseMatrix::Scale(double ScalarA) {

  int i, j;

  double * ptr;

  if (!Upper()) {
    for (j=0; j<N_; j++) {
      ptr = A_ + j + j*LDA_;
      for (i=j; i<N_; i++) {*ptr = *ptr * ScalarA; ptr++;}
    }
  }
  else {
    for (j=0; j<N_; j++) {
      ptr = A_ + j*LDA_;
      for (i=0; i<j; i++) {*ptr = *ptr * ScalarA; ptr++;}
    }
  }
  UpdateFlops(N_*(N_+1)/2);
  return(0);
}

