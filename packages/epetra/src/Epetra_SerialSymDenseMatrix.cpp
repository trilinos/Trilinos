
//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright 2001 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// Export of this program may require a license from the United States
// Government.
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

