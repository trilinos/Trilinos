
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
double Epetra_SerialSymDenseMatrix::OneNorm(void) {

    return(InfNorm());
}
//=============================================================================
double Epetra_SerialSymDenseMatrix::InfNorm(void) {

  int i, j;

  double anorm = 0.0;
  double * ptr;

  if (!Upper()) {
    for (j=0; j<N_; j++) {
      double sum = 0.0;
      ptr = A_ + j + j*LDA_;
      for (i=j; i<N_; i++) sum += fabs(*ptr++);
      ptr = A_ + j;
      for (i=0; i<j; i++) {
	sum += fabs(*ptr);
	ptr += LDA_;
      }
      anorm = EPETRA_MAX(anorm, sum);
    }
  }
  else {
    for (j=0; j<N_; j++) {
      double sum = 0.0;
      ptr = A_ + j*LDA_;
      for (i=0; i<j; i++) sum += fabs(*ptr++);
      ptr = A_ + j + j*LDA_;
      for (i=j; i<N_; i++) {
	sum += fabs(*ptr);
	ptr += LDA_;
      }
      anorm = EPETRA_MAX(anorm, sum);
    }
  }
  UpdateFlops(N_*N_);
  return(anorm);
}

