#ifndef _ML_ANASAZI_H_
#define _ML_ANASAZI_H_

#include "ml_include.h"

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_ANASAZI)

// anasazi interface

extern int ML_Anasazi_Interface(const Epetra_RowMatrix * RowMatrix, Epetra_MultiVector & EigenVectors,
				double RealEigenvalues[], double ImagEigenvalues[],
				Teuchos::ParameterList & List);


extern int ML_Anasazi_Get_FiledOfValuesBox(const Epetra_RowMatrix * RowMatrix, 
					   double & MaxReal, double & MaxImag );

#endif

#endif /* #ifndef _ML_ANASAZI_H_ */
