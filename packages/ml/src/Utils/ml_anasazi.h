#ifndef _ML_ANASAZI_H_
#define _ML_ANASAZI_H_

#include "ml_include.h"

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_ANASAZI)

// anasazi interface

#include "Teuchos_ParameterList.hpp"

namespace ML_Anasazi {

extern int Interface(const Epetra_RowMatrix * RowMatrix, Epetra_MultiVector & EigenVectors,
				double RealEigenvalues[], double ImagEigenvalues[],
				Teuchos::ParameterList & List);


extern int GetFieldOfValuesBox(const Epetra_RowMatrix * RowMatrix, 
					   double & MaxReal, double & MaxImag,
					   Teuchos::ParameterList & AnasaziList);

}

#endif

extern "C" {
  
  extern int ML_Anasazi_Get_FiledOfValuesBox_Interface(ML_Operator * Amat,
						       struct ML_Field_Of_Values * fov );
  extern int ML_Anasazi_Get_SpectralNorm_Anasazi(ML_Operator * Amat,
						 int MaxIters, double Tolerance,
						 int IsProblemSymmetric,
						 int UseDiagonalScaling,
						 double * LambdaMax );
}

#endif /* #ifndef _ML_ANASAZI_H_ */
