#ifndef _ML_ANASAZI_H_
#define _ML_ANASAZI_H_

#include "ml_include.h"

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_ANASAZI)

// anasazi interface

class Teuchos::ParameterList;

namespace ML_Anasazi {

extern int Interface(const Epetra_RowMatrix * RowMatrix, Epetra_MultiVector & EigenVectors,
				double RealEigenvalues[], double ImagEigenvalues[],
				Teuchos::ParameterList & List);


extern int GetFieldOfValuesBox(const Epetra_RowMatrix * RowMatrix, 
					   double & MaxReal, double & MaxImag,
					   Teuchos::ParameterList & AnasaziList);

}

extern "C" {
  
extern int ML_Anasazi_Get_FiledOfValuesBox_Interface(ML_Operator * Amat,
						     struct ML_Field_Of_Values * fov );
  
}

#endif

#endif /* #ifndef _ML_ANASAZI_H_ */
