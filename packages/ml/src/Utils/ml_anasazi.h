#ifndef _ML_ANASAZI_H_
#define _ML_ANASAZI_H_

#include "ml_include.h"

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS)

// anasazi interface

extern int ML_Anasazi_Interface(const Epetra_RowMatrix * RowMatrix, int & NullSpaceDim,
				 double * & NullSpacePtr, ParameterList & List);

#endif

#endif /* #ifndef _ML_ANASAZI_H_ */
