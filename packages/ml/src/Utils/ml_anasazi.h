#ifndef _ML_ANASAZI_H_
#define _ML_ANASAZI_H_

#include "ml_include.h"
#include "ml_struct.h"

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_ANASAZI)

// anasazi interface

#include "Teuchos_ParameterList.hpp"

//! ML_Anasazi: default namespace for all Anasazi interfaces.

namespace ML_Anasazi {

  //! \brief Interface to Anasazi to eigen-computations.
  /*!
    This function creates an interface from an Epetra_RowMatrix and Anasazi.
    Parameters are specified using a Teuchos parameters list.

    \param \in RowMatrix: Epetra_RowMatrix;

    \param \out EigenVectors: Epetra_MultiVectors, that will contain the
    requested eigenvectors;

    \param \out RealEigenvalues: double array that will contain the real
    part of the computed eigenvalues;
    
    \param \out ImagEigenvalues: double array that will contain the
    imaginary part of the computed eigenvalues;

    \param \in List: Teuchos parameters' list containing the required options.
    
  */
extern int Interface(const Epetra_RowMatrix * RowMatrix, Epetra_MultiVector & EigenVectors,
		     double RealEigenvalues[], double ImagEigenvalues[],
		     Teuchos::ParameterList & List,
		     double RealEigenvectors[] = 0, double ImagEigenvectors[] = 0,
		     int * NumRealEigenvectors = 0, int * NumImagEigenvectors = 0,
		     ML * ml = 0);

  //! Computes the size of a box containing the field of values.
extern int GetFieldOfValuesBox(const Epetra_RowMatrix * RowMatrix, 
			       double & MaxReal, double & MaxImag,
			       Teuchos::ParameterList & AnasaziList);

}

#endif

extern "C" {

  //! Interface from C code to Anasazi to compute the field of values. 
  extern int ML_Anasazi_Get_FieldOfValuesBox_Interface(ML_Operator * Amat,
						       struct ML_Field_Of_Values * fov );

  //! Interface from C code to Anasazi to compute the maximum eigenvalue.
  extern int ML_Anasazi_Get_SpectralNorm_Anasazi(ML_Operator * Amat,
						 int MaxIters, double Tolerance,
						 int IsProblemSymmetric,
						 int UseDiagonalScaling,
						 double * LambdaMax );
}

#endif /* #ifndef _ML_ANASAZI_H_ */
