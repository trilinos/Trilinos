/* Copyright (2004) Sandia Corportation. Under the terms of Contract
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

/*!
 *  \file ml_anasazi.h
 *
 *  \brief Interface to the Trilinos package Anasazi.
 *
 *  Anasazi can be used by ML to eigenvalue computations. The user can:
 *  - use Anasazi to estimate the maximum eigenvalue of a given level matrix
 *  (for smoothed aggregation);
 *  - use Anasazi  to compute the low-convergence modes, and filter them;
 *  - use Anasazi to compute the field-of-values of a non-symmetric operator,
 *    and use this information to improve smoothed aggregation for highly
 *    non-symmetric systems.
 *    
 *  \author Marzio Sala, SNL (9214)
 *  
 *  \date Last update to Doxygen: 22-Jul-04
 *
 */

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

    \param RowMatrix (In) : Epetra_RowMatrix;

    \param EigenVectors (Out) : Epetra_MultiVectors, that will contain the
           requested eigenvectors;

    \param RealEigenvalues (Out) : double array that will contain the real
           part of the computed eigenvalues;
    
    \param ImagEigenvalues (Out) : double array that will contain the
           imaginary part of the computed eigenvalues;

    \param RealEigenvectors (Out) : pointer to allocated space to store the 
           real eigenvectors

    \param ImagEigenvectors (Out) : pointer to allocated space to store the 
           imaginary eigenvectors

    \param NumRealEigenvectors (Out) : number of computed real eigenvectors

    \param NumImagEigenvectors (Out) : number of computed imaginary eigenvectors

    \param List (InOut) : Teuchos parameters' list containing the required options.

    \param ml (In) : already filled ML hierarchy (only for the eigen-analysis 
           of "I-ML^{-1}A").
    
    The following parameters parsed from List:
    - "eigen-analysis: use diagonal scaling". Enables 
      scaling by the diagonal. Default: \c true.
    - "eigen-analysis: symmetric problem" . Instructs Anasazi to use an 
      algorithm for symmetric problems. Default: \c false.
    - "eigen-analysis: matrix operation". Defined the matrix-vector product. 
      Possible values are:
      - "A"
      - "I-A"
      - "A+A^T"
      - "A-A^T"
      - "I-ML^{-1}A" (this defines no scaling, and "eigen-analysis: symmetric problem" = \c false.
      Default: "A"
    - "eigen-analysis: length", 20)
    - "eigen-analysis: block-size", 1)
    - "eigen-analysis: tolerance", 1.0e-5)
    - "eigen-analysis: action", "LM"
    - "eigen-analysis: restart", 100)
    - "eigen-analysis: output". Defined the output level, from 0 to 10 (1- being verbose).
  */
 int Interface(const Epetra_RowMatrix * RowMatrix, Epetra_MultiVector & EigenVectors,
		     double RealEigenvalues[], double ImagEigenvalues[],
		     Teuchos::ParameterList & List,
		     double RealEigenvectors[] = 0, double ImagEigenvectors[] = 0,
		     int * NumRealEigenvectors = 0, int * NumImagEigenvectors = 0,
		     ML * ml = 0);

  //! Computes the size of a box containing the field of values.
 int GetFieldOfValuesBox(const Epetra_RowMatrix * RowMatrix, 
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
