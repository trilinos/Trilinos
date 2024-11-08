#ifndef ML_OPERATOR_UTILS_H
#define ML_OPERATOR_UTILS_H

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif
/*!
\file MLAPI_Operator_Utils

\brief Functions to create and compose MLAPI::Operator's.

\author Marzio Sala, SNL 9214.

\date Last updated on Feb-05.
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
#include "ml_common.h"

#include "ml_include.h"
#include <iostream>
#include "ml_operator.h"
#include "ml_epetra.h"
#include "ml_amesos.h"
#include "ml_epetra_utils.h"
#include "ml_amesos_wrap.h"
#ifdef HAVE_ML_ANASAZI
#include "ml_anasazi.h"
#endif
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "MLAPI_Error.h"
#include "MLAPI_Space.h"
#include "MLAPI_MultiVector.h"
#include "MLAPI_Operator.h"

namespace Teuchos {
  class ParameterList;
}

namespace MLAPI {

//! Performs a triple matrix-matrix product, res = R * A *P.
Operator GetRAP(const Operator& R, const Operator& A,
                const Operator& P);

//! Returns a newly created transpose of \c A.
Operator GetTranspose(const Operator& A, const bool byrow = true);

//! Returns the identity matrix.
Operator GetIdentity(const Space& DomainSpace, const Space& RangeSpace);

//! Returns a vector containing the diagonal elements of \c A.
MultiVector GetDiagonal(const Operator& A);

//! Returns a vector containing the diagonal elements of \c A.
MultiVector GetDiagonal(const Operator& A, const int offset);

//! Returns a newly created operator, containing D on the diagonal.
Operator GetDiagonal(const MultiVector& D);

//! Returns an operator defined as (I - Damping A).
Operator GetJacobiIterationOperator(const Operator& Amat, double Damping);

//! Returns a newly created operator, containing D on the diagonal.
Operator GetPtent1D(const MultiVector& D, const int offset = 0);

// Creates C = scalarA * A + scalarB * B.
int ML_Operator_Add2(ML_Operator *A, ML_Operator *B, ML_Operator *C,
		    int matrix_type, double scalarA, double scalarB);

//! Performs a cheap analysis of the properties of the input operator.
void AnalyzeCheap(const Operator& A);

//! Prints on file the sparsity structure of input operator.
void PrintSparsity(const Operator& A, int NumPDEEquations = 1);

//! Multiply A by a double value, \c alpha.
Operator GetScaledOperator(const Operator& A, const double alpha);

//! Duplicates a given operator.
Operator Duplicate(const Operator& A);

} // namespace MLAPI

#endif // ML_OPERATOR_UTILS_H
