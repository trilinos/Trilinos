#ifndef ML_EXPRESSIONS_H
#define ML_EXPRESSIONS_H

#include "MLAPI_Error.h"
#include "MLAPI_Workspace.h"
#include "MLAPI_Space.h"
#include "MLAPI_MultiVector.h"
#include "MLAPI_Operator.h"
#include "MLAPI_InverseOperator.h"

namespace MLAPI {

/*!
\file MLAPI_Expressions.h

\brief Overloaded operators for MultiVector's, Operator's, and InverseOpereator's.

\author Marzio Sala, SNL 9214.

\date Last updated on Feb-05.
*/

// ====================================================================== 
//! Creates a new MultiVector, defined as x + y
// ====================================================================== 

MultiVector operator+(const MultiVector& x, const MultiVector& y);

// ====================================================================== 
//! Creates a new MultiVector, defined as x - y
// ====================================================================== 

MultiVector operator-(const MultiVector& x, const MultiVector& y);

// ====================================================================== 
//! Creates a new MultiVector, defined as x + y
// ====================================================================== 

MultiVector operator+(const MultiVector& x, const double alpha);

// ====================================================================== 
//! Creates a new MultiVector, defined as x - y
// ====================================================================== 

MultiVector operator-(const MultiVector& x, const double alpha);

// ====================================================================== 
//! Creates a new MultiVector, defined as x + y
// ====================================================================== 

MultiVector operator+(const double alpha, const MultiVector& x);

// ====================================================================== 
//! Creates a new MultiVector, defined as x - y
// ====================================================================== 

MultiVector operator-(const double alpha, const MultiVector& x);

#if 0
// ====================================================================== 
//! Adds a constant to a vector
// ====================================================================== 

MultiVector operator+= (const double alpha);

// ====================================================================== 
//! Subtracts a constant to a vector
// ====================================================================== 

MultiVector operator-= (const double alpha);
#endif

// ====================================================================== 
//! Creates a new Operator, defined as A + B
// ====================================================================== 

Operator operator+(const Operator& A, const Operator& B);

// ====================================================================== 
//! Creates a new Operator, defined as A - B
// ====================================================================== 

Operator operator-(const Operator& A, const Operator& B);

// ====================================================================== 
//! Creates a new Operator, defined as A * B
// ====================================================================== 

Operator operator*(const Operator& A, const Operator& B);

// ====================================================================== 
//! Creates a new Operator, defined as A * alpha
// ====================================================================== 

Operator operator*(const Operator& A, const double alpha);

// ====================================================================== 
//! Creates a new Operator, defined as alpha * A
// ====================================================================== 

Operator operator*(const double alpha, const Operator& A);

// ====================================================================== 
//! Creates a new MultiVector, defined as x * alpha
// ====================================================================== 

MultiVector operator*(const MultiVector& x, const double alpha);

// ====================================================================== 
//! Creates a new MultiVector y, such that y = x / alpha
// ====================================================================== 

MultiVector operator/(const MultiVector& x, const double alpha);

// ====================================================================== 
//! Creates a new MultiVector y, such that y = A * x.
// ====================================================================== 

MultiVector operator*(const BaseOperator& A, const MultiVector& x);

// ====================================================================== 
//! Computes the dot product between the first vector in x and y
// ====================================================================== 

double operator* (const MultiVector& x, const MultiVector& y);


} // namespace MLAPI

#endif // if ML_EXPRESSIONS_H
