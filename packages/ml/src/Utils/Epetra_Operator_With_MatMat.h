//@HEADER
/*!
 * \file Epetra_Operator_With_MatMat.h
 *
 * \class Epetra_Operator_With_MatMat
 *
 * \brief  A pure virtual class, derived from Epetra_Operator that adds a
 * matmat routine (using Epetra_CrsMatrix).
 *
 * \date Last update to Doxygen: 25-Jan-07
 *
 */
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
#ifndef EPETRA_OPERATOR_WITH_MATMAT_H
#define EPETRA_OPERATOR_WITH_MATMAT_H

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif
#if defined(HAVE_ML_EPETRA)
#include "Epetra_Operator.h"
#include "Epetra_CrsMatrix.h"
#include "ml_include.h"

namespace ML_Epetra{

/*! The Epetra_Operator_With_MatMat class is a pure virtual class (specifies interface only) that
  enable the use of real-valued double-precision operators.  It's only addition
  to the Epetra_Operator is that it adds matmat operations with
  Epetra_CrsMatrices, returning either an EpetraCrsMatrix or an ML_Operator.
*/
class Epetra_Operator_With_MatMat: public Epetra_Operator {

public:

  //! @name Destructor
  //@{
  //! Destructor
  virtual ~Epetra_Operator_With_MatMat() {};
  //@}


  //@{ Mathematical Functions
  //! Computes C= <me> * A
  virtual int MatrixMatrix_Multiply(const Epetra_CrsMatrix & A, Epetra_CrsMatrix **C) const=0;

  //! Computes C= <me> * A
  virtual int MatrixMatrix_Multiply(const Epetra_CrsMatrix & A, ML_Comm *comm,  ML_Operator **C) const=0;

  //@}
};

}/*end namespace*/

#endif
#endif /* EPETRA_OPERATOR_WITH_MATMAT_H */
