// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

//-----------------------------------------------------
// Ifpack2::BorderedOperator is a translation of the
// LOCA_BorderedSolver_EpetraHouseholder
// implementation written by Eric Phipps.
// DMD.
//------------------------------------------------------


#ifndef IFPACK2_BORDEREDOPERATOR_DECL_HPP
#define IFPACK2_BORDEREDOPERATOR_DECL_HPP

#include "Ifpack2_ConfigDefs.hpp"
#include "Tpetra_Operator.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include <iostream>

namespace Ifpack2 {
//! Ifpack2 bordered operator

/*!
  Ifpack2::BorderedOperator is a concrete class
extending Tpetra::Operator and defining an interface.

  This class extends Tpetra::Operator, providing the additional methods:
  - compute() computes everything required to apply the
    bordered operator, using matrix values (and assuming that the
    sparsity of the matrix has not been changed);
  - isComputed() should return true if the bordered operator
    has been successfully computed, false otherwise.
  - getRHS() returns a reference to the matrix to be preconditioned.
  - getLHS() returns a reference to the matrix to be preconditioned.

The bordered operator is applied by apply()
(which returns if isComputed() is false).
Each time compute() is called, the object re-computes the actual values of
the bordered operator.

<b>Title of Method Description</b>
*/

template<class Scalar,
         class LocalOrdinal = typename Tpetra::Operator<Scalar>::local_ordinal_type,
         class GlobalOrdinal = typename Tpetra::Operator<Scalar, LocalOrdinal>::global_ordinal_type,
         class Node = typename Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal>::node_type>
class BorderedOperator :
    virtual public Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node > {
public:
  //! Constructor with Tpetra::Operator input.
  BorderedOperator (const Teuchos::RCP<const Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& A);

  //! Destructor.
  virtual ~BorderedOperator() {}

  //! \name Implementation of Tpetra::Operator.
  //@{

  //! The domain Map of this operator.  It must be compatible with X.getMap().
  virtual Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > getDomainMap() const;

  //! The range Map of this operator.  It must be compatible with Y.getMap().
  virtual Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > getRangeMap() const;

  //! Whether this operator can apply the transpose or conjugate transpose.
  bool hasTransposeApply() const;

  //! Apply the bordered operator.
  void
  apply (const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X,
         Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y,
         Teuchos::ETransp mode = Teuchos::NO_TRANS,
         Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
         Scalar beta = Teuchos::ScalarTraits<Scalar>::zero()) const;
  //@}

private:
  //! The Operator with which this instance was constructed.
  Teuchos::RCP<const Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> > A_;
};

}//namespace Ifpack2

#endif // IFPACK2_BORDEREDOPERATOR_DECL_HPP
