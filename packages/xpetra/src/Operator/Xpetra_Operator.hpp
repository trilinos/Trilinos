// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef XPETRA_OPERATOR_HPP
#define XPETRA_OPERATOR_HPP

#include "Xpetra_ConfigDefs.hpp"

#include <Teuchos_Describable.hpp>
#include <Teuchos_BLAS_types.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include "Xpetra_Map.hpp"
#include "Xpetra_MultiVector.hpp"

namespace Xpetra {

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
class Operator : virtual public Teuchos::Describable {
  typedef Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> map_type;
  typedef Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> mv_type;

 public:
  virtual ~Operator() {}
  /** \name Typedefs that give access to the template parameters. */
  //@{

  //! The type of the entries of the input and output multivectors.
  typedef Scalar scalar_type;

  //! The local index type.
  typedef LocalOrdinal local_ordinal_type;

  //! The global index type.
  typedef GlobalOrdinal global_ordinal_type;

  //! The Kokkos Node type.
  typedef Node node_type;

  //@}
  /** \name Pure virtual functions to be overridden by subclasses. */
  //@{

  //! The Map associated with the domain of this operator, which must be compatible with X.getMap().
  virtual const Teuchos::RCP<const map_type> getDomainMap() const = 0;

  //! The Map associated with the range of this operator, which must be compatible with Y.getMap().
  virtual const Teuchos::RCP<const map_type> getRangeMap() const = 0;

  //! \brief Computes the operator-multivector application.
  /*! Loosely, performs \f$Y = \alpha \cdot A^{\textrm{mode}} \cdot X + \beta \cdot Y\f$. However, the details of operation
      vary according to the values of \c alpha and \c beta. Specifically
      - if <tt>beta == 0</tt>, apply() <b>must</b> overwrite \c Y, so that any values in \c Y (including NaNs) are ignored.
      - if <tt>alpha == 0</tt>, apply() <b>may</b> short-circuit the operator, so that any values in \c X (including NaNs) are ignored.
   */
  virtual void
  apply(const mv_type& X, mv_type& Y,
        Teuchos::ETransp mode = Teuchos::NO_TRANS,
        Scalar alpha          = Teuchos::ScalarTraits<Scalar>::one(),
        Scalar beta           = Teuchos::ScalarTraits<Scalar>::zero()) const = 0;

  /// \brief Whether this operator supports applying the transpose or conjugate transpose.
  ///
  /// By default, this returns false.  Subclasses must override this
  /// method if they can support apply() with
  /// <tt>mode=Teuchos::TRANS</tt> or
  /// <tt>mode=Teuchos::CONJ_TRANS</tt>.
  virtual bool hasTransposeApply() const { return false; }

  //@}

  virtual void removeEmptyProcessesInPlace(const RCP<const map_type>& /* newMap */) {}

  //! Compute a residual R = B - (*this) * X
  virtual void residual(const mv_type& X,
                        const mv_type& B,
                        mv_type& R) const = 0;
};

}  // namespace Xpetra

#define XPETRA_OPERATOR_SHORT
#endif  // XPETRA_OPERATOR_HPP
