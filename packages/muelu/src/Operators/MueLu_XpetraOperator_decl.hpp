// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_XPETRAOPERATOR_DECL_HPP
#define MUELU_XPETRAOPERATOR_DECL_HPP

#include "MueLu_ConfigDefs.hpp"

#include <Xpetra_Operator_fwd.hpp>
#include <Xpetra_MultiVector_fwd.hpp>
#include "MueLu_Level.hpp"
#include "MueLu_Hierarchy_decl.hpp"

namespace MueLu {

/*!  @brief Wraps an existing MueLu::Hierarchy as a Xpetra::Operator.
 */
template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class XpetraOperator : public Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
 protected:
  XpetraOperator() {}

 public:
  //! @name Constructor/Destructor
  //@{

  //! Constructor
  XpetraOperator(const RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& H);

  //! Destructor.
  virtual ~XpetraOperator();

  //@}

  //! Returns the Tpetra::Map object associated with the domain of this operator.
  const Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > getDomainMap() const;

  //! Returns the Tpetra::Map object associated with the range of this operator.
  const Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > getRangeMap() const;

  //! Returns in Y the result of a Xpetra::Operator applied to a Xpetra::MultiVector X.
  /*!
    \param[in]  X - Xpetra::MultiVector of dimension NumVectors to multiply with matrix.
    \param[out] Y - Xpetra::MultiVector of dimension NumVectors containing result.
  */
  void apply(const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& X,
             Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Y,
             Teuchos::ETransp mode = Teuchos::NO_TRANS,
             Scalar /* alpha */    = Teuchos::ScalarTraits<Scalar>::one(),
             Scalar /* beta */     = Teuchos::ScalarTraits<Scalar>::one()) const;

  //! Indicates whether this operator supports applying the adjoint operator.
  bool hasTransposeApply() const;

  //! Compute a residual R = B - (*this) * X
  void residual(const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& X,
                const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& B,
                Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& R) const;

  //! @name MueLu specific
  //@{

  //! Direct access to the underlying MueLu::Hierarchy.
  RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node> > GetHierarchy() const;

  //@}

 private:
  RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Hierarchy_;
};

}  // namespace MueLu

#endif  // MUELU_XPETRAOPERATOR_DECL_HPP
