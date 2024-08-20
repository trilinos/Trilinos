// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_TPETRAOPERATOR_DECL_HPP
#define MUELU_TPETRAOPERATOR_DECL_HPP

#include "MueLu_ConfigDefs.hpp"

#include <Tpetra_Operator.hpp>
#include <Tpetra_MultiVector_decl.hpp>
#include "MueLu_Level.hpp"
#include "MueLu_Hierarchy_decl.hpp"

namespace MueLu {

/*!  @brief Wraps an existing MueLu::Hierarchy as a Tpetra::Operator.
 */
template <class Scalar        = Tpetra::Operator<>::scalar_type,
          class LocalOrdinal  = typename Tpetra::Operator<Scalar>::local_ordinal_type,
          class GlobalOrdinal = typename Tpetra::Operator<Scalar, LocalOrdinal>::global_ordinal_type,
          class Node          = typename Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal>::node_type>
class TpetraOperator : public Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
 protected:
  TpetraOperator() = delete;

 public:
  //! @name Constructor/Destructor
  //@{

  //! Constructor
  TpetraOperator(const RCP<Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& Op)
    : Operator_(Op) {}

  //! Constructor
  TpetraOperator(const RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& H)
    : Hierarchy_(H) {}

  //! Destructor.
  virtual ~TpetraOperator();

  //@}

  //! Returns the Tpetra::Map object associated with the domain of this operator.
  Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > getDomainMap() const;

  //! Returns the Tpetra::Map object associated with the range of this operator.
  Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > getRangeMap() const;

  //! Returns in Y the result of a Tpetra::Operator applied to a Tpetra::MultiVector X.
  /*!
    \param[in]  X - Tpetra::MultiVector of dimension NumVectors to multiply with matrix.
    \param[out] Y -Tpetra::MultiVector of dimension NumVectors containing result.
  */
  void apply(const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& X,
             Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Y,
             Teuchos::ETransp mode = Teuchos::NO_TRANS,
             Scalar alpha          = Teuchos::ScalarTraits<Scalar>::one(),
             Scalar beta           = Teuchos::ScalarTraits<Scalar>::one()) const;

  //! Indicates whether this operator supports applying the adjoint operator.
  bool hasTransposeApply() const;

  //! @name MueLu specific
  //@{

  //! Direct access to the underlying MueLu::Hierarchy.
  RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node> > GetHierarchy() const;

  //! Direct access to the underlying MueLu::Operator
  RCP<Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> > GetOperator() const;

  //@}

 private:
  RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Hierarchy_;
  RCP<Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Operator_;
};

}  // namespace MueLu

#endif  // MUELU_TPETRAOPERATOR_DECL_HPP
