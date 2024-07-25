// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_XPETRA_LINEAR_OP_DECL_HPP
#define THYRA_XPETRA_LINEAR_OP_DECL_HPP

#include "Thyra_LinearOpDefaultBase.hpp"
#include "Xpetra_Operator.hpp"
#include "Teuchos_ConstNonconstObjectContainer.hpp"

namespace Thyra {

/** \brief Concrete Thyra::LinearOpBase subclass for Xpetra::Operator.
 *
 * \todo Move this to Thyra??
 *
 * \ingroup Xpetra_Thyra_Op_Vec_adapters_grp
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal = LocalOrdinal,
          class Node = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
class XpetraLinearOp
  : virtual public Thyra::LinearOpDefaultBase<Scalar> {
 public:
  /** \name Constructors/initializers. */
  //@{

  /** \brief Construct to uninitialized. */
  XpetraLinearOp();

  ~XpetraLinearOp();

  /** \brief Initialize. */
  void initialize(
      const RCP<const VectorSpaceBase<Scalar> > &rangeSpace,
      const RCP<const VectorSpaceBase<Scalar> > &domainSpace,
      const RCP<Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> > &xpetraOperator);

  /** \brief Initialize. */
  void constInitialize(
      const RCP<const VectorSpaceBase<Scalar> > &rangeSpace,
      const RCP<const VectorSpaceBase<Scalar> > &domainSpace,
      const RCP<const Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> > &xpetraOperator);

  /** \brief Get embedded non-const Xpetra::Operator. */
  RCP<Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  getXpetraOperator();

  /** \brief Get embedded const Xpetra::Operator. */
  RCP<const Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  getConstXpetraOperator() const;

  //@}

  /** \name Public Overridden functions from LinearOpBase. */
  //@{

  /** \brief . */
  RCP<const Thyra::VectorSpaceBase<Scalar> > range() const;

  /** \brief . */
  RCP<const Thyra::VectorSpaceBase<Scalar> > domain() const;

  //@}

 protected:
  /** \name Protected Overridden functions from LinearOpBase. */
  //@{

  /** \brief . */
  bool opSupportedImpl(Thyra::EOpTransp M_trans) const;

  /** \brief . */
  void applyImpl(
      const Thyra::EOpTransp M_trans,
      const Thyra::MultiVectorBase<Scalar> &X_in,
      const Teuchos::Ptr<Thyra::MultiVectorBase<Scalar> > &Y_inout,
      const Scalar alpha,
      const Scalar beta) const;

  //@}

 private:
  RCP<const VectorSpaceBase<Scalar> >
      rangeSpace_;

  RCP<const VectorSpaceBase<Scalar> >
      domainSpace_;

  Teuchos::ConstNonconstObjectContainer<Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
      xpetraOperator_;

  template <class XpetraOperator_t>
  void initializeImpl(
      const RCP<const VectorSpaceBase<Scalar> > &rangeSpace,
      const RCP<const VectorSpaceBase<Scalar> > &domainSpace,
      const RCP<XpetraOperator_t> &xpetraOperator);
};

/** \brief Nonmmeber constructor for XpetraLinearOp.
 *
 * \relates XpetraLinearOp
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<XpetraLinearOp<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
xpetraLinearOp(
    const RCP<const VectorSpaceBase<Scalar> > &rangeSpace,
    const RCP<const VectorSpaceBase<Scalar> > &domainSpace,
    const RCP<Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> > &xpetraOperator) {
  const RCP<XpetraLinearOp<Scalar, LocalOrdinal, GlobalOrdinal, Node> > op =
      Teuchos::rcp(new XpetraLinearOp<Scalar, LocalOrdinal, GlobalOrdinal, Node>);
  op->initialize(rangeSpace, domainSpace, xpetraOperator);
  return op;
}

/** \brief Nonmmeber constructor for XpetraLinearOp.
 *
 * \relates XpetraLinearOp
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const XpetraLinearOp<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
constXpetraLinearOp(
    const RCP<const VectorSpaceBase<Scalar> > &rangeSpace,
    const RCP<const VectorSpaceBase<Scalar> > &domainSpace,
    const RCP<const Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> > &xpetraOperator) {
  const RCP<XpetraLinearOp<Scalar, LocalOrdinal, GlobalOrdinal, Node> > op =
      Teuchos::rcp(new XpetraLinearOp<Scalar, LocalOrdinal, GlobalOrdinal, Node>);
  op->constInitialize(rangeSpace, domainSpace, xpetraOperator);
  return op;
}

}  // namespace Thyra

#endif  // THYRA_XPETRA_LINEAR_OP_DECL_HPP
