// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teko_PreconditionerLinearOp_hpp__
#define __Teko_PreconditionerLinearOp_hpp__

#include "Teko_PreconditionerLinearOpDecl.hpp"

#include "Thyra_LinearOpBase.hpp"
#include "Thyra_PreconditionerBase.hpp"

namespace Teko {

template <typename ScalarT>
PreconditionerLinearOp<ScalarT>::PreconditionerLinearOp() {}

template <typename ScalarT>
PreconditionerLinearOp<ScalarT>::PreconditionerLinearOp(
    const Teuchos::RCP<Thyra::PreconditionerBase<ScalarT> >& prec) {
  preconditioner_.initialize(prec);
}

template <typename ScalarT>
PreconditionerLinearOp<ScalarT>::PreconditionerLinearOp(
    const Teuchos::RCP<const Thyra::PreconditionerBase<ScalarT> >& prec) {
  preconditioner_.initialize(prec);
}

//! build a linear operator using this preconditioner, this initialization permits changes
template <typename ScalarT>
void PreconditionerLinearOp<ScalarT>::initialize(
    const Teuchos::RCP<Thyra::PreconditionerBase<ScalarT> >& prec) {
  uninitialize();
  preconditioner_.initialize(prec);
}

//! build a linear operator using this preconditioner, this initialization refuses changes
template <typename ScalarT>
void PreconditionerLinearOp<ScalarT>::initialize(
    const Teuchos::RCP<const Thyra::PreconditionerBase<ScalarT> >& prec) {
  uninitialize();
  preconditioner_.initialize(prec);
}

//! Disassociate this object with the currently owned preconditioner
template <typename ScalarT>
void PreconditionerLinearOp<ScalarT>::uninitialize() {
  preconditioner_.uninitialize();
}

/** @brief Range space of this operator */
template <typename ScalarT>
Teuchos::RCP<const Thyra::VectorSpaceBase<ScalarT> > PreconditionerLinearOp<ScalarT>::range()
    const {
  return getOperator_cnoc()->range();
}

/** @brief Domain space of this operator */
template <typename ScalarT>
Teuchos::RCP<const Thyra::VectorSpaceBase<ScalarT> > PreconditionerLinearOp<ScalarT>::domain()
    const {
  return getOperator_cnoc()->domain();
}

template <typename ScalarT>
bool PreconditionerLinearOp<ScalarT>::opSupportedImpl(const Thyra::EOpTransp M_trans) const {
  return getOperator_cnoc()->opSupported(M_trans);
}

template <typename ScalarT>
void PreconditionerLinearOp<ScalarT>::applyImpl(
    const Thyra::EOpTransp M_trans, const Thyra::MultiVectorBase<ScalarT>& x,
    const Teuchos::Ptr<Thyra::MultiVectorBase<ScalarT> >& y, const ScalarT alpha,
    const ScalarT beta) const {
  getOperator_cnoc()->apply(M_trans, x, y, alpha, beta);
}

//! Get a nonconstant <code>PreconditionerBase</code> object
template <typename ScalarT>
Teuchos::RCP<Thyra::PreconditionerBase<ScalarT> >
PreconditionerLinearOp<ScalarT>::getNonconstPreconditioner() {
  return preconditioner_.getNonconstObj();
}

//! Get a constant <code>PreconditionerBase</code> object
template <typename ScalarT>
Teuchos::RCP<const Thyra::PreconditionerBase<ScalarT> >
PreconditionerLinearOp<ScalarT>::getPreconditioner() const {
  return preconditioner_.getConstObj();
}

//! get operator associated with the preconditioner
template <typename ScalarT>
Teuchos::ConstNonconstObjectContainer<Thyra::LinearOpBase<ScalarT> >
PreconditionerLinearOp<ScalarT>::getOperator_cnoc() const {
  Teuchos::ConstNonconstObjectContainer<Thyra::LinearOpBase<ScalarT> > oper;
  oper.initialize(preconditioner_.getConstObj()->getUnspecifiedPrecOp());

  return oper;
}

//! get operator associated with the preconditioner
template <typename ScalarT>
Teuchos::ConstNonconstObjectContainer<Thyra::LinearOpBase<ScalarT> >
PreconditionerLinearOp<ScalarT>::getOperator_cnoc() {
  Teuchos::ConstNonconstObjectContainer<Thyra::LinearOpBase<ScalarT> > oper;
  oper.initialize(preconditioner_.getNonconstObj()->getNonconstUnspecifiedPrecOp());

  return oper;
}

template <typename ScalarT>
void PreconditionerLinearOp<ScalarT>::describe(Teuchos::FancyOStream& out_arg,
                                               const Teuchos::EVerbosityLevel verbLevel) const {
  using Teuchos::OSTab;

  Teuchos::RCP<Teuchos::FancyOStream> out = rcp(&out_arg, false);
  OSTab tab0(out);
  switch (verbLevel) {
    case Teuchos::VERB_DEFAULT:
    case Teuchos::VERB_LOW:
      *out << this->description() << " ( [Operator] = " << getOperator_cnoc()->description() << " )"
           << std::endl;
      break;
    case Teuchos::VERB_MEDIUM:
    case Teuchos::VERB_HIGH:
    case Teuchos::VERB_EXTREME: {
      *out << Teuchos::Describable::description() << "{"
           << "rangeDim=" << this->range()->dim() << ",domainDim=" << this->domain()->dim()
           << "}\n";
      {
        OSTab tab1(out);
        *out << "[Operator] = ";
        *out << Teuchos::describe(*getOperator_cnoc(), verbLevel);
      }
      break;
    }
    default: TEUCHOS_TEST_FOR_EXCEPT(true);  // Should never get here!
  }
}

inline Teko::LinearOp extractOperatorFromPrecOp(const Teko::LinearOp& lo) {
  Teuchos::RCP<const Teko::PreconditionerLinearOp<double> > plo =
      Teuchos::rcp_dynamic_cast<const Teko::PreconditionerLinearOp<double> >(lo);
  if (plo != Teuchos::null)
    return plo->getOperator();
  else
    return lo;
}

}  // end namespace Teko

#endif
