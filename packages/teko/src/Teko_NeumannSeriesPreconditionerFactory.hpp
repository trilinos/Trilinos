// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teko_NeumannSeriesPreconditionerFactory_hpp__
#define __Teko_NeumannSeriesPreconditionerFactory_hpp__

#include "Teko_NeumannSeriesPreconditionerFactoryDecl.hpp"

#include "Thyra_DefaultPreconditioner.hpp"
#include "Thyra_DefaultPreconditioner.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"
#include "Thyra_DefaultAddedLinearOp.hpp"
#include "Thyra_DefaultMultipliedLinearOp.hpp"
#include "Thyra_DefaultIdentityLinearOp.hpp"

#include "Teuchos_Array.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"

namespace Teko {

using Teuchos::RCP;
using Teuchos::rcp;

static RCP<Teuchos::StringToIntegralParameterEntryValidator<Teko::DiagonalType> > scalingTypeVdtor;

template <typename ScalarT>
NeumannSeriesPreconditionerFactory<ScalarT>::NeumannSeriesPreconditionerFactory()
    : numberOfTerms_(1), scalingType_(Teko::NotDiag) {}

//! is this operator compatiable with the preconditioner factory?
template <typename ScalarT>
bool NeumannSeriesPreconditionerFactory<ScalarT>::isCompatible(
    const Thyra::LinearOpSourceBase<ScalarT> & /* fwdOpSrc */) const {
  return true;
}

//! create an instance of the preconditioner
template <typename ScalarT>
RCP<Thyra::PreconditionerBase<ScalarT> > NeumannSeriesPreconditionerFactory<ScalarT>::createPrec()
    const {
  return rcp(new Thyra::DefaultPreconditioner<ScalarT>());
}

/** \brief initialize a newly created preconditioner object
 *
 * Initialize a newly created preconditioner object.
 *
 * \param[in] fwdOpSrc Forward operator to be preconditioned
 * \param[in,out] precOp Return location for the preconditioner
 * \param[in] supportSolveUse Thyra information (?)
 */
template <typename ScalarT>
void NeumannSeriesPreconditionerFactory<ScalarT>::initializePrec(
    const RCP<const Thyra::LinearOpSourceBase<ScalarT> > &fwdOpSrc,
    Thyra::PreconditionerBase<ScalarT> *prec,
    const Thyra::ESupportSolveUse /* supportSolveUse */) const {
  using Thyra::add;
  using Thyra::multiply;
  using Thyra::scale;

  RCP<const Thyra::LinearOpBase<ScalarT> > M;  // left-preconditioner
  RCP<const Thyra::LinearOpBase<ScalarT> > A = fwdOpSrc->getOp();
  if (scalingType_ != Teko::NotDiag) {
    M = Teko::getInvDiagonalOp(A, scalingType_);
    A = Thyra::multiply(M, A);
  }

  RCP<const Thyra::LinearOpBase<ScalarT> > id   = Thyra::identity<ScalarT>(A->range());  // I
  RCP<const Thyra::LinearOpBase<ScalarT> > idMA = add(id, scale(-1.0, A));               // I - A

  RCP<const Thyra::LinearOpBase<ScalarT> > precOp;
  if (numberOfTerms_ == 1) {
    // no terms requested, just return identity matrix
    precOp = id;
  } else {
    int iters = numberOfTerms_ - 1;
    // use Horner's method to computed higher order polynomials
    precOp = add(scale(2.0, id), scale(-1.0, A));                              // I + (I - A)
    for (int i = 0; i < iters; i++) precOp = add(id, multiply(idMA, precOp));  // I + (I - A) * p
  }

  // multiply by the preconditioner if it exists
  if (M != Teuchos::null) precOp = Thyra::multiply(precOp, M);

  // must first cast that to be initialized
  Thyra::DefaultPreconditioner<ScalarT> &dPrec =
      Teuchos::dyn_cast<Thyra::DefaultPreconditioner<ScalarT> >(*prec);

  // this const-cast is unfortunate...needs to be fixed (larger than it seems!) ECC FIXME!
  dPrec.initializeUnspecified(Teuchos::rcp_const_cast<Thyra::LinearOpBase<ScalarT> >(precOp));
}

//! wipe clean a already initialized preconditioner object
template <typename ScalarT>
void NeumannSeriesPreconditionerFactory<ScalarT>::uninitializePrec(
    Thyra::PreconditionerBase<ScalarT> *prec,
    RCP<const Thyra::LinearOpSourceBase<ScalarT> > * /* fwdOpSrc */,
    Thyra::ESupportSolveUse * /* supportSolveUse */) const {
  Thyra::DefaultPreconditioner<ScalarT> &dPrec =
      Teuchos::dyn_cast<Thyra::DefaultPreconditioner<ScalarT> >(*prec);

  // wipe out any old preconditioner
  dPrec.uninitialize();
}

// for ParameterListAcceptor

//! \brief Set parameters from a parameter list
template <typename ScalarT>
void NeumannSeriesPreconditionerFactory<ScalarT>::setParameterList(
    const RCP<Teuchos::ParameterList> &paramList) {
  TEUCHOS_TEST_FOR_EXCEPT(paramList == Teuchos::null);

  // check the parameters are correct
  paramList->validateParametersAndSetDefaults(*getValidParameters(), 0);

  // store the parameter list
  paramList_ = paramList;

  numberOfTerms_ = paramList_->get<int>("Number of Terms");

  // get the scaling type
  scalingType_                         = Teko::NotDiag;
  const Teuchos::ParameterEntry *entry = paramList_->getEntryPtr("Scaling Type");
  if (entry != NULL) scalingType_ = scalingTypeVdtor->getIntegralValue(*entry);
}

/** \brief Get the valid parameters */
template <typename ScalarT>
RCP<const Teuchos::ParameterList> NeumannSeriesPreconditionerFactory<ScalarT>::getValidParameters()
    const {
  static RCP<Teuchos::ParameterList> validPL;

  // only fill valid parameter list once
  if (validPL == Teuchos::null) {
    RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList());

    // build the validator for scaling type
    scalingTypeVdtor = Teuchos::stringToIntegralParameterEntryValidator<DiagonalType>(
        Teuchos::tuple<std::string>("Diagonal", "Lumped", "AbsRowSum", "None"),
        Teuchos::tuple<Teko::DiagonalType>(Teko::Diagonal, Teko::Lumped, Teko::AbsRowSum,
                                           Teko::NotDiag),
        "Scaling Type");

    pl->set<int>("Number of Terms", 1,
                 "The number of terms to use in the Neumann series expansion.");
    pl->set("Scaling Type", "None", "The number of terms to use in the Neumann series expansion.",
            scalingTypeVdtor);

    validPL = pl;
  }

  return validPL;
}

//! Unset the parameter list that was set using setParameterList().
template <typename ScalarT>
RCP<Teuchos::ParameterList> NeumannSeriesPreconditionerFactory<ScalarT>::unsetParameterList() {
  Teuchos::RCP<Teuchos::ParameterList> oldList = paramList_;
  paramList_                                   = Teuchos::null;
  return oldList;
}

//! Get the parameter list that was set using setParameterList().
template <typename ScalarT>
RCP<const Teuchos::ParameterList> NeumannSeriesPreconditionerFactory<ScalarT>::getParameterList()
    const {
  return paramList_;
}

//! Get the parameter list that was set using setParameterList().
template <typename ScalarT>
RCP<Teuchos::ParameterList>
NeumannSeriesPreconditionerFactory<ScalarT>::getNonconstParameterList() {
  return paramList_;
}

template <typename ScalarT>
std::string NeumannSeriesPreconditionerFactory<ScalarT>::description() const {
  std::ostringstream oss;
  oss << "Teko::NeumannSeriesPreconditionerFactory";
  return oss.str();
}

}  // end namespace Teko

#endif
