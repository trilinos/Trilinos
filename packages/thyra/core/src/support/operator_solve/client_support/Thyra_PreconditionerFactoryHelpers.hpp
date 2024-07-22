// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_PRECONDITIONER_FACTORY_HELPERS_DECL_HPP
#define THYRA_PRECONDITIONER_FACTORY_HELPERS_DECL_HPP


#include "Thyra_PreconditionerFactoryBase.hpp"
#include "Thyra_DefaultLinearOpSource.hpp"
#include "Thyra_DefaultPreconditioner.hpp"


namespace Thyra {


/** \brief Initialize a preconditioner from a forward linear operator.
 *
 * \relates PreconditionerFactoryBase
 */
template <class Scalar>
void initializePrec(
  const PreconditionerFactoryBase<Scalar> &precFactory,
  const Teuchos::RCP<const LinearOpBase<Scalar> > &fwdOp,
  const Teuchos::Ptr<PreconditionerBase<Scalar> > &prec,
  const ESupportSolveUse supportSolveUse = SUPPORT_SOLVE_UNSPECIFIED
  )
{
  precFactory.initializePrec(defaultLinearOpSource(fwdOp), prec.get(),
    supportSolveUse);
}


/** \brief Uninitialize a preconditioner and optionally extra what was used to
 * create it.
 *
 * \relates PreconditionerFactoryBase
 */
template <class Scalar>
void uninitializePrec(
  const PreconditionerFactoryBase<Scalar> &precFactory,
  const Teuchos::Ptr<PreconditionerBase<Scalar> > &prec,
  const Teuchos::Ptr<Teuchos::RCP<const LinearOpBase<Scalar> > > &fwdOp = Teuchos::null,
  const Teuchos::Ptr<ESupportSolveUse> &supportSolveUse = Teuchos::null
  )
{
  Teuchos::RCP<const LinearOpSourceBase<Scalar> > fwdOpSrc;
  precFactory.uninitializePrec(prec, Teuchos::outArg(fwdOpSrc), supportSolveUse);
  if (nonnull(fwdOp)) {*fwdOp = fwdOpSrc->getOp();}
}


/** \brief Create and initialize a preconditioner from a forward linear operator.
 *
 * \relates PreconditionerFactoryBase
 */
template <class Scalar>
Teuchos::RCP<PreconditionerBase<Scalar> >
prec(
  const PreconditionerFactoryBase<Scalar> &precFactory,
  const Teuchos::RCP<const LinearOpBase<Scalar> > &fwdOp,
  const ESupportSolveUse supportSolveUse = SUPPORT_SOLVE_UNSPECIFIED
  )
{
  Teuchos::RCP<PreconditionerBase<Scalar> > prec =
    precFactory.createPrec();
  precFactory.initializePrec(defaultLinearOpSource(fwdOp), &*prec, supportSolveUse);
  return prec;
}


} // namespace Thyra


#endif // THYRA_PRECONDITIONER_FACTORY_HELPERS_DECL_HPP
