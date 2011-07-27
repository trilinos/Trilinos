// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
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


//
// Deprecated
//


/** \brief Deprecated.
 *
 * \defgroup Thyra_deprecated_grp
 */
template <class Scalar>
THYRA_DEPRECATED
void initializePrec(
  const PreconditionerFactoryBase<Scalar> &precFactory,
  const Teuchos::RCP<const LinearOpBase<Scalar> > &fwdOp,
  PreconditionerBase<Scalar> *prec,
  const ESupportSolveUse supportSolveUse = SUPPORT_SOLVE_UNSPECIFIED
  )
{
  initializePrec<Scalar>(precFactory, fwdOp, Teuchos::ptr(prec),
    supportSolveUse);
}


/** \brief Deprecated.
 *
 * \defgroup Thyra_deprecated_grp
 */
template <class Scalar>
THYRA_DEPRECATED
void uninitializePrec(
  const PreconditionerFactoryBase<Scalar> &precFactory,
  PreconditionerBase<Scalar> *prec,
  Teuchos::RCP<const LinearOpBase<Scalar> > *fwdOp = NULL,
  ESupportSolveUse *supportSolveUse = NULL
  )
{
  using Teuchos::ptr;
  uninitializePrec<Scalar>(precFactory, ptr(prec), ptr(fwdOp),
    ptr(supportSolveUse));
}


} // namespace Thyra


#endif // THYRA_PRECONDITIONER_FACTORY_HELPERS_DECL_HPP
