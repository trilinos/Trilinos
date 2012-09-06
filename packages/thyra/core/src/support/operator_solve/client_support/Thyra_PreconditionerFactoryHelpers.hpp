// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roscoe A. Bartlett (bartlettra@ornl.gov) 
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


#ifndef THYRA_HIDE_DEPRECATED_CODE
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


#endif // THYRA_HIDE_DEPRECATED_CODE
} // namespace Thyra


#endif // THYRA_PRECONDITIONER_FACTORY_HELPERS_DECL_HPP
