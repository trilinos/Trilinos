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

#ifndef THYRA_LINEAR_OP_WITH_SOLVE_FACTORY_HELPERS_HPP
#define THYRA_LINEAR_OP_WITH_SOLVE_FACTORY_HELPERS_HPP

#include "Thyra_LinearOpWithSolveFactoryBase.hpp"
#include "Thyra_DefaultLinearOpSource.hpp"
#include "Thyra_DefaultPreconditioner.hpp"
#include "Thyra_DefaultInverseLinearOp.hpp"


namespace Thyra {


/** \brief Return if the forward operator is a compatible source for a LOWSFB
 * object.
 *
 * \relates LinearOpWithSolveFactoryBase
 */
template<class Scalar>
bool isCompatible(
  const LinearOpWithSolveFactoryBase<Scalar> &lowsFactory,
  const LinearOpBase<Scalar> &fwdOp
  )
{
  return lowsFactory.isCompatible(*defaultLinearOpSource(fwdOp));
}


/** \brief Set default label on a LOWSB object.
 *
 * \relates LinearOpWithSolveFactoryBase
 */
template<class Scalar>
void setDefaultObjectLabel(
  const LinearOpBase<Scalar> &fwdOp,
  const Ptr<LinearOpWithSolveBase<Scalar> > &Op
  )
{
  const std::string OpLabel = Op->getObjectLabel();
  const std::string fwdOpLabel = fwdOp.getObjectLabel();
  if ( !OpLabel.length() && fwdOpLabel.length() )
    Op->setObjectLabel(fwdOpLabel);
}


/** \brief Initialize a pre-created LOWSB object given a forward operator.
 *
 * \relates LinearOpWithSolveFactoryBase
 */
template<class Scalar>
void initializeOp(
  const LinearOpWithSolveFactoryBase<Scalar> &lowsFactory,
  const RCP<const LinearOpBase<Scalar> > &fwdOp,
  const Ptr<LinearOpWithSolveBase<Scalar> > &Op,
  const ESupportSolveUse supportSolveUse = SUPPORT_SOLVE_UNSPECIFIED
  )
{
  lowsFactory.initializeOp(defaultLinearOpSource(fwdOp), &*Op, supportSolveUse);
  setDefaultObjectLabel(*fwdOp, Op);
}


/** \brief Reinitialize a pre-created LOWSB object given a forward operator,
 * reusing a much as possible from the prior LOWSB object.
 *
 * \relates LinearOpWithSolveFactoryBase
 */
template<class Scalar>
void initializeAndReuseOp(
  const LinearOpWithSolveFactoryBase<Scalar> &lowsFactory,
  const RCP<const LinearOpBase<Scalar> > &fwdOp,
  const Ptr<LinearOpWithSolveBase<Scalar> > &Op
  )
{
  lowsFactory.initializeAndReuseOp(defaultLinearOpSource(fwdOp), &*Op);
  setDefaultObjectLabel(*fwdOp, Op);
}


/** \brief Initialize a preconditioned LOWSB object given an external
 * preconditioner.
 *
 * \relates LinearOpWithSolveFactoryBase
 */
template<class Scalar>
void initializePreconditionedOp(
  const LinearOpWithSolveFactoryBase<Scalar> &lowsFactory,
  const RCP<const LinearOpBase<Scalar> > &fwdOp,
  const RCP<const PreconditionerBase<Scalar> > &prec,
  const Ptr<LinearOpWithSolveBase<Scalar> > &Op,
  const ESupportSolveUse supportSolveUse = SUPPORT_SOLVE_UNSPECIFIED
  )
{
  lowsFactory.initializePreconditionedOp(defaultLinearOpSource(fwdOp),
    prec, &*Op, supportSolveUse);
  setDefaultObjectLabel(*fwdOp ,Op);
}


/** \brief Initialize a preconditioned LOWSB object given an external operator
 * to be used to generate the preconditioner internally.
 *
 * \relates LinearOpWithSolveFactoryBase
 */
template<class Scalar>
void initializeApproxPreconditionedOp(
  const LinearOpWithSolveFactoryBase<Scalar> &lowsFactory,
  const RCP<const LinearOpBase<Scalar> > &fwdOp,
  const RCP<const LinearOpBase<Scalar> > &approxFwdOp,
  const Ptr<LinearOpWithSolveBase<Scalar> > &Op,
  const ESupportSolveUse supportSolveUse = SUPPORT_SOLVE_UNSPECIFIED
  )
{
  lowsFactory.initializeApproxPreconditionedOp(defaultLinearOpSource(fwdOp),
    defaultLinearOpSource(approxFwdOp), &*Op, supportSolveUse);
  setDefaultObjectLabel(*fwdOp,Op);
}


/** \brief Create and initialize a <tt>LinearOpWithSolveBase</tt> object from
 * a <tt>LinearOpBase</tt> object using a
 * <tt>LinearOpWithSolveFactoryBase</tt> strategy object.
 * \relates LinearOpWithSolveFactoryBase */
template<class Scalar>
RCP<LinearOpWithSolveBase<Scalar> >
linearOpWithSolve(
  const LinearOpWithSolveFactoryBase<Scalar> &lowsFactory,
  const RCP<const LinearOpBase<Scalar> > &fwdOp,
  const ESupportSolveUse supportSolveUse = SUPPORT_SOLVE_UNSPECIFIED
  )
{
  RCP<LinearOpWithSolveBase<Scalar> > Op = lowsFactory.createOp();
  Thyra::initializeOp<Scalar>( lowsFactory, fwdOp, Op.ptr(), supportSolveUse);
  return Op;
}


/** \brief Form a const implicit inverse operator <tt>M = inv(A)</tt> given a
 * factory.
 *
 * \relates LinearOpWithSolveFactoryBase
 */
template<class Scalar>
RCP<LinearOpBase<Scalar> >
inverse(
  const LinearOpWithSolveFactoryBase<Scalar> &LOWSF,
  const RCP<const LinearOpBase<Scalar> > &fwdOp,
  const ESupportSolveUse supportSolveUse = SUPPORT_SOLVE_UNSPECIFIED,
  const Ptr<const SolveCriteria<Scalar> > &fwdSolveCriteria = Teuchos::null,
  const EThrowOnSolveFailure throwOnFwdSolveFailure = THROW_ON_SOLVE_FAILURE,
  const Ptr<const SolveCriteria<Scalar> > &adjSolveCriteria = Teuchos::null,
  const EThrowOnSolveFailure throwOnAdjSolveFailure = THROW_ON_SOLVE_FAILURE
  )
{
  return inverse<Scalar>(linearOpWithSolve<Scalar>(LOWSF, fwdOp, supportSolveUse),
    fwdSolveCriteria, throwOnFwdSolveFailure, adjSolveCriteria, throwOnAdjSolveFailure);
}

  
/** \brief Uninitialized a pre-created LOWSB object, returning input objects
 * used to initialize it.
 *
 * \relates LinearOpWithSolveFactoryBase
 */
template<class Scalar>
void uninitializeOp(
  const LinearOpWithSolveFactoryBase<Scalar> &lowsFactory,
  const Ptr<LinearOpWithSolveBase<Scalar> > &Op,
  const Ptr<RCP<const LinearOpBase<Scalar> > > &fwdOp = Teuchos::null,
  const Ptr<RCP<const PreconditionerBase<Scalar> > > &prec = Teuchos::null,
  const Ptr<RCP<const LinearOpBase<Scalar> > > &approxFwdOp = Teuchos::null,
  const Ptr<ESupportSolveUse> &supportSolveUse = Teuchos::null
  )
{
  RCP<const LinearOpSourceBase<Scalar> > fwdOpSrc;
  RCP<const LinearOpSourceBase<Scalar> > approxFwdOpSrc;
  lowsFactory.uninitializeOp(Op.get(), &fwdOpSrc, prec.get(), &approxFwdOpSrc,
    supportSolveUse.get());
  if (nonnull(fwdOp)) {
    *fwdOp = ( nonnull(fwdOpSrc) ? fwdOpSrc->getOp() : Teuchos::null );
  }
  if (nonnull(approxFwdOp)) {
    *approxFwdOp = ( nonnull(approxFwdOpSrc) ? approxFwdOpSrc->getOp() : Teuchos::null );
  }
}


#ifndef THYRA_HIDE_DEPRECATED_CODE
//
// Deprecated
//


/** \brief Deprecated.
 *
 * \defgroup Thyra_deprecated_grp
 */
template<class Scalar>
THYRA_DEPRECATED
void setDefaultObjectLabel(
  const LinearOpBase<Scalar> &fwdOp,
  LinearOpWithSolveBase<Scalar> *Op
  )
{
  setDefaultObjectLabel<Scalar>(fwdOp, Teuchos::ptr(Op));
}


/** \brief Deprecated.
 *
 * \defgroup Thyra_deprecated_grp
 */
template<class Scalar>
THYRA_DEPRECATED
void initializeOp(
  const LinearOpWithSolveFactoryBase<Scalar> &lowsFactory,
  const RCP<const LinearOpBase<Scalar> > &fwdOp,
  LinearOpWithSolveBase<Scalar> *Op,
  const ESupportSolveUse supportSolveUse = SUPPORT_SOLVE_UNSPECIFIED
  )
{
  initializeOp<Scalar>(lowsFactory, fwdOp, Teuchos::ptr(Op), supportSolveUse);
}


/** \brief Deprecated.
 *
 * \defgroup Thyra_deprecated_grp
 */
template<class Scalar>
THYRA_DEPRECATED
void initializeAndReuseOp(
  const LinearOpWithSolveFactoryBase<Scalar> &lowsFactory,
  const RCP<const LinearOpBase<Scalar> > &fwdOp,
  LinearOpWithSolveBase<Scalar> *Op
  )
{
  initializeAndReuseOp<Scalar>(lowsFactory, fwdOp, Teuchos::ptr(Op));
}


/** \brief Deprecated.
 *
 * \defgroup Thyra_deprecated_grp
 */
template<class Scalar>
THYRA_DEPRECATED
void initializePreconditionedOp(
  const LinearOpWithSolveFactoryBase<Scalar> &lowsFactory,
  const RCP<const LinearOpBase<Scalar> > &fwdOp,
  const RCP<const PreconditionerBase<Scalar> > &prec,
  LinearOpWithSolveBase<Scalar> *Op,
  const ESupportSolveUse supportSolveUse = SUPPORT_SOLVE_UNSPECIFIED
  )
{
  initializePreconditionedOp<Scalar>(lowsFactory, fwdOp, prec, Teuchos::ptr(Op),
    supportSolveUse);
}


/** \brief Deprecated.
 *
 * \defgroup Thyra_deprecated_grp
 */
template<class Scalar>
THYRA_DEPRECATED
void initializeApproxPreconditionedOp(
  const LinearOpWithSolveFactoryBase<Scalar> &lowsFactory,
  const RCP<const LinearOpBase<Scalar> > &fwdOp,
  const RCP<const LinearOpBase<Scalar> > &approxFwdOp,
  LinearOpWithSolveBase<Scalar> *Op,
  const ESupportSolveUse supportSolveUse = SUPPORT_SOLVE_UNSPECIFIED
  )
{
  initializeApproxPreconditionedOp<Scalar>(lowsFactory, fwdOp, approxFwdOp,
    Teuchos::ptr(Op), supportSolveUse);
}


/** \brief Deprecated.
 *
 * \defgroup Thyra_deprecated_grp
 */
template<class Scalar>
THYRA_DEPRECATED
RCP<LinearOpBase<Scalar> >
inverse(
  const LinearOpWithSolveFactoryBase<Scalar> &LOWSF,
  const RCP<const LinearOpBase<Scalar> > &fwdOp,
  const ESupportSolveUse supportSolveUse,
  const SolveCriteria<Scalar> *fwdSolveCriteria = NULL,
  const EThrowOnSolveFailure throwOnFwdSolveFailure = THROW_ON_SOLVE_FAILURE,
  const SolveCriteria<Scalar> *adjSolveCriteria = NULL,
  const EThrowOnSolveFailure throwOnAdjSolveFailure = THROW_ON_SOLVE_FAILURE
  )
{
  return inverse<Scalar>(LOWSF, fwdOp, supportSolveUse,
    Teuchos::ptr(fwdSolveCriteria), throwOnFwdSolveFailure,
    Teuchos::ptr(adjSolveCriteria), throwOnAdjSolveFailure);
}

  
/** \brief Deprecated.
 *
 * \defgroup Thyra_deprecated_grp
 */
template<class Scalar>
THYRA_DEPRECATED
void uninitializeOp(
  const LinearOpWithSolveFactoryBase<Scalar> &lowsFactory,
  LinearOpWithSolveBase<Scalar> *Op,
  RCP<const LinearOpBase<Scalar> > *fwdOp = NULL,
  RCP<const PreconditionerBase<Scalar> > *prec = NULL,
  RCP<const LinearOpBase<Scalar> > *approxFwdOp = NULL,
  ESupportSolveUse *supportSolveUse = NULL
  )
{
  using Teuchos::ptr;
  uninitializeOp<Scalar>(lowsFactory, ptr(Op), ptr(fwdOp), ptr(prec),
    ptr(approxFwdOp), ptr(supportSolveUse));
}


#endif // THYRA_HIDE_DEPRECATED_CODE
} // namespace Thyra


#endif // THYRA_LINEAR_OP_WITH_SOLVE_FACTORY_HELPERS_HPP
