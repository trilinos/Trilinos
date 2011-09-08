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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_DEFAULT_SERIAL_DENSE_LINEAR_OP_WITH_SOLVE_FACTORY_HPP
#define THYRA_DEFAULT_SERIAL_DENSE_LINEAR_OP_WITH_SOLVE_FACTORY_HPP


#include "Thyra_DefaultSerialDenseLinearOpWithSolveFactory_decl.hpp"
#include "Thyra_DefaultSerialDenseLinearOpWithSolve.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_ScaledAdjointLinearOpBase.hpp"
#include "Thyra_DefaultLinearOpSource.hpp"


namespace Thyra {


// Overridden from ParameterListAcceptor


template<class Scalar>
void DefaultSerialDenseLinearOpWithSolveFactory<Scalar>::setParameterList(
  RCP<ParameterList> const& paramList
  )
{
  paramList->validateParameters(*this->getValidParameters());
  // Nothing to set because we have not parameters!
}


template<class Scalar>
RCP<const ParameterList>
DefaultSerialDenseLinearOpWithSolveFactory<Scalar>::getValidParameters() const
{
  static RCP<const ParameterList> validPL = Teuchos::parameterList();
  return validPL;
}


// Overridden from LinearOpWithSolveFactoyBase


template<class Scalar>
bool DefaultSerialDenseLinearOpWithSolveFactory<Scalar>::acceptsPreconditionerFactory() const
{
  return false;
}


template<class Scalar>
void DefaultSerialDenseLinearOpWithSolveFactory<Scalar>::setPreconditionerFactory(
  const RCP<PreconditionerFactoryBase<Scalar> > &precFactory,
  const std::string &precFactoryName
  )
{
  TEST_FOR_EXCEPT_MSG(true, "Error, we don't support a preconditioner factory!");
}


template<class Scalar>
RCP<PreconditionerFactoryBase<Scalar> >
DefaultSerialDenseLinearOpWithSolveFactory<Scalar>::getPreconditionerFactory() const
{
  return Teuchos::null;
}


template<class Scalar>
void DefaultSerialDenseLinearOpWithSolveFactory<Scalar>::unsetPreconditionerFactory(
  RCP<PreconditionerFactoryBase<Scalar> > *precFactory,
  std::string *precFactoryName
  )
{
  TEST_FOR_EXCEPT_MSG(true, "Error, we don't support a preconditioner factory!");
}


template<class Scalar>
bool DefaultSerialDenseLinearOpWithSolveFactory<Scalar>::isCompatible(
  const LinearOpSourceBase<Scalar> &fwdOpSrc
  ) const
{
  return !is_null(
    Teuchos::rcp_dynamic_cast<const MultiVectorBase<Scalar> >(fwdOpSrc.getOp()));
}


template<class Scalar>
RCP<LinearOpWithSolveBase<Scalar> >
DefaultSerialDenseLinearOpWithSolveFactory<Scalar>::createOp() const
{
  return defaultSerialDenseLinearOpWithSolve<Scalar>();
}


template<class Scalar>
void DefaultSerialDenseLinearOpWithSolveFactory<Scalar>::initializeOp(
  const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
  LinearOpWithSolveBase<Scalar> *Op,
  const ESupportSolveUse supportSolveUse
  ) const
{

  using Teuchos::dyn_cast;
  using Teuchos::rcp_dynamic_cast;

#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(0==Op);
#endif

  const RCP<const LinearOpBase<Scalar> > tmpFwdOp = fwdOpSrc->getOp();
  RCP<const LinearOpBase<Scalar> > fwdOp;
  Scalar fwdOp_scalar = 0.0;
  EOpTransp fwdOp_transp;
  unwrap<Scalar>(tmpFwdOp, &fwdOp_scalar, &fwdOp_transp, &fwdOp);

  const RCP<const MultiVectorBase<Scalar> > fwdMv =
    rcp_dynamic_cast<const MultiVectorBase<Scalar> >(fwdOp, true);

  dyn_cast<DefaultSerialDenseLinearOpWithSolve<Scalar> >(*Op).initialize(fwdMv);

}


template<class Scalar>
void DefaultSerialDenseLinearOpWithSolveFactory<Scalar>::initializeAndReuseOp(
  const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
  LinearOpWithSolveBase<Scalar> *Op
  ) const
{
  initializeOp(fwdOpSrc, Op, SUPPORT_SOLVE_UNSPECIFIED);
}


template<class Scalar>
void DefaultSerialDenseLinearOpWithSolveFactory<Scalar>::uninitializeOp(
  LinearOpWithSolveBase<Scalar> *Op,
  RCP<const LinearOpSourceBase<Scalar> > *fwdOpSrc,
  RCP<const PreconditionerBase<Scalar> > *prec,
  RCP<const LinearOpSourceBase<Scalar> > *approxFwdOpSrc,
  ESupportSolveUse *supportSolveUse
  ) const
{
  using Teuchos::dyn_cast;
  using Teuchos::is_null;
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(0==Op);
#endif // TEUCHOS_DEBUG
  typedef DefaultSerialDenseLinearOpWithSolve<Scalar> DSDLOWS;
  DSDLOWS &dsdlows = dyn_cast<DSDLOWS>(*Op);
  if (fwdOpSrc) {
    // find a valid fwdOp
    const RCP<const LinearOpBase<Scalar> > fwdOp = dsdlows.getFwdOp();
    // pass out a valid fwsOpSrc
    if (!is_null(fwdOp)) {
      *fwdOpSrc = defaultLinearOpSource<Scalar>(fwdOp);
    } else {
      *fwdOpSrc = Teuchos::null;
    }
  }
  if (prec) *prec = Teuchos::null;
  if (approxFwdOpSrc) *approxFwdOpSrc = Teuchos::null;
}


template<class Scalar>
bool DefaultSerialDenseLinearOpWithSolveFactory<Scalar>::supportsPreconditionerInputType(
  const EPreconditionerInputType precOpType
  ) const
{
  // LAPACK does not support any external preconditioners!
  return false;
}


template<class Scalar>
void DefaultSerialDenseLinearOpWithSolveFactory<Scalar>::initializePreconditionedOp(
  const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
  const RCP<const PreconditionerBase<Scalar> > &prec,
  LinearOpWithSolveBase<Scalar> *Op,
  const ESupportSolveUse supportSolveUse
  ) const
{
  TEST_FOR_EXCEPT_MSG(true, "Error, we don't support an external preconditioner!");
}


template<class Scalar>
void DefaultSerialDenseLinearOpWithSolveFactory<Scalar>::initializeApproxPreconditionedOp(
  const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
  const RCP<const LinearOpSourceBase<Scalar> > &approxFwdOpSrc,
  LinearOpWithSolveBase<Scalar> *Op,
  const ESupportSolveUse supportSolveUse
  ) const
{
  TEST_FOR_EXCEPT_MSG(true, "Error, we don't support an external preconditioner!");
}


} // namespace Thyra


#endif // THYRA_DEFAULT_SERIAL_DENSE_LINEAR_OP_WITH_SOLVE_FACTORY_HPP
