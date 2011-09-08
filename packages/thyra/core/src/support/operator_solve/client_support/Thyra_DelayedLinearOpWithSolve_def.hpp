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


#ifndef THYRA_DELAYED_LINEAR_OP_WITH_SOLVE_HPP
#define THYRA_DELAYED_LINEAR_OP_WITH_SOLVE_HPP


#include "Thyra_DelayedLinearOpWithSolve_decl.hpp"
#include "Thyra_LinearOpWithSolveFactoryBase.hpp"
#include "Thyra_LinearOpWithSolveBase.hpp"
#include "Thyra_SolveSupportTypes.hpp"
#include "Teuchos_VerboseObject.hpp"


namespace Thyra {


// Constructor/Initializers


template <class Scalar>
DelayedLinearOpWithSolve<Scalar>::DelayedLinearOpWithSolve()
  :lows_is_valid_(false)
{}


template <class Scalar>
void DelayedLinearOpWithSolve<Scalar>::initialize(
  const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
  const RCP<const PreconditionerBase<Scalar> > &prec,
  const RCP<const LinearOpSourceBase<Scalar> > &approxFwdOpSrc,
  const ESupportSolveUse supportSolveUse,
  const RCP<LinearOpWithSolveFactoryBase<Scalar> > &lowsf
  )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(is_null(fwdOpSrc));
  TEST_FOR_EXCEPT(!is_null(prec) && !is_null(approxFwdOpSrc));
  TEST_FOR_EXCEPT(is_null(lowsf));
#endif  
  fwdOpSrc_ = fwdOpSrc;
  prec_ = prec;
  approxFwdOpSrc_ = approxFwdOpSrc;
  lowsf_ = lowsf;
  supportSolveUse_ = supportSolveUse;
  fwdOp_ = fwdOpSrc_->getOp().assert_not_null();
  lows_is_valid_ = false;
}


template <class Scalar>
RCP<const LinearOpSourceBase<Scalar> >
DelayedLinearOpWithSolve<Scalar>::getFwdOpSrc() const
{
  return fwdOpSrc_;
}


template <class Scalar>
RCP<const PreconditionerBase<Scalar> >
DelayedLinearOpWithSolve<Scalar>::getPrec() const
{
  return prec_;
}
 

template <class Scalar>
RCP<const LinearOpSourceBase<Scalar> >
DelayedLinearOpWithSolve<Scalar>::getApproxFwdOpSrc() const
{
  return approxFwdOpSrc_;
}


template <class Scalar>
ESupportSolveUse
DelayedLinearOpWithSolve<Scalar>::getSupportSolveUse() const
{
  return supportSolveUse_;
}


// Overridden from Teuchos::Describable


template<class Scalar>
std::string DelayedLinearOpWithSolve<Scalar>::description() const
{
  std::ostringstream oss;
  oss << this->Teuchos::Describable::description()
      << "{";
  oss << "fwdOp_=";
  if (!is_null(fwdOp_))
    oss << fwdOp_->description();
  else
    oss << "NULL,";
  oss << "lows=";
  if (!is_null(lows_))
    oss << lows_->description();
  else
    oss << "NULL";
  oss << "}";
  return oss.str();
}


// Overridden from LinearOpBase


template <class Scalar>
RCP< const VectorSpaceBase<Scalar> >
DelayedLinearOpWithSolve<Scalar>::range() const
{
  if (!is_null(fwdOp_))
    return fwdOp_->range();
  return Teuchos::null;
}


template <class Scalar>
RCP< const VectorSpaceBase<Scalar> >
DelayedLinearOpWithSolve<Scalar>::domain() const
{
  if (!is_null(fwdOp_))
    return fwdOp_->domain();
  return Teuchos::null;
}


template <class Scalar>
RCP<const LinearOpBase<Scalar> >
DelayedLinearOpWithSolve<Scalar>::clone() const
{
  return Teuchos::null; // ToDo: Implement if needed!
}


// protected


template<class Scalar>
void DelayedLinearOpWithSolve<Scalar>::informUpdatedVerbosityState() const
{
  if (!is_null(lowsf_)) {
    lowsf_->setVerbLevel(this->getVerbLevel());
    lowsf_->setOStream(this->getOStream());
  }
  if (!is_null(lows_)) {
    lows_->setVerbLevel(this->getVerbLevel());
    lows_->setOStream(this->getOStream());
  }
}


// Overridden from LinearOpBase

 
template <class Scalar>
bool DelayedLinearOpWithSolve<Scalar>::opSupportedImpl(EOpTransp M_trans) const
{
  return Thyra::opSupported(*fwdOp_,M_trans);
}


template <class Scalar>
void DelayedLinearOpWithSolve<Scalar>::applyImpl(
  const EOpTransp M_trans,
  const MultiVectorBase<Scalar> &X,
  const Ptr<MultiVectorBase<Scalar> > &Y,
  const Scalar alpha,
  const Scalar beta
  ) const
{
  Thyra::apply(*fwdOp_, M_trans, X, Y, alpha, beta);
}


// Overridden from LinearOpWithSolveBase


template <class Scalar>
bool DelayedLinearOpWithSolve<Scalar>::solveSupportsImpl(EOpTransp M_trans) const
{
  updateSolver();
  return Thyra::solveSupports(*lows_, M_trans);
}


template <class Scalar>
bool DelayedLinearOpWithSolve<Scalar>::solveSupportsSolveMeasureTypeImpl(
  EOpTransp M_trans, const SolveMeasureType& solveMeasureType
  ) const
{
  updateSolver();
  return Thyra::solveSupportsSolveMeasureType(*lows_, M_trans, solveMeasureType);
}


template <class Scalar>
SolveStatus<Scalar>
DelayedLinearOpWithSolve<Scalar>::solveImpl(
  const EOpTransp transp,
  const MultiVectorBase<Scalar> &B,
  const Ptr<MultiVectorBase<Scalar> > &X,
  const Ptr<const SolveCriteria<Scalar> > solveCriteria
  ) const
{
  updateSolver();
  return Thyra::solve(*lows_, transp, B, X, solveCriteria);
}


// private


template <class Scalar>
void DelayedLinearOpWithSolve<Scalar>::updateSolver() const
{
  if (is_null(lows_))
    lows_ = lowsf_->createOp();
  if (!lows_is_valid_) {
    if (!is_null(prec_))
      lowsf_->initializePreconditionedOp(
        fwdOpSrc_,prec_,&*lows_,supportSolveUse_);
    else if (!is_null(approxFwdOpSrc_))
      lowsf_->initializeApproxPreconditionedOp(
        fwdOpSrc_,approxFwdOpSrc_,&*lows_,supportSolveUse_);
    else
      lowsf_->initializeOp(
        fwdOpSrc_,&*lows_,supportSolveUse_);
    lows_is_valid_ = true;
  }
}


} // namespace Thyra


#endif // THYRA_DELAYED_LINEAR_OP_WITH_SOLVE_HPP





