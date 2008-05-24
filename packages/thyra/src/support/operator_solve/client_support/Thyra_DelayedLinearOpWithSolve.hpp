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


#ifndef THYRA_DELAYED_LINEAR_OP_WITH_SOLVE_HPP
#define THYRA_DELAYED_LINEAR_OP_WITH_SOLVE_HPP


#include "Thyra_DelayedLinearOpWithSolveDecl.hpp"
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
    oss << "NULL";
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


// Overridden from SingleScalarLinearOpBase

 
template <class Scalar>
bool DelayedLinearOpWithSolve<Scalar>::opSupported(EOpTransp M_trans) const
{
  return Thyra::opSupported(*fwdOp_,M_trans);
}


template <class Scalar>
void DelayedLinearOpWithSolve<Scalar>::apply(
  const EOpTransp M_trans,
  const MultiVectorBase<Scalar> &X,
    MultiVectorBase<Scalar> *Y,
  const Scalar alpha,
  const Scalar beta
  ) const
{
  Thyra::apply(*fwdOp_,M_trans,X,Y,alpha,beta);
}


template <class Scalar>
bool DelayedLinearOpWithSolve<Scalar>::solveSupportsTrans(
  EOpTransp M_trans
  ) const
{
  updateSolver();
  return Thyra::solveSupports(*lows_,M_trans);
}


template <class Scalar>
bool DelayedLinearOpWithSolve<Scalar>::solveSupportsSolveMeasureType(
  EOpTransp M_trans, const SolveMeasureType& solveMeasureType
  ) const
{
  updateSolver();
  return Thyra::solveSupportsSolveMeasureType(*lows_,M_trans,solveMeasureType);
}


template <class Scalar>
void DelayedLinearOpWithSolve<Scalar>::solve(
  const EOpTransp M_trans,
  const MultiVectorBase<Scalar> &B,
  MultiVectorBase<Scalar> *X,
  const int numBlocks,
  const BlockSolveCriteria<Scalar> blockSolveCriteria[],
  SolveStatus<Scalar> blockSolveStatus[]
  ) const
{
  updateSolver();
  Thyra::solve(
    *lows_, M_trans, B, X, numBlocks, blockSolveCriteria, blockSolveStatus );
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





