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

#ifndef THYRA_DELAYED_LINEAR_OP_WITH_SOLVE_FACTORY_HPP
#define THYRA_DELAYED_LINEAR_OP_WITH_SOLVE_FACTORY_HPP


#include "Thyra_DelayedLinearOpWithSolveFactory_decl.hpp"
#include "Thyra_LinearOpWithSolveBase.hpp"
#include "Thyra_DelayedLinearOpWithSolve.hpp"


namespace Thyra {


// Overridden from Constructors/Initializers/Accessors

  
template<class Scalar>
DelayedLinearOpWithSolveFactory<Scalar>::DelayedLinearOpWithSolveFactory(
  const RCP<LinearOpWithSolveFactoryBase<Scalar> > &lowsf
  )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(is_null(lowsf));
#endif
  lowsf_ = lowsf;
}

template<class Scalar>
RCP<LinearOpWithSolveFactoryBase<Scalar> >
DelayedLinearOpWithSolveFactory<Scalar>::getUnderlyingLOWSF()
{
  return lowsf_;
}


template<class Scalar>
RCP<const LinearOpWithSolveFactoryBase<Scalar> >
DelayedLinearOpWithSolveFactory<Scalar>::getUnderlyingLOWSF() const
{
  return lowsf_;
}


// Overridden from Teuchos::Describable


template<class Scalar>
std::string DelayedLinearOpWithSolveFactory<Scalar>::description() const
{
  std::ostringstream oss;
  oss << this->Teuchos::Describable::description()
      << "{"
      << "lowsf=";
  if (!is_null(lowsf_))
    oss << lowsf_->description();
  else
    oss << "NULL";
  oss << "}";
  return oss.str();
}


// Overridden from ParameterListAcceptor


template<class Scalar>
void DelayedLinearOpWithSolveFactory<Scalar>::setParameterList(
  RCP<ParameterList> const& paramList
  )
{
  lowsf_->setParameterList(paramList);
}


template<class Scalar>
RCP<ParameterList>
DelayedLinearOpWithSolveFactory<Scalar>::getNonconstParameterList()
{
  return lowsf_->getNonconstParameterList();
}


template<class Scalar>
RCP<ParameterList> 
DelayedLinearOpWithSolveFactory<Scalar>::unsetParameterList()
{
  return lowsf_->unsetParameterList();
}


template<class Scalar>
RCP<const ParameterList>
DelayedLinearOpWithSolveFactory<Scalar>::getParameterList() const
{
  return lowsf_->getParameterList();
}


template<class Scalar>
RCP<const ParameterList>
DelayedLinearOpWithSolveFactory<Scalar>::getValidParameters() const
{
  return lowsf_->getValidParameters();
}


// Overridden from LinearOpWithSolveFactoyBase


template<class Scalar>
bool DelayedLinearOpWithSolveFactory<Scalar>::acceptsPreconditionerFactory() const
{
  return lowsf_->acceptsPreconditionerFactory();
}


template<class Scalar>
void DelayedLinearOpWithSolveFactory<Scalar>::setPreconditionerFactory(
  const RCP<PreconditionerFactoryBase<Scalar> > &precFactory,
  const std::string &precFactoryName
  )
{
  lowsf_->setPreconditionerFactory(precFactory,precFactoryName);
}


template<class Scalar>
RCP<PreconditionerFactoryBase<Scalar> >
DelayedLinearOpWithSolveFactory<Scalar>::getPreconditionerFactory() const
{
  return lowsf_->getPreconditionerFactory();
}


template<class Scalar>
void DelayedLinearOpWithSolveFactory<Scalar>::unsetPreconditionerFactory(
  RCP<PreconditionerFactoryBase<Scalar> > *precFactory,
  std::string *precFactoryName
  )
{
  lowsf_->unsetPreconditionerFactory(precFactory);
}


template<class Scalar>
bool DelayedLinearOpWithSolveFactory<Scalar>::isCompatible(
  const LinearOpSourceBase<Scalar> &fwdOpSrc
  ) const
{
  return lowsf_->isCompatible(fwdOpSrc);
}


template<class Scalar>
RCP<LinearOpWithSolveBase<Scalar> >
DelayedLinearOpWithSolveFactory<Scalar>::createOp() const
{
  RCP<LinearOpWithSolveBase<Scalar> >
    dlows = Teuchos::rcp(new DelayedLinearOpWithSolve<Scalar>());
  dlows->setVerbLevel(this->getVerbLevel());
  dlows->setOStream(this->getOStream());
  return dlows;
}


template<class Scalar>
void DelayedLinearOpWithSolveFactory<Scalar>::initializeOp(
  const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
  LinearOpWithSolveBase<Scalar> *Op,
  const ESupportSolveUse supportSolveUse
  ) const
{
  using Teuchos::null;
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(is_null(fwdOpSrc));
  TEST_FOR_EXCEPT(0==Op);
#endif
  DelayedLinearOpWithSolve<Scalar>
    &dlows = Teuchos::dyn_cast<DelayedLinearOpWithSolve<Scalar> >(*Op);
  dlows.initialize( fwdOpSrc, null, null, supportSolveUse, lowsf_ );
}


template<class Scalar>
void DelayedLinearOpWithSolveFactory<Scalar>::initializeAndReuseOp(
  const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
  LinearOpWithSolveBase<Scalar> *Op
  ) const
{
  TEST_FOR_EXCEPT(true);
}


template<class Scalar>
void DelayedLinearOpWithSolveFactory<Scalar>::uninitializeOp(
  LinearOpWithSolveBase<Scalar> *Op,
  RCP<const LinearOpSourceBase<Scalar> > *fwdOpSrc,
  RCP<const PreconditionerBase<Scalar> > *prec,
  RCP<const LinearOpSourceBase<Scalar> > *approxFwdOpSrc,
  ESupportSolveUse *supportSolveUse
  ) const
{

  using Teuchos::dyn_cast;

#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(0==Op);
#endif

  DelayedLinearOpWithSolve<Scalar>
    &dlows = dyn_cast<DelayedLinearOpWithSolve<Scalar> >(*Op);
  
  if (fwdOpSrc)
    *fwdOpSrc = dlows.getFwdOpSrc();
  if (prec)
    *prec = dlows.getPrec();
  if (approxFwdOpSrc)
    *approxFwdOpSrc = dlows.getApproxFwdOpSrc();
  if (supportSolveUse)
    *supportSolveUse = dlows.getSupportSolveUse();

  // ToDo: 2007/08/16: rabartl: Consider uninitalizing dlows?

}


template<class Scalar>
bool DelayedLinearOpWithSolveFactory<Scalar>::supportsPreconditionerInputType(
  const EPreconditionerInputType precOpType
  ) const
{
  return lowsf_->supportsPreconditionerInputType(precOpType);
}


template<class Scalar>
void DelayedLinearOpWithSolveFactory<Scalar>::initializePreconditionedOp(
  const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
  const RCP<const PreconditionerBase<Scalar> > &prec,
  LinearOpWithSolveBase<Scalar> *Op,
  const ESupportSolveUse supportSolveUse
  ) const
{
  using Teuchos::null;
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(is_null(fwdOpSrc));
  TEST_FOR_EXCEPT(0==Op);
#endif
  DelayedLinearOpWithSolve<Scalar>
    &dlows = Teuchos::dyn_cast<DelayedLinearOpWithSolve<Scalar> >(*Op);
  dlows.initialize( fwdOpSrc, prec, null, supportSolveUse, lowsf_ );
}


template<class Scalar>
void DelayedLinearOpWithSolveFactory<Scalar>::initializeApproxPreconditionedOp(
  const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
  const RCP<const LinearOpSourceBase<Scalar> > &approxFwdOpSrc,
  LinearOpWithSolveBase<Scalar> *Op,
  const ESupportSolveUse supportSolveUse
  ) const
{
  using Teuchos::null;
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(is_null(fwdOpSrc));
  TEST_FOR_EXCEPT(0==Op);
#endif
  DelayedLinearOpWithSolve<Scalar>
    &dlows = Teuchos::dyn_cast<DelayedLinearOpWithSolve<Scalar> >(*Op);
  dlows.initialize( fwdOpSrc, null, approxFwdOpSrc, supportSolveUse, lowsf_ );
}


// protected


template<class Scalar>
void DelayedLinearOpWithSolveFactory<Scalar>::informUpdatedVerbosityState() const
{
  lowsf_->setVerbLevel(this->getVerbLevel());
  lowsf_->setOStream(this->getOStream());
}


} // namespace Thyra


#endif // THYRA_DELAYED_LINEAR_OP_WITH_SOLVE_FACTORY_HPP
