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

#ifndef THYRA_LINEAR_OP_WITH_SOLVE_FACTORY_HELPERS_HPP
#define THYRA_LINEAR_OP_WITH_SOLVE_FACTORY_HELPERS_HPP

#include "Thyra_LinearOpWithSolveFactoryBase.hpp"
#include "Thyra_DefaultLinearOpSource.hpp"
#include "Thyra_DefaultPreconditioner.hpp"

namespace Thyra {

template<class Scalar>
bool isCompatible(
  const LinearOpWithSolveFactoryBase<Scalar>    &lowsFactory
  ,const LinearOpBase<Scalar>                   &fwdOp
  )
{
  return lowsFactory.isCompatible(*defaultLinearOpSource(fwdOp));
}

/** \brief Create and initialize a <tt>LinearOpWithSolveBase</tt> object from
 * a <tt>LinearOpBase</tt> object using a
 * <tt>LinearOpWithSolveFactoryBase</tt> strategy object.
 */
template<class Scalar>
Teuchos::RefCountPtr<LinearOpWithSolveBase<Scalar> >
linearOpWithSolve(
  const LinearOpWithSolveFactoryBase<Scalar>                  &lowsFactory
  ,const Teuchos::RefCountPtr<const LinearOpBase<Scalar> >    &fwdOp
  ,const ESupportSolveUse                                     supportSolveUse = SUPPORT_SOLVE_UNSPECIFIED
  )
{
  Teuchos::RefCountPtr<LinearOpWithSolveBase<Scalar> >
    Op = lowsFactory.createOp();
  lowsFactory.initializeOp(defaultLinearOpSource(fwdOp),&*Op,supportSolveUse);
  return Op;
}

template<class Scalar>
void initializeOp(
  const LinearOpWithSolveFactoryBase<Scalar>                  &lowsFactory
  ,const Teuchos::RefCountPtr<const LinearOpBase<Scalar> >    &fwdOp
  ,LinearOpWithSolveBase<Scalar>                              *Op
  ,const ESupportSolveUse                                     supportSolveUse = SUPPORT_SOLVE_UNSPECIFIED
  )
{
  lowsFactory.initializeOp(defaultLinearOpSource(fwdOp),Op,supportSolveUse);
}

template<class Scalar>
void initializeAndReuseOp(
  const LinearOpWithSolveFactoryBase<Scalar>                  &lowsFactory
  ,const Teuchos::RefCountPtr<const LinearOpBase<Scalar> >    &fwdOp
  ,LinearOpWithSolveBase<Scalar>                              *Op
  )
{
  lowsFactory.initializeAndReuseOp(defaultLinearOpSource(fwdOp),Op);
}
  
template<class Scalar>
void uninitializeOp(
  const LinearOpWithSolveFactoryBase<Scalar>                  &lowsFactory
  ,LinearOpWithSolveBase<Scalar>                              *Op
  ,Teuchos::RefCountPtr<const LinearOpBase<Scalar> >          *fwdOp           = NULL
  ,Teuchos::RefCountPtr<const PreconditionerBase<Scalar> >    *prec            = NULL
  ,Teuchos::RefCountPtr<const LinearOpBase<Scalar> >          *approxFwdOp     = NULL
  ,ESupportSolveUse                                           *supportSolveUse = NULL
  )
{
  Teuchos::RefCountPtr<const LinearOpSourceBase<Scalar> > fwdOpSrc;
  Teuchos::RefCountPtr<const LinearOpSourceBase<Scalar> > approxFwdOpSrc;
  lowsFactory.uninitializeOp(Op,&fwdOpSrc,prec,&approxFwdOpSrc,supportSolveUse);
  if(fwdOp)
    *fwdOp = ( fwdOpSrc.get() ? fwdOpSrc->getOp() : Teuchos::null );
  if(approxFwdOp)
    *approxFwdOp = ( approxFwdOpSrc.get() ? approxFwdOpSrc->getOp() : Teuchos::null );
}

template<class Scalar>
void initializePreconditionedOp(
  const LinearOpWithSolveFactoryBase<Scalar>                       &lowsFactory
  ,const Teuchos::RefCountPtr<const LinearOpBase<Scalar> >         &fwdOp
  ,const Teuchos::RefCountPtr<const PreconditionerBase<Scalar> >   &prec
  ,LinearOpWithSolveBase<Scalar>                                   *Op
  ,const ESupportSolveUse                                          supportSolveUse = SUPPORT_SOLVE_UNSPECIFIED
  )
{
  lowsFactory.initializePreconditionedOp(defaultLinearOpSource(fwdOp),prec,Op,supportSolveUse);
}

template<class Scalar>
void initializeApproxPreconditionedOp(
  const LinearOpWithSolveFactoryBase<Scalar>                    &lowsFactory
  ,const Teuchos::RefCountPtr<const LinearOpBase<Scalar> >      &fwdOp
  ,const Teuchos::RefCountPtr<const LinearOpBase<Scalar> >      &approxFwdOp
  ,LinearOpWithSolveBase<Scalar>                                *Op
  ,const ESupportSolveUse                                       supportSolveUse = SUPPORT_SOLVE_UNSPECIFIED
  )
{
  lowsFactory.initializeApproxPreconditionedOp(
    defaultLinearOpSource(fwdOp),defaultLinearOpSource(approxFwdOp),Op,supportSolveUse
    );
}

} // namespace Thyra

#endif // THYRA_LINEAR_OP_WITH_SOLVE_FACTORY_HELPERS_HPP
