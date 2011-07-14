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

#ifndef THYRA_LINEAR_OP_WITH_SOLVE_FACTORY_BASE_DEF_HPP
#define THYRA_LINEAR_OP_WITH_SOLVE_FACTORY_BASE_DEF_HPP


#include "Thyra_LinearOpWithSolveFactoryBase_decl.hpp"


namespace Thyra {


template<class Scalar>
bool LinearOpWithSolveFactoryBase<Scalar>::acceptsPreconditionerFactory() const
{
  return false;
}


template<class Scalar>
void LinearOpWithSolveFactoryBase<Scalar>::setPreconditionerFactory(
  const RCP<PreconditionerFactoryBase<Scalar> > &precFactory
  ,const std::string &precFactoryName
  )
{
  TEST_FOR_EXCEPTION(
    true,std::logic_error
    ,"Error, the concrete implementation described as \'"<<this->description()<<"\' did not override this "
    "setPreconditionerFactory(...) function and the default implementation throws this exception!"
    );
}


template<class Scalar>
RCP<PreconditionerFactoryBase<Scalar> >
LinearOpWithSolveFactoryBase<Scalar>::getPreconditionerFactory() const
{
  return Teuchos::null;
}


template<class Scalar>
void LinearOpWithSolveFactoryBase<Scalar>::unsetPreconditionerFactory(
  RCP<PreconditionerFactoryBase<Scalar> > *precFactory
  ,std::string *precFactoryName
  )
{
  if(precFactory) *precFactory = Teuchos::null;
  if(precFactoryName) *precFactoryName = "";
}


template<class Scalar>
void LinearOpWithSolveFactoryBase<Scalar>::initializeAndReuseOp(
  const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc
  ,LinearOpWithSolveBase<Scalar> *Op
  ) const
{
  this->initializeOp(fwdOpSrc,Op);
}


template<class Scalar>
bool LinearOpWithSolveFactoryBase<Scalar>::supportsPreconditionerInputType(
  const EPreconditionerInputType precOpType
  ) const
{
  return false;
}


template<class Scalar>
void LinearOpWithSolveFactoryBase<Scalar>::initializePreconditionedOp(
  const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc
  ,const RCP<const PreconditionerBase<Scalar> > &prec
  ,LinearOpWithSolveBase<Scalar> *Op
  ,const ESupportSolveUse supportSolveUse
  ) const
{
  TEST_FOR_EXCEPTION(
    true,std::logic_error
    ,"Error, the concrete implementation described as \'"<<this->description()<<"\' did not override this "
    "initializePreconditionedOp(...) function and the default implementation throws this exception!"
    );
}


template<class Scalar>
void LinearOpWithSolveFactoryBase<Scalar>::initializeApproxPreconditionedOp(
  const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc
  ,const RCP<const LinearOpSourceBase<Scalar> > &approxFwdOpSrc
  ,LinearOpWithSolveBase<Scalar> *Op
  ,const ESupportSolveUse supportSolveUse
  ) const
{
  TEST_FOR_EXCEPTION(
    true,std::logic_error
    ,"Error, the concrete implementation described as \'"<<this->description()<<"\' did not override this "
    "initializePreconditionedOp(...) function and the default implementation throws this exception!"
    );
}


} // namespace Thyra


#endif // THYRA_LINEAR_OP_WITH_SOLVE_FACTORY_BASE_DEF_HPP
