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
