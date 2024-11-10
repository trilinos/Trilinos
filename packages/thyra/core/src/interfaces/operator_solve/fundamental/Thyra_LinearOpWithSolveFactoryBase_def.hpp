// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
  const RCP<PreconditionerFactoryBase<Scalar> > &/* precFactory */
  ,const std::string &/* precFactoryName */
  )
{
  TEUCHOS_TEST_FOR_EXCEPTION(
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
  const EPreconditionerInputType /* precOpType */
  ) const
{
  return false;
}


template<class Scalar>
void LinearOpWithSolveFactoryBase<Scalar>::initializePreconditionedOp(
  const RCP<const LinearOpSourceBase<Scalar> > &/* fwdOpSrc */
  ,const RCP<const PreconditionerBase<Scalar> > &/* prec */
  ,LinearOpWithSolveBase<Scalar> * /* Op */
  ,const ESupportSolveUse /* supportSolveUse */
  ) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    true,std::logic_error
    ,"Error, the concrete implementation described as \'"<<this->description()<<"\' did not override this "
    "initializePreconditionedOp(...) function and the default implementation throws this exception!"
    );
}


template<class Scalar>
void LinearOpWithSolveFactoryBase<Scalar>::initializeApproxPreconditionedOp(
  const RCP<const LinearOpSourceBase<Scalar> > &/* fwdOpSrc */
  ,const RCP<const LinearOpSourceBase<Scalar> > &/* approxFwdOpSrc */
  ,LinearOpWithSolveBase<Scalar> * /* Op */
  ,const ESupportSolveUse /* supportSolveUse */
  ) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    true,std::logic_error
    ,"Error, the concrete implementation described as \'"<<this->description()<<"\' did not override this "
    "initializePreconditionedOp(...) function and the default implementation throws this exception!"
    );
}


} // namespace Thyra


#endif // THYRA_LINEAR_OP_WITH_SOLVE_FACTORY_BASE_DEF_HPP
