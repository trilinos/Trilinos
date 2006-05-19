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

#ifndef THYRA_MODEL_EVALUATOR_HELPERS_HPP
#define THYRA_MODEL_EVALUATOR_HELPERS_HPP

#include "Thyra_ModelEvaluator.hpp"

namespace Thyra {

/** \relates ModelEvaluator */
template<class Scalar>
ModelEvaluatorBase::DerivativeMultiVector<Scalar>
create_DfDp_mv( const ModelEvaluator<Scalar>& model, int l, ModelEvaluatorBase::EDerivativeMultiVectorOrientation orientation )
{
  TEST_FOR_EXCEPT(!(orientation==ModelEvaluatorBase::DERIV_MV_BY_COL));
  return createMembers( model.get_f_space(), model.get_p_space(l)->dim() );
}

/** \relates ModelEvaluator */
template<class Scalar>
ModelEvaluatorBase::DerivativeMultiVector<Scalar>
create_DgDx_mv( const ModelEvaluator<Scalar>& model, int j, ModelEvaluatorBase::EDerivativeMultiVectorOrientation orientation )
{
  typedef ModelEvaluatorBase MEB;
  switch(orientation) {
    case MEB::DERIV_MV_BY_COL:
      return createMembers( model.get_g_space(j), model.get_x_space()->dim() );
    case MEB::DERIV_TRANS_MV_BY_ROW:
      return createMembers( model.get_x_space(j), model.get_g_space()->dim() );
    default:
      TEST_FOR_EXCEPT(true);
  }
  return ModelEvaluatorBase::DerivativeMultiVector<Scalar>(); // Never be executed!
}

/** \relates ModelEvaluator */
template<class Scalar>
ModelEvaluatorBase::DerivativeMultiVector<Scalar>
create_DgDp_mv( const ModelEvaluator<Scalar>& model, int j, int l, ModelEvaluatorBase::EDerivativeMultiVectorOrientation orientation )
{
  typedef ModelEvaluatorBase MEB;
  switch(orientation) {
    case MEB::DERIV_MV_BY_COL:
      return createMembers( model.get_g_space(j), model.get_p_space(l)->dim() );
    case MEB::DERIV_TRANS_MV_BY_ROW:
      return createMembers( model.get_p_space(l), model.get_g_space(j)->dim() );
    default:
      TEST_FOR_EXCEPT(true);
  }
  return ModelEvaluatorBase::DerivativeMultiVector<Scalar>(); // Never be executed!
}

/** \relates ModelEvaluator */
template<class Scalar>
ModelEvaluatorBase::DerivativeMultiVector<Scalar>
get_dmv(
  const ModelEvaluatorBase::Derivative<Scalar>             &deriv
  ,const std::string                                       &derivName
  )
{
  TEST_FOR_EXCEPTION(
    deriv.getLinearOp().get()!=NULL, std::logic_error
    ,"Error, LinearOpBase type not expected for " << derivName <<"!"
    );
  return deriv.getDerivativeMultiVector();
}

/** \relates ModelEvaluator */
template<class Scalar>
Teuchos::RefCountPtr<MultiVectorBase<Scalar> >
get_mv(
  const ModelEvaluatorBase::Derivative<Scalar>             &deriv
  ,const std::string                                       &derivName
  ,ModelEvaluatorBase::EDerivativeMultiVectorOrientation   orientation
  )
{
  typedef ModelEvaluatorBase MEB;
  TEST_FOR_EXCEPTION(
    deriv.getLinearOp().get()!=NULL, std::logic_error
    ,"Error, LinearOpBase type not expected for " << derivName <<"!"
    );
  MEB::DerivativeMultiVector<Scalar>
    dmv = deriv.getDerivativeMultiVector();
  Teuchos::RefCountPtr<MultiVectorBase<Scalar> >
    mv = dmv.getMultiVector();
  if( mv.get() ) {
    TEST_FOR_EXCEPTION(
      dmv.getOrientation() != orientation, std::logic_error
      ,"Error, the orientation " << toString(dmv.getOrientation()) << " is not the"
      " expected orientation of " << toString(orientation)
      << " for " << derivName << "!"
      );
  }
  return mv;
}

/** \brief Evaluate <tt>f(x)</tt>. */
template<class Scalar>
void eval_f(
  const ModelEvaluator<Scalar>                                    &model
  ,const VectorBase<Scalar>                                       &x
  ,VectorBase<Scalar>                                             *f
  )
{
  typedef ModelEvaluatorBase MEB;
  MEB::InArgs<Scalar>   inArgs  = model.createInArgs();
  MEB::OutArgs<Scalar>  outArgs = model.createOutArgs();
  inArgs.set_x(Teuchos::rcp(&x,false));
  outArgs.set_f(Teuchos::rcp(f,false));
  model.evalModel(inArgs,outArgs);
}

/** \brief Evaluate <tt>f(x)</tt> and <tt>W(x) = beta*DfDx(x)</tt>. */
template<class Scalar>
void eval_f_W(
  const ModelEvaluator<Scalar>                                    &model
  ,const VectorBase<Scalar>                                       &x
  ,VectorBase<Scalar>                                             *f
  ,LinearOpWithSolveBase<Scalar>                                  *W
  )
{

  typedef ModelEvaluatorBase MEB;

  MEB::InArgs<Scalar>   inArgs  = model.createInArgs();
  MEB::OutArgs<Scalar>  outArgs = model.createOutArgs();

  inArgs.set_x(Teuchos::rcp(&x,false));

  outArgs.set_f(Teuchos::rcp(f,false));
  if(W) outArgs.set_W(Teuchos::rcp(W,false));

  model.evalModel(inArgs,outArgs);

}

/** \brief Evaluate <tt>f(x)</tt> and <tt>W(x) = beta*DfDx(x)</tt>. */
template<class Scalar>
void eval_f_W(
  const ModelEvaluator<Scalar>                                    &model
  ,const VectorBase<Scalar>                                       &x
  ,const Scalar                                                   &beta
  ,VectorBase<Scalar>                                             *f
  ,LinearOpWithSolveBase<Scalar>                                  *W
  )
{

  typedef ModelEvaluatorBase MEB;

  MEB::InArgs<Scalar>   inArgs  = model.createInArgs();
  MEB::OutArgs<Scalar>  outArgs = model.createOutArgs();

  inArgs.set_x(Teuchos::rcp(&x,false));
  inArgs.set_beta(beta);

  outArgs.set_f(Teuchos::rcp(f,false));
  if(W) outArgs.set_W(Teuchos::rcp(W,false));

  model.evalModel(inArgs,outArgs);

}

/** \brief Evaluate <tt>f(x,t)</tt>. */
template<class Scalar>
void eval_f(
  const ModelEvaluator<Scalar>                                    &model
  ,const VectorBase<Scalar>                                       &x
  ,const Scalar                                                   &t
  ,VectorBase<Scalar>                                             *f
  )
{
  typedef ModelEvaluatorBase MEB;
  MEB::InArgs<Scalar>   inArgs  = model.createInArgs();
  MEB::OutArgs<Scalar>  outArgs = model.createOutArgs();
  inArgs.set_x(Teuchos::rcp(&x,false));
  if(inArgs.supports(MEB::IN_ARG_t)) inArgs.set_t(t);
  outArgs.set_f(Teuchos::rcp(f,false));
  model.evalModel(inArgs,outArgs);
}


/** \brief Evaluate <tt>f(x_dot,x,t)</tt>. */
template<class Scalar>
void eval_f(
  const ModelEvaluator<Scalar>                                    &model
  ,const VectorBase<Scalar>                                       &x_dot
  ,const VectorBase<Scalar>                                       &x
  ,const typename ModelEvaluatorBase::InArgs<Scalar>::ScalarMag   &t
  ,VectorBase<Scalar>                                             *f
  )
{

  typedef ModelEvaluatorBase MEB;

  MEB::InArgs<Scalar>   inArgs  = model.createInArgs();
  MEB::OutArgs<Scalar>  outArgs = model.createOutArgs();

  inArgs.set_x_dot(Teuchos::rcp(&x_dot,false));
  inArgs.set_x(Teuchos::rcp(&x,false));
  if(inArgs.supports(MEB::IN_ARG_t))
    inArgs.set_t(t);

  outArgs.set_f(Teuchos::rcp(f,false));

  model.evalModel(inArgs,outArgs);

}

/** \brief Evaluate <tt>f(x_dot,x,t)</tt> and
 * <tt>W(x_dot,x,t,alpha,beta) = beta*DfDx_dot(x_dot,x,t)beta*DfDx(x_dot,x,t)</tt>. */
template<class Scalar>
void eval_f_W(
  const ModelEvaluator<Scalar>                                    &model
  ,const VectorBase<Scalar>                                       &x_dot
  ,const VectorBase<Scalar>                                       &x
  ,const typename ModelEvaluatorBase::InArgs<Scalar>::ScalarMag   &t
  ,const Scalar                                                   &alpha
  ,const Scalar                                                   &beta
  ,VectorBase<Scalar>                                             *f
  ,LinearOpWithSolveBase<Scalar>                                  *W
  )
{

  typedef ModelEvaluatorBase MEB;

  MEB::InArgs<Scalar>   inArgs  = model.createInArgs();
  MEB::OutArgs<Scalar>  outArgs = model.createOutArgs();

  inArgs.set_x_dot(Teuchos::rcp(&x_dot,false));
  inArgs.set_x(Teuchos::rcp(&x,false));
  if(inArgs.supports(MEB::IN_ARG_t))
    inArgs.set_t(t);
  inArgs.set_alpha(alpha);
  inArgs.set_beta(beta);

  if(f) outArgs.set_f(Teuchos::rcp(f,false));
  if(W) outArgs.set_W(Teuchos::rcp(W,false));

  model.evalModel(inArgs,outArgs);

}

/** \brief . */
template<class Scalar>
void eval_f_poly(
  const ModelEvaluator<Scalar>                                    &model
  ,const Teuchos::Polynomial< VectorBase<Scalar> >                &x_poly
  ,const typename ModelEvaluatorBase::InArgs<Scalar>::ScalarMag   &t
  ,Teuchos::Polynomial< VectorBase<Scalar> >                      *f_poly
  )
{

  typedef ModelEvaluatorBase MEB;

  MEB::InArgs<Scalar>   inArgs  = model.createInArgs();
  MEB::OutArgs<Scalar>  outArgs = model.createOutArgs();

  inArgs.set_x_poly(Teuchos::rcp(&x_poly,false));
  if(inArgs.supports(MEB::IN_ARG_t))
    inArgs.set_t(t);

  outArgs.set_f_poly(Teuchos::rcp(f_poly,false));

  model.evalModel(inArgs,outArgs);

}

/** \brief . */
template<class Scalar>
void eval_f_poly(
  const ModelEvaluator<Scalar>                                    &model
  ,const Teuchos::Polynomial< VectorBase<Scalar> >                &x_dot_poly
  ,const VectorBase<Scalar>                                       &x_poly
  ,const typename ModelEvaluatorBase::InArgs<Scalar>::ScalarMag   &t
  ,Teuchos::Polynomial< VectorBase<Scalar> >                      *f_poly
  )
{

  typedef ModelEvaluatorBase MEB;

  MEB::InArgs<Scalar>   inArgs  = model.createInArgs();
  MEB::OutArgs<Scalar>  outArgs = model.createOutArgs();

  inArgs.set_x_dot_poly(Teuchos::rcp(&x_dot_poly,false));
  inArgs.set_x_poly(Teuchos::rcp(&x_poly,false));
  if(inArgs.supports(MEB::IN_ARG_t))
    inArgs.set_t(t);

  outArgs.set_f_poly(Teuchos::rcp(f_poly,false));

  model.evalModel(inArgs,outArgs);

}

} // namespace Thyra


#endif // THYRA_MODEL_EVALUATOR_HELPERS_HPP
