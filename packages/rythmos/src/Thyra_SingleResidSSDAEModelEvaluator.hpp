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

#ifndef THYRA_SINGLE_RESID_SS_DAE_MODEL_EVALUATOR_HPP
#define THYRA_SINGLE_RESID_SS_DAE_MODEL_EVALUATOR_HPP

#include "Thyra_ModelEvaluator.hpp"
#include "Thyra_VectorStdOps.hpp"

#ifdef THYRA_RYTHMOS_DEBUG
#include "Thyra_TestingTools.hpp"
#endif // THYRA_RYTHMOS_DEBUG

namespace Thyra {

/** \brief Decorator subclass for a steady-state version of a DAE for single-residual
 * time stepper methods.
 *
 * This class provides a set of nonlinear equations of the form:
 *
 \verbatim

 f_bar(x_bar) = f( coeff_x_dot * x_bar + x_dot_base, coeff_x * x_bar + x_base, t_base ) = 0

 \endverbatim
 *
 * with initial starting guess <tt>x_bar_init</tt>.
 *
 * The Jacobian of this model is:

 \verbatim

 beta * d(f_bar)/d(x_bar) = (beta*coeff_x_dot) * d(f)/d(x_dot) + (beta*coeff_x) * d(f)/d(x)

 \endverbatim
 *
 * ToDo: Finish Documentation!
 */
template<class Scalar>
class SingleResidSSDAEModelEvaluator : public ModelEvaluator<Scalar> {
public:

  /** \name Constructors/initializers/accessors */
  //@{

  /** \brief . */
  SingleResidSSDAEModelEvaluator();

  /** \brief . */
  SingleResidSSDAEModelEvaluator(
    const Teuchos::RefCountPtr<const ModelEvaluator<Scalar> >     &daeModel
    ,const Scalar                                                 &coeff_x_dot
    ,const Teuchos::RefCountPtr<const VectorBase<Scalar> >        &x_dot_base
    ,const Scalar                                                 &coeff_x
    ,const Teuchos::RefCountPtr<const VectorBase<Scalar> >        &x_base
    ,const Scalar                                                 &t_base
    ,const Teuchos::RefCountPtr<const VectorBase<Scalar> >        &x_bar_init
    );

  /** \brief . */
  void initialize(
    const Teuchos::RefCountPtr<const ModelEvaluator<Scalar> >     &daeModel
    ,const Scalar                                                 &coeff_x_dot
    ,const Teuchos::RefCountPtr<const VectorBase<Scalar> >        &x_dot_base
    ,const Scalar                                                 &coeff_x
    ,const Teuchos::RefCountPtr<const VectorBase<Scalar> >        &x_base
    ,const Scalar                                                 &t_base
    ,const Teuchos::RefCountPtr<const VectorBase<Scalar> >        &x_bar_init
    );

  //@}

  /** \name Overridden from ModelEvaluator */
  //@{

  /** \breif . */
  Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > get_x_space() const;

  /** \breif . */
  Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > get_f_space() const;

  /** \breif . */
  Teuchos::RefCountPtr<const VectorBase<Scalar> > get_x_init() const;

  /** \breif . */
  Teuchos::RefCountPtr<LinearOpWithSolveBase<Scalar> > create_W() const;

  /** \brief . */
  ModelEvaluatorBase::InArgs<Scalar> createInArgs() const;

  /** \brief . */
  ModelEvaluatorBase::OutArgs<Scalar> createOutArgs() const;

  /** \brief . */
  void evalModel( const ModelEvaluatorBase::InArgs<Scalar>& inArgs, const ModelEvaluatorBase::OutArgs<Scalar>& outArgs ) const;

  //@}

private:

  Teuchos::RefCountPtr<const ModelEvaluator<Scalar> >    daeModel_;
  Scalar                                                 coeff_x_dot_;
  Teuchos::RefCountPtr<const VectorBase<Scalar> >        x_dot_base_;
  Scalar                                                 coeff_x_;
  Teuchos::RefCountPtr<const VectorBase<Scalar> >        x_base_;
  Scalar                                                 t_base_;
  Teuchos::RefCountPtr<const VectorBase<Scalar> >        x_bar_init_;

  // cache
  Teuchos::RefCountPtr<VectorBase<Scalar> >        x_;
  Teuchos::RefCountPtr<VectorBase<Scalar> >        x_dot_;


};

// ///////////////////////
// Definition

// Constructors/initializers/accessors

template<class Scalar>
SingleResidSSDAEModelEvaluator<Scalar>::SingleResidSSDAEModelEvaluator()
{
  // Compiler makes me write this!!!! (gcc 3.4.3)
}

template<class Scalar>
SingleResidSSDAEModelEvaluator<Scalar>::SingleResidSSDAEModelEvaluator(
  const Teuchos::RefCountPtr<const ModelEvaluator<Scalar> >     &daeModel
  ,const Scalar                                                 &coeff_x_dot
  ,const Teuchos::RefCountPtr<const VectorBase<Scalar> >        &x_dot_base
  ,const Scalar                                                 &coeff_x
  ,const Teuchos::RefCountPtr<const VectorBase<Scalar> >        &x_base
  ,const Scalar                                                 &t_base
  ,const Teuchos::RefCountPtr<const VectorBase<Scalar> >        &x_bar_init
  )
{
  initialize(daeModel,coeff_x_dot,x_dot_base,coeff_x,x_base,t_base,x_bar_init);
}

template<class Scalar>
void SingleResidSSDAEModelEvaluator<Scalar>::initialize(
  const Teuchos::RefCountPtr<const ModelEvaluator<Scalar> >     &daeModel
  ,const Scalar                                                 &coeff_x_dot
  ,const Teuchos::RefCountPtr<const VectorBase<Scalar> >        &x_dot_base
  ,const Scalar                                                 &coeff_x
  ,const Teuchos::RefCountPtr<const VectorBase<Scalar> >        &x_base
  ,const Scalar                                                 &t_base
  ,const Teuchos::RefCountPtr<const VectorBase<Scalar> >        &x_bar_init
  )
{
  daeModel_      = daeModel;
  coeff_x_dot_   = coeff_x_dot;
  x_dot_base_    = x_dot_base;
  coeff_x_       = coeff_x;
  x_base_        = x_base;
  t_base_        = t_base;
  x_bar_init_    = x_bar_init;

  x_dot_ = createMember( daeModel_->get_x_space() );
  x_ = createMember( daeModel_->get_x_space() );

  // ToDo: Check that daeModel supports x_dot, x and maybe t

#ifdef THYRA_RYTHMOS_DEBUG
  std::cout << "----------------------------------------------------------------------" << std::endl;
  std::cout << "Thyra::SingleResidSSDAEModelEvaluator::initialize" << std::endl;
  std::cout << "coeff_x_dot_ = " << coeff_x_dot_ << std::endl;
  std::cout << "x_dot_base_ = ";                  
  if ( x_dot_base_.get() ) 
    std::cout << "\n" << *x_dot_base_            << std::endl;
  else
    std::cout << "null"                          << std::endl;
  std::cout << "coeff_x_ = " << coeff_x_         << std::endl;
  std::cout << "x_base_ = ";                      
  if ( x_base_.get() )
    std::cout << "\n" << *x_base_                << std::endl;
  else
    std::cout << "null"                          << std::endl;
  std::cout << "t_base_ = " << t_base_           << std::endl;
  std::cout << "x_bar_init_ = ";                  
  if ( x_bar_init_.get() )
    std::cout << "\n" <<  *x_bar_init_           << std::endl;
  else
    std::cout << "null"                          << std::endl;
  std::cout << "x_dot_ = ";                       
  if ( x_dot_.get() )
    std::cout << "\n" << *x_dot_                 << std::endl;
  else
    std::cout << "null"                          << std::endl;
  std::cout << "x_ = ";                           
  if ( x_.get() )
    std::cout << "\n" << *x_                     << std::endl;
  else
    std::cout << "null"                          << std::endl;
  std::cout << "----------------------------------------------------------------------" << std::endl;
#endif // THYRA_RYTHMOS_DEBUG
}

// Overridden from ModelEvaluator

template<class Scalar>
Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> >
SingleResidSSDAEModelEvaluator<Scalar>::get_x_space() const
{
  return daeModel_->get_x_space();
}

template<class Scalar>
Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> >
SingleResidSSDAEModelEvaluator<Scalar>::get_f_space() const
{
  return daeModel_->get_f_space();
}

template<class Scalar>
Teuchos::RefCountPtr<const VectorBase<Scalar> >
SingleResidSSDAEModelEvaluator<Scalar>::get_x_init() const
{
  return x_bar_init_;
}

template<class Scalar>
Teuchos::RefCountPtr<LinearOpWithSolveBase<Scalar> >
SingleResidSSDAEModelEvaluator<Scalar>::create_W() const
{
  return daeModel_->create_W();
}

template<class Scalar>
ModelEvaluatorBase::InArgs<Scalar>
SingleResidSSDAEModelEvaluator<Scalar>::createInArgs() const
{
  typedef ModelEvaluatorBase MEB;
  MEB::InArgsSetup<Scalar> inArgs;
  inArgs.setSupports(MEB::IN_ARG_x);
  inArgs.setSupports(MEB::IN_ARG_beta);
  return inArgs;
}

template<class Scalar>
ModelEvaluatorBase::OutArgs<Scalar>
SingleResidSSDAEModelEvaluator<Scalar>::createOutArgs() const
{
  typedef ModelEvaluatorBase MEB;
  MEB::OutArgsSetup<Scalar> outArgs;
  outArgs.setSupports(MEB::OUT_ARG_f);
  return outArgs;
}

template<class Scalar>
void SingleResidSSDAEModelEvaluator<Scalar>::evalModel(
  const ModelEvaluatorBase::InArgs<Scalar>&   inArgs_bar
  ,const ModelEvaluatorBase::OutArgs<Scalar>& outArgs_bar
  ) const
{
#ifdef THYRA_RYTHMOS_DEBUG
  std::cout << "----------------------------------------------------------------------" << std::endl;
  std::cout << "Thyra::SingleResidSSDAEModelEvaluator::evalModel" << std::endl;
#endif // THYRA_RYTHMOS_DEBUG
  const VectorBase<Scalar> &x_bar = *inArgs_bar.get_x(); 
  // x_dot = coeff_x_dot * x_bar + x_dot_base
  if (x_dot_base_.get())
    Thyra::V_StVpV( &*x_dot_, coeff_x_dot_, x_bar, *x_dot_base_ );
  else
    Thyra::V_StV( &*x_dot_, coeff_x_dot_, x_bar);
#ifdef THYRA_RYTHMOS_DEBUG
  std::cout << "x_dot_ = coeff_x_dot_ * x_bar + x_dot_base_" << std::endl;
  std::cout << "coeff_x_dot_ = " << coeff_x_dot_             << std::endl;
  std::cout << "x_bar = "                                    << std::endl;
  std::cout <<  x_bar                                        << std::endl;
  std::cout << "x_dot_base_ = ";
  if ( x_dot_base_.get() )
    std::cout << "\n" << *x_dot_base_                        << std::endl;
  else
    std::cout << "null"                                      << std::endl;
  std::cout << "x_dot_ = ";
  if ( x_dot_.get() )
    std::cout << "\n" << *x_dot_                             << std::endl;
  else
    std::cout << "null"                                      << std::endl;
#endif // THYRA_RYTHMOS_DEBUG

  // x = coeff_x * x_bar + x_base
  if (x_base_.get())
    Thyra::V_StVpV( &*x_, coeff_x_, x_bar, *x_base_ );
  else
    Thyra::V_StV( &*x_, coeff_x_, x_bar);
#ifdef THYRA_RYTHMOS_DEBUG
  std::cout << "x_ = coeff_x_ * x_bar + x_base_" << std::endl;
  std::cout << "coeff_x_ = " << coeff_x_         << std::endl;
  std::cout << "x_bar = "                        << std::endl;
  std::cout <<  x_bar                            << std::endl;
  std::cout << "x_base_ = ";
  if ( x_base_.get() )
    std::cout << "\n" << *x_base_                << std::endl;
  else
    std::cout << "null"                          << std::endl;
  std::cout << "x_ = ";
  if ( x_.get() )
    std::cout << "\n" << *x_                     << std::endl;
  else
    std::cout << "null"                          << std::endl;
#endif // THYRA_RYTHMOS_DEBUG

  // Compute W and f
  Teuchos::RefCountPtr<LinearOpWithSolveBase<Scalar> > W;
  if( (W = outArgs_bar.get_W()).get() ) {
    // Compute Jacobian and the residual
    const Scalar beta = inArgs_bar.get_beta();
    eval_f_W(
      *daeModel_
      ,*x_dot_, *x_, t_base_, Scalar(beta*coeff_x_dot_), Scalar(beta*coeff_x_)
      ,outArgs_bar.get_f().get(), &*W
      );
#ifdef THYRA_RYTHMOS_DEBUG
    std::cout << "f = "                 << std::endl;
    std::cout << *(outArgs_bar.get_f()) << std::endl;
    std::cout << "W = "                 << std::endl;
    std::cout << *W                     << std::endl;
#endif // THYRA_RYTHMOS_DEBUG
  }
  else {
    // Compute only the residual
    eval_f( *daeModel_, *x_dot_, *x_, t_base_, &*outArgs_bar.get_f() );
#ifdef THYRA_RYTHMOS_DEBUG
    std::cout << "f = "                 << std::endl;
    std::cout << *(outArgs_bar.get_f()) << std::endl;
#endif // THYRA_RYTHMOS_DEBUG
  }
#ifdef THYRA_RYTHMOS_DEBUG
  std::cout << "----------------------------------------------------------------------" << std::endl;
#endif // THYRA_RYTHMOS_DEBUG
}

} // namespace Thyra

#endif // THYRA_SINGLE_RESID_SS_DAE_MODEL_EVALUATOR_HPP
