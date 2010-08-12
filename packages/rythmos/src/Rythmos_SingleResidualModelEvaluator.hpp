//@HEADER
// ***********************************************************************
//
//                           Rythmos Package
//                 Copyright (2006) Sandia Corporation
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
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER


#ifndef RYTHMOS_SINGLE_RESIDUAL_MODEL_EVALUATOR_HPP
#define RYTHMOS_SINGLE_RESIDUAL_MODEL_EVALUATOR_HPP


#include "Rythmos_SingleResidualModelEvaluatorBase.hpp"
#include "Thyra_ModelEvaluatorDelegatorBase.hpp"
#include "Thyra_ModelEvaluatorHelpers.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_TestingTools.hpp"
#include "Teuchos_as.hpp"


namespace Rythmos {


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
class SingleResidualModelEvaluator
  : virtual public SingleResidualModelEvaluatorBase<Scalar>,
    virtual public Thyra::ModelEvaluatorDelegatorBase<Scalar>
{
public:

  /** \name Constructors/initializers/accessors */
  //@{

  /** \brief . */
  SingleResidualModelEvaluator();

  //@}

  /** \name Overridden from SingleResidualModelEvaluatorBase */
  //@{

  /** \brief . */
  void initializeSingleResidualModel(
    const RCP<const Thyra::ModelEvaluator<Scalar> > &daeModel,
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &basePoint,
    const Scalar &coeff_x_dot,
    const RCP<const Thyra::VectorBase<Scalar> > &x_dot_base,
    const Scalar &coeff_x,
    const RCP<const Thyra::VectorBase<Scalar> > &x_base,
    const Scalar &t_base,
    const RCP<const Thyra::VectorBase<Scalar> > &x_bar_init
    );

  /** \brief . */
  Scalar get_coeff_x_dot() const;

  /** \brief . */
  RCP<const Thyra::VectorBase<Scalar> >
  get_x_dot_base() const;

  /** \brief . */
  Scalar get_coeff_x() const;

  /** \brief . */
  RCP<const Thyra::VectorBase<Scalar> >
  get_x_base() const;

  /** \brief . */
  Scalar get_t_base() const;

  //@}

  /** \name Public functions overridden from ModelEvaluator */
  //@{

  /** \breif . */
  Thyra::ModelEvaluatorBase::InArgs<Scalar> getNominalValues() const;
  /** \brief . */
  Thyra::ModelEvaluatorBase::InArgs<Scalar> createInArgs() const;

  //@}

private:

  /** \name Private functions overridden from ModelEvaluatorDefaultBase */
  //@{

  /** \brief . */
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> createOutArgsImpl() const;
  /** \brief . */
  void evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs
    ) const;

  //@}


private:

  Thyra::ModelEvaluatorBase::InArgs<Scalar> basePoint_;
  Scalar coeff_x_dot_;
  RCP<const Thyra::VectorBase<Scalar> > x_dot_base_;
  Scalar coeff_x_;
  RCP<const Thyra::VectorBase<Scalar> > x_base_;
  Scalar t_base_;

  Thyra::ModelEvaluatorBase::InArgs<Scalar> nominalValues_;

  // cache
  RCP<Thyra::VectorBase<Scalar> > x_;
  RCP<Thyra::VectorBase<Scalar> > x_dot_;

};


// ///////////////////////
// Definition


// Constructors/initializers/accessors


template<class Scalar>
SingleResidualModelEvaluator<Scalar>::SingleResidualModelEvaluator()
{}


// Overridden from SingleResidualModelEvaluatorBase


template<class Scalar>
void SingleResidualModelEvaluator<Scalar>::initializeSingleResidualModel(
  const RCP<const Thyra::ModelEvaluator<Scalar> > &daeModel,
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &basePoint,
  const Scalar &coeff_x_dot,
  const RCP<const Thyra::VectorBase<Scalar> > &x_dot_base,
  const Scalar &coeff_x,
  const RCP<const Thyra::VectorBase<Scalar> > &x_base,
  const Scalar &t_base,
  const RCP<const Thyra::VectorBase<Scalar> > &x_bar_init
  )
{
  this->Thyra::ModelEvaluatorDelegatorBase<Scalar>::initialize(daeModel);
  basePoint_ = basePoint;
  coeff_x_dot_ = coeff_x_dot;
  x_dot_base_ = x_dot_base;
  coeff_x_ = coeff_x;
  x_base_ = x_base;
  t_base_ = t_base;

  nominalValues_ = daeModel->getNominalValues();
  nominalValues_.set_x(x_bar_init);

  x_dot_ = createMember( daeModel->get_x_space() );
  x_ = createMember( daeModel->get_x_space() );

  // ToDo: Check that daeModel supports x_dot, x and maybe t

#ifdef THYRA_RYTHMOS_DEBUG
  std::cout << "----------------------------------------------------------------------" << std::endl;
  std::cout << "Rythmos::SingleResidualModelEvaluator::initialize" << std::endl;
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
  std::cout << "x_bar_init = ";                  
  if ( x_bar_init.get() )
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


template<class Scalar>
Scalar SingleResidualModelEvaluator<Scalar>::get_coeff_x_dot() const
{
  return coeff_x_dot_;
}


template<class Scalar>
RCP<const Thyra::VectorBase<Scalar> >
SingleResidualModelEvaluator<Scalar>::get_x_dot_base() const
{
  return x_dot_base_;
}


template<class Scalar>
Scalar SingleResidualModelEvaluator<Scalar>::get_coeff_x() const
{
  return coeff_x_;
}


template<class Scalar>
RCP<const Thyra::VectorBase<Scalar> >
SingleResidualModelEvaluator<Scalar>::get_x_base() const
{
  return x_base_;
}


template<class Scalar>
Scalar SingleResidualModelEvaluator<Scalar>::get_t_base() const
{
  return t_base_;
}


// Overridden from ModelEvaluator


template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
SingleResidualModelEvaluator<Scalar>::getNominalValues() const
{
  return nominalValues_;
}


template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
SingleResidualModelEvaluator<Scalar>::createInArgs() const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  MEB::InArgsSetup<Scalar> inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.setSupports(MEB::IN_ARG_x);
  return inArgs;
}


// Private functions overridden from ModelEvaluatorDefaultBase


template<class Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
SingleResidualModelEvaluator<Scalar>::createOutArgsImpl() const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  MEB::OutArgsSetup<Scalar> outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.setSupports(MEB::OUT_ARG_f);
  outArgs.setSupports(MEB::OUT_ARG_W);
  return outArgs;
}


template<class Scalar>
void SingleResidualModelEvaluator<Scalar>::evalModelImpl(
  const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs_bar,
  const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs_bar
  ) const
{

  using std::endl;
  using Teuchos::as;
  typedef Thyra::ModelEvaluatorBase MEB;

  const RCP<const Thyra::ModelEvaluator<Scalar> >
    daeModel = this->getUnderlyingModel();

  THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_BEGIN(
    "Rythmos::SingleResidualModelEvaluator",inArgs_bar,outArgs_bar
    );

  const bool dumpAll = ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_EXTREME) );

  const Thyra::VectorBase<Scalar> &x_bar = *inArgs_bar.get_x(); 

  // x_dot = coeff_x_dot * x_bar + x_dot_base

  if (x_dot_base_.get())
    Thyra::V_StVpV( x_dot_.ptr(), coeff_x_dot_, x_bar, *x_dot_base_ );
  else
    Thyra::V_StV( x_dot_.ptr(), coeff_x_dot_, x_bar);

  if (dumpAll) {
    *out << "\nx_dot_ = coeff_x_dot_ * x_bar + x_dot_base_\n";
    *out << "\ncoeff_x_dot_ = " << coeff_x_dot_ << endl;
    *out << "\nx_bar = " << x_bar;
    *out << "\nx_dot_base_ = ";
    if ( x_dot_base_.get() )
      *out << *x_dot_base_;
    else
      *out << "null" << endl;
    *out << "\nx_dot_ = ";
    if ( x_dot_.get() )
      *out << *x_dot_;
    else
      *out << "null" << endl;
  }

  // x = coeff_x * x_bar + x_base

  if (x_base_.get())
    Thyra::V_StVpV( x_.ptr(), coeff_x_, x_bar, *x_base_ );
  else
    Thyra::V_StV( x_.ptr(), coeff_x_, x_bar);

  if (dumpAll) {
    *out << "\nx_ = coeff_x_ * x_bar + x_base_\n";
    *out << "\ncoeff_x_ = " << coeff_x_ << endl;
    *out << "\nx_bar = " << x_bar;
    *out << "\nx_base_ = ";
    if ( x_base_.get() )
      *out << *x_base_;
    else
      *out << "null" << endl;
    *out << "\nx_ = ";
    if ( x_.get() )
      *out << *x_;
    else
      *out << "null" << endl;
  }

  // Compute W and f

  if (as<int>(verbLevel) >= as<int>(Teuchos::VERB_LOW))
    *out << "\nEvaluating the underlying DAE model at (x_bar_dot,x_bar,t) ...\n";

  RCP<Thyra::LinearOpWithSolveBase<Scalar> > W;

  MEB::InArgs<Scalar> daeInArgs = daeModel->createInArgs();
  daeInArgs.setArgs(basePoint_);
  daeInArgs.set_x_dot(x_dot_);
  daeInArgs.set_x(x_);
  daeInArgs.set_t(t_base_);
  daeInArgs.set_alpha(coeff_x_dot_);
  daeInArgs.set_beta(coeff_x_);
  MEB::OutArgs<Scalar> daeOutArgs = daeModel->createOutArgs();
  daeOutArgs.set_f(outArgs_bar.get_f()); // can be null
  daeOutArgs.set_W(outArgs_bar.get_W()); // can be null
  daeModel->evalModel(daeInArgs,daeOutArgs);

  THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_END();

}


} // namespace Rythmos


#endif // RYTHMOS_SINGLE_RESIDUAL_MODEL_EVALUATOR_HPP
