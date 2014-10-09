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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef RYTHMOS_STATE_AND_FORWARD_SENSITIVITY_MODEL_EVALUATOR_HPP
#define RYTHMOS_STATE_AND_FORWARD_SENSITIVITY_MODEL_EVALUATOR_HPP


#include "Rythmos_ForwardSensitivityModelEvaluatorBase.hpp"
#include "Thyra_ModelEvaluator.hpp" // Interface
#include "Thyra_StateFuncModelEvaluatorBase.hpp" // Implementation
#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_ModelEvaluatorDelegatorBase.hpp"
#include "Thyra_DefaultMultiVectorProductVectorSpace.hpp"
#include "Thyra_DefaultMultiVectorProductVector.hpp"
#include "Thyra_DefaultMultiVectorLinearOpWithSolve.hpp"
#include "Teuchos_implicit_cast.hpp"


namespace Rythmos {


/** \brief Combined State and Forward sensitivity transient
 * <tt>ModelEvaluator</tt> subclass.
 *
 * This class provides an implemenation of a combined state and forward
 * sensitivity model evaluator for a DAE.
 *
 * The form of the parameterized state equation is:

 \verbatim

   f(x_dot(t),x(t),p) = 0, over t = [t0,tf]

   x(t0) = x_init(p)

 \endverbatim

 * The forward sensitivity equations, as written in multi-vector form, are:

 \verbatim

   d(f)/d(x_dot)*S_dot + d(f)/d(x)*S + d(f)/d(p) = 0, over t = [t0,tf]

   S(t0) = d(x_init)/d(p)

 \endverbatim

 * where <tt>S</tt> is a multi-vector with <tt>np</tt> columns where each
 * column <tt>S(:,j) = d(x)/d(p_j)</tt> is the sensitivity of <tt>x(t)</tt>
 * with respect to the <tt>p_j</tt> parameter.
 *
 * This model evaluator class represents the full state plus forward
 * sensitivity system given as:

 \verbatim

   f_bar(x_bar_dot(t),x_bar(t)) = 0, over t = [t0,tf]

   x_bar(t0) = x_bar_init

 \endverbatim

 * where

 \verbatim

   x_bar = [ x; s_bar ] 

   s_bar = [ S(:,0); S(:,0); ...; S(:,np-1) ]

 \endverbatim

 * and <tt>f_bar(...)</tt> is the obvious concatenated state and sensitivity
 * systems expresses in product vector form.
 *
 * The vector <tt>x_bar</tt> is represented as a
 * <tt>Thyra::ProductVectorBase</tt> object with two vector blocks <tt>x</tt>
 * and <tt>s_bar</tt>.  The flattened out sensitivity vector <tt>s_bar</tt> is
 * then represented as a specialized product vector of type
 * <tt>Thyra::DefaultMultiVectorProductVector</tt>.
 *
 * If <tt>x_bar</tt> is an <tt>RCP<Thyra::VectorBase<Scalar> ></tt> object, then
 * <tt>x</tt> can be access as

 \code

   RCP<Thyra::VectorBase<Scalar> >
     x = Thyra::nonconstProductVectorBase<Scalar>(x_bar)->getVectorBlock(0);

 \endcode

 * If <tt>x_bar</tt> is an <tt>RCP<const Thyra::VectorBase<Scalar> ></tt>
 * object, then <tt>x</tt> can be access as

 \code

   RCP<const Thyra::VectorBase<Scalar> >
     x = Thyra::productVectorBase<Scalar>(x_bar)->getVectorBlock(0);

 \endcode

 * Likewise, <tt>s_bar</tt> can be access as


 \code

   RCP<Thyra::VectorBase<Scalar> >
     s_bar = Thyra::nonconstProductVectorBase<Scalar>(x_bar)->getVectorBlock(1);

 \endcode

 * when non-const and when const as:

 \code

   RCP<const Thyra::VectorBase<Scalar> >
     s_bar = Thyra::productVectorBase<Scalar>(x_bar)->getVectorBlock(1);

 \endcode

 * Given the flattened out vector form <tt>s_bar</tt>, one can get the
 * underlying mulit-vector <tt>S</tt> as:

 \code
 
  typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;

   RCP<Thyra::MultiVectorBase<Scalar> >
     S = rcp_dynamic_cast<DMVPV>(s_bar)->getNonconstMultiVector();

 \endcode

 * for a nonconst vector/multi-vector and for a const vector/multi-vector as:

 \code
 
  typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;

   RCP<const Thyra::MultiVectorBase<Scalar> >
     S = rcp_dynamic_cast<const DMVPV>(s_bar)->getMultiVector();

 \endcode
 
 * ToDo: Replace the above documentation with the helper functions that will
 * do all of this!
 *
 * Currently, this class does not implement the full ModelEvaluator interface
 * and it really just provides the spaces for x_bar and f_bar and the InArgs
 * and OutArgs creation functions to allow for the full specification of the
 * <tt>ForwardSensitivityStepper</tt> class.  This is especially important in
 * order to correctly set the full initial condition for the state and the
 * forward sensitivities.  Later this class can be completely finished in
 * which case it would necessarly implement the evaluations and the linear
 * solves which would automatically support the simultaneous corrector method.
 *
 * ToDo: Finish documentation!
 */
template<class Scalar>
class StateAndForwardSensitivityModelEvaluator
  : virtual public Thyra::StateFuncModelEvaluatorBase<Scalar>
{
public:

  /** \name Constructors/Intializers/Accessors */
  //@{

  /** \brief . */
  StateAndForwardSensitivityModelEvaluator();

  /** \brief Set up the structure of the state and sensitivity equations.
   *
   * \param  sensModel
   *           [in,persisting] The structure-initialized forward sensitivity
   *           model.  From this object, the state model and other information
   *           will be extracted.
   */
  void initializeStructure(
    const Teuchos::RCP<const ForwardSensitivityModelEvaluatorBase<Scalar> > &sensModel
    );

  // 2007/05/30: rabartl: ToDo: Add function to set the nominal values etc.

  /** \brief Create a wrapped product vector of the form <tt>x_bar = [ x; s_bar ]</tt>.
   *
   * Note: This does not copy any vector data, it only creates the wrapped
   * product vector.
   */
  Teuchos::RCP<const Thyra::DefaultProductVector<Scalar> >
  create_x_bar_vec(
    const Teuchos::RCP<const Thyra::VectorBase<Scalar> > &x_vec,
    const Teuchos::RCP<const Thyra::VectorBase<Scalar> > &s_bar_vec
    ) const;

  //@}

  /** \name Public functions overridden from ModelEvaulator. */
  //@{

  /** \brief Returns 0 . */
  int Np() const;
  /** \brief . */
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_p_space(int l) const;
  /** \brief . */
  Teuchos::RCP<const Teuchos::Array<std::string> > get_p_names(int l) const;
  /** \brief . */
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_x_space() const;
  /** \brief . */
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_f_space() const;
  /** \brief . */
  Thyra::ModelEvaluatorBase::InArgs<Scalar> getNominalValues() const;
  /** \brief . */
  Teuchos::RCP<Thyra::LinearOpWithSolveBase<Scalar> > create_W() const;
  /** \brief . */
  Thyra::ModelEvaluatorBase::InArgs<Scalar> createInArgs() const;

  //@}

private:

  /** \name Private functions overridden from ModelEvaulatorDefaultBase. */
  //@{

  /** \brief . */
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> createOutArgsImpl() const;
  /** \brief . */
  void evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs
    ) const;

  //@}

private:

  // /////////////////////////
  // Private data members

  Teuchos::RCP<const ForwardSensitivityModelEvaluatorBase<Scalar> > sensModel_;

  int Np_;
  Teuchos::RCP<const Thyra::DefaultProductVectorSpace<Scalar> > x_bar_space_;
  Teuchos::RCP<const Thyra::DefaultProductVectorSpace<Scalar> > f_bar_space_;
  
};


// /////////////////////////////////
// Implementations


// Constructors/Intializers/Accessors


template<class Scalar>
StateAndForwardSensitivityModelEvaluator<Scalar>::StateAndForwardSensitivityModelEvaluator()
  :Np_(0)
{}


template<class Scalar>
void StateAndForwardSensitivityModelEvaluator<Scalar>::initializeStructure(
  const Teuchos::RCP<const ForwardSensitivityModelEvaluatorBase<Scalar> > &sensModel
  )
{

  using Teuchos::tuple; using Teuchos::RCP;

  TEUCHOS_TEST_FOR_EXCEPT( is_null(sensModel) );
  
  sensModel_ = sensModel;

  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >
    stateModel = sensModel_->getStateModel();
  
  x_bar_space_ = Thyra::productVectorSpace<Scalar>(
    tuple<RCP<const Thyra::VectorSpaceBase<Scalar> > >(
      stateModel->get_x_space(), sensModel_->get_x_space()
      )
    );

  f_bar_space_ = Thyra::productVectorSpace<Scalar>(
    tuple<RCP<const Thyra::VectorSpaceBase<Scalar> > >(
      stateModel->get_f_space(), sensModel_->get_f_space()
      )
    );

  Np_ = stateModel->Np();

}


template<class Scalar> 
Teuchos::RCP<const Thyra::DefaultProductVector<Scalar> >
StateAndForwardSensitivityModelEvaluator<Scalar>::create_x_bar_vec(
  const Teuchos::RCP<const Thyra::VectorBase<Scalar> > &x_vec,
  const Teuchos::RCP<const Thyra::VectorBase<Scalar> > &s_bar_vec
  ) const
{

  using Teuchos::tuple;
  using Teuchos::RCP;
  typedef RCP<const Thyra::VectorBase<Scalar> > RCPCV;

  return Thyra::defaultProductVector<Scalar>(
    x_bar_space_, tuple<RCPCV>(x_vec,s_bar_vec)
    );

}


// Public functions overridden from ModelEvaulator


template<class Scalar>
int StateAndForwardSensitivityModelEvaluator<Scalar>::Np() const
{
  return Np_;
}


template<class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
StateAndForwardSensitivityModelEvaluator<Scalar>::get_p_space(int l) const
{
  return sensModel_->getStateModel()->get_p_space(l);
}


template<class Scalar>
Teuchos::RCP<const Teuchos::Array<std::string> >
StateAndForwardSensitivityModelEvaluator<Scalar>::get_p_names(int l) const
{
  return sensModel_->getStateModel()->get_p_names(l);
}


template<class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
StateAndForwardSensitivityModelEvaluator<Scalar>::get_x_space() const
{
  return x_bar_space_;
}


template<class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
StateAndForwardSensitivityModelEvaluator<Scalar>::get_f_space() const
{
  return f_bar_space_;
}


template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
StateAndForwardSensitivityModelEvaluator<Scalar>::getNominalValues() const
{
  return this->createInArgs();
}


template<class Scalar>
Teuchos::RCP<Thyra::LinearOpWithSolveBase<Scalar> >
StateAndForwardSensitivityModelEvaluator<Scalar>::create_W() const
{
  TEUCHOS_TEST_FOR_EXCEPT("ToDo: Implement create_W() when needed!");
  return Teuchos::null;
}


template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
StateAndForwardSensitivityModelEvaluator<Scalar>::createInArgs() const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  MEB::InArgs<Scalar>
    stateModelInArgs = sensModel_->getStateModel()->createInArgs();
  MEB::InArgsSetup<Scalar> inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.set_Np(Np_);
  inArgs.setSupports( MEB::IN_ARG_x_dot,
    stateModelInArgs.supports(MEB::IN_ARG_x_dot) );
  inArgs.setSupports( MEB::IN_ARG_x );
  inArgs.setSupports( MEB::IN_ARG_t );
  inArgs.setSupports( MEB::IN_ARG_alpha,
    stateModelInArgs.supports(MEB::IN_ARG_alpha) );
  inArgs.setSupports( MEB::IN_ARG_beta,
    stateModelInArgs.supports(MEB::IN_ARG_beta) );
  return inArgs;
}


// Private functions overridden from ModelEvaulatorDefaultBase


template<class Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
StateAndForwardSensitivityModelEvaluator<Scalar>::createOutArgsImpl() const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  MEB::OutArgs<Scalar>
    stateModelOutArgs = sensModel_->getStateModel()->createOutArgs();
  MEB::OutArgsSetup<Scalar> outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.set_Np_Ng(Np_,0);
  outArgs.setSupports(MEB::OUT_ARG_f);
  if (stateModelOutArgs.supports(MEB::OUT_ARG_W) ) {
    outArgs.setSupports(MEB::OUT_ARG_W);
    outArgs.set_W_properties(stateModelOutArgs.get_W_properties());
  }
  return outArgs;
}


template<class Scalar>
void StateAndForwardSensitivityModelEvaluator<Scalar>::evalModelImpl(
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
  const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs
  ) const
{
  TEUCHOS_TEST_FOR_EXCEPT("ToDo: Implement evalModel(...) when needed!");
}


} // namespace Rythmos


#endif // RYTHMOS_STATE_AND_FORWARD_SENSITIVITY_MODEL_EVALUATOR_HPP
