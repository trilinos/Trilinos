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


#ifndef RYTHMOS_IMPLICITRK_MODEL_EVALUATOR_HPP
#define RYTHMOS_IMPLICITRK_MODEL_EVALUATOR_HPP


#include "Rythmos_Types.hpp"
#include "Rythmos_RKButcherTableau.hpp"
#include "Thyra_StateFuncModelEvaluatorBase.hpp"
#include "Thyra_ModelEvaluatorHelpers.hpp"
#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_DefaultBlockedTriangularLinearOpWithSolveFactory.hpp" // Default implementation
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_TestingTools.hpp"
#include "Teuchos_as.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_Assert.hpp"


namespace Rythmos {


/** \brief . */
template<class Scalar>
class ImplicitRKModelEvaluator : virtual public Thyra::StateFuncModelEvaluatorBase<Scalar>
{
public:

  /** \name Constructors/initializers/accessors */
  //@{

  /** \brief . */
  ImplicitRKModelEvaluator();

  /** \brief . */
  void initializeIRKModel(
    const RCP<const Thyra::ModelEvaluator<Scalar> > &daeModel,
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &basePoint,
    const RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> > &irk_W_factory,
    const RKButcherTableau<Scalar> &irkButcherTableau
    );

  /** \brief . */
  void setTimeStepPoint(
    const RCP<const Thyra::VectorBase<Scalar> > &x_old,
    const Scalar &t_old,
    const Scalar &delta_t
    );

  //@}

  /** \name Public functions overridden from ModelEvaluator */
  //@{

  /** \brief . */
  RCP<const Thyra::VectorSpaceBase<Scalar> > get_x_space() const;
  /** \brief . */
  RCP<const Thyra::VectorSpaceBase<Scalar> > get_f_space() const;
  /** \brief . */
  RCP<Thyra::LinearOpBase<Scalar> > create_W_op() const;
  /** \brief . */
  RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> > get_W_factory() const;
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

  RCP<const Thyra::ModelEvaluator<Scalar> > daeModel_;
  Thyra::ModelEvaluatorBase::InArgs<Scalar> basePoint_;
  RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> > irk_W_factory_;
  RKButcherTableau<Scalar> irkButcherTableau_;

  RCP<const Thyra::ProductVectorSpaceBase<Scalar> > x_bar_space_;
  RCP<const Thyra::ProductVectorSpaceBase<Scalar> > f_bar_space_;
  RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> > W_bar_factory_;

  Thyra::ModelEvaluatorBase::InArgs<Scalar> nominalValues_;

  RCP<const Thyra::VectorBase<Scalar> > x_old_;
  Scalar t_old_;
  Scalar delta_t_;

};


/** \brief Nonmember constructor function.
 *
 * \relates ???
 */
template<class Scalar>
RCP<ImplicitRKModelEvaluator<Scalar> >
implicitRKModelEvaluator(
  const RCP<const Thyra::ModelEvaluator<Scalar> > &daeModel,
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &basePoint,
  const RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> > &irk_W_factory,
  const RKButcherTableau<Scalar> &irkButcherTableau
  )
{
  RCP<ImplicitRKModelEvaluator<Scalar> >
    irkModel = Teuchos::rcp(new ImplicitRKModelEvaluator<Scalar>());
  irkModel->initializeIRKModel(daeModel,basePoint,irk_W_factory,irkButcherTableau);
  return irkModel;
}




// ///////////////////////
// Definition


// Constructors/initializers/accessors


template<class Scalar>
ImplicitRKModelEvaluator<Scalar>::ImplicitRKModelEvaluator()
{}


// Overridden from ImplicitRKModelEvaluatorBase


template<class Scalar>
void ImplicitRKModelEvaluator<Scalar>::initializeIRKModel(
  const RCP<const Thyra::ModelEvaluator<Scalar> > &daeModel,
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &basePoint,
  const RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> > &irk_W_factory,
  const RKButcherTableau<Scalar> &irkButcherTableau
  )
{

  // ToDo: Assert input arguments!

  daeModel_ = daeModel;
  basePoint_ = basePoint;
  irk_W_factory_ = irk_W_factory;
  irkButcherTableau_ = irkButcherTableau;

  const int numStages = irkButcherTableau_.numStages();

  x_bar_space_ = productVectorSpace(daeModel_->get_x_space(),numStages);
  f_bar_space_ = productVectorSpace(daeModel_->get_f_space(),numStages);

  // ToDo: create the preconditioner factory and set this on irk_W_factory_!
  
}


template<class Scalar>
void ImplicitRKModelEvaluator<Scalar>::setTimeStepPoint(
  const RCP<const Thyra::VectorBase<Scalar> > &x_old,
  const Scalar &t_old,
  const Scalar &delta_t
  )
{
  typedef ScalarTraits<Scalar> ST;
  TEST_FOR_EXCEPT( is_null(x_old) );
  TEST_FOR_EXCEPT( delta_t <= ST::zero() );
  // ToDo: Validate input
  x_old_ = x_old;
  t_old_ = t_old;
  delta_t_ = delta_t;
}


// Overridden from ModelEvaluator


template<class Scalar>
RCP<const Thyra::VectorSpaceBase<Scalar> >
ImplicitRKModelEvaluator<Scalar>::get_x_space() const
{
  return x_bar_space_;
}


template<class Scalar>
RCP<const Thyra::VectorSpaceBase<Scalar> >
ImplicitRKModelEvaluator<Scalar>::get_f_space() const
{
  return f_bar_space_;
}


template<class Scalar>
RCP<Thyra::LinearOpBase<Scalar> >
ImplicitRKModelEvaluator<Scalar>::create_W_op() const
{
  return daeModel_->create_W_op();
  // ToDo: Create a blocked version!
}


template<class Scalar>
RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> >
ImplicitRKModelEvaluator<Scalar>::get_W_factory() const
{
  return irk_W_factory_;
}  


template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
ImplicitRKModelEvaluator<Scalar>::getNominalValues() const
{
  return nominalValues_;
}


template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
ImplicitRKModelEvaluator<Scalar>::createInArgs() const
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
ImplicitRKModelEvaluator<Scalar>::createOutArgsImpl() const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  MEB::OutArgsSetup<Scalar> outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.setSupports(MEB::OUT_ARG_f);
  outArgs.setSupports(MEB::OUT_ARG_W_op);
  return outArgs;
}


template<class Scalar>
void ImplicitRKModelEvaluator<Scalar>::evalModelImpl(
  const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs_bar,
  const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs_bar
  ) const
{

  using std::endl;
  using Teuchos::as;
  using Teuchos::rcp_dynamic_cast;
  typedef ScalarTraits<Scalar> ST;
  typedef Thyra::ModelEvaluatorBase MEB;

  THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_GEN_BEGIN(
    "Rythmos::ImplicitRKModelEvaluator",inArgs_bar,outArgs_bar,daeModel_
    );

  //
  // A) Unwrap the inArgs
  //

  const RCP<const Thyra::ProductVectorBase<Scalar> >
    x_bar = rcp_dynamic_cast<const Thyra::ProductVectorBase<Scalar> >(
      inArgs_bar.get_x(), true
      );

  //
  // B) Unwrap the outArgs
  //

  const RCP<Thyra::ProductVectorBase<Scalar> >
    f_bar = rcp_dynamic_cast<Thyra::ProductVectorBase<Scalar> >(
      outArgs_bar.get_f(), true
      );
  
  // ToDo: Change this to a block operator!
  RCP<Thyra::LinearOpBase<Scalar> >
    W_op = outArgs_bar.get_W_op();

  //
  // C) Assemble f_bar and W_op_bar
  //

  MEB::InArgs<Scalar> daeInArgs = daeModel_->createInArgs();
  MEB::OutArgs<Scalar> daeOutArgs = daeModel_->createOutArgs();
  daeInArgs.setArgs(basePoint_);
  
  const int numStages = irkButcherTableau_.numStages();
  for ( int i = 0; i < numStages; ++i ) {

    // C.1) Setup the DAE's inArgs for this stage function

    const RCP<const Thyra::VectorBase<Scalar> >
      x_dot_i = x_bar->getVectorBlock(i);

    const RCP<Thyra::VectorBase<Scalar> >
      x_i = createMember(daeModel_->get_x_space());

    assembleIRKState(
      i, irkButcherTableau_.A(), delta_t_, *x_old_, *x_bar,
      &*x_i
      );
    
    daeInArgs.set_x_dot( x_dot_i );
    daeInArgs.set_x( x_i );
    daeInArgs.set_t( t_old_ + irkButcherTableau_.c()(i) * delta_t_ );
    daeInArgs.set_alpha(ST::one());
    daeInArgs.set_beta( delta_t_ * irkButcherTableau_.A()(i,0) );

    // C.2) Setup the DAE's outArgs for this stage function

    if (!is_null(f_bar))
      daeOutArgs.set_f( f_bar->getNonconstVectorBlock(i) );

    if (!is_null(W_op))
      daeOutArgs.set_W_op(W_op);
    // ToDo: Modify this for general block W_op_bar!

    // C.3) Compute f_bar(i) and W_op_bar(i,0) ...

    daeModel_->evalModel( daeInArgs, daeOutArgs );

    daeOutArgs.set_f(Teuchos::null);
    
    // C.4) Evaluate the rest of the W_op_bar(i,j) ...

    if (numStages > 1) {
      TEST_FOR_EXCEPT(true); // Can't handle this yet
      // ToDo: Set beta to delta_t_ * A(i,j)!
    }

  }

  THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_END();

}


} // namespace Rythmos


#endif // RYTHMOS_IMPLICITRK_MODEL_EVALUATOR_HPP
