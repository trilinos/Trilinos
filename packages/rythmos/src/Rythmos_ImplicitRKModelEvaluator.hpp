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
#include "Rythmos_RKButcherTableauHelpers.hpp"
#include "Thyra_StateFuncModelEvaluatorBase.hpp"
#include "Thyra_ModelEvaluatorHelpers.hpp"
#include "Thyra_ModelEvaluatorDelegatorBase.hpp"
#include "Thyra_DefaultProductVectorSpace.hpp"
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
    const RCP<const Thyra::ModelEvaluator<Scalar> >& daeModel,
    const Thyra::ModelEvaluatorBase::InArgs<Scalar>& basePoint,
    const RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> >& irk_W_factory,
    const RCP<const RKButcherTableauBase<Scalar> >& irkButcherTableau
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
  RCP<const RKButcherTableauBase<Scalar> > irkButcherTableau_;

  Thyra::ModelEvaluatorBase::InArgs<Scalar> inArgs_;
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> outArgs_;
  Thyra::ModelEvaluatorBase::InArgs<Scalar> nominalValues_;

  RCP<const Thyra::ProductVectorSpaceBase<Scalar> > x_bar_space_;
  RCP<const Thyra::ProductVectorSpaceBase<Scalar> > f_bar_space_;
  RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> > W_bar_factory_;

  bool setTimeStepPointCalled_;
  RCP<const Thyra::VectorBase<Scalar> > x_old_;
  Scalar t_old_;
  Scalar delta_t_;

  bool isInitialized_;

};


/** \brief Nonmember constructor function.
 *
 * \relates ImplicitRKModelEvaluator
 */
template<class Scalar>
RCP<ImplicitRKModelEvaluator<Scalar> >
implicitRKModelEvaluator(
  const RCP<const Thyra::ModelEvaluator<Scalar> >& daeModel,
  const Thyra::ModelEvaluatorBase::InArgs<Scalar>& basePoint,
  const RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> >& irk_W_factory,
  const RCP<const RKButcherTableauBase<Scalar> >& irkButcherTableau
  )
{
  RCP<ImplicitRKModelEvaluator<Scalar> >
    irkModel = rcp(new ImplicitRKModelEvaluator<Scalar>());
  irkModel->initializeIRKModel(daeModel,basePoint,irk_W_factory,irkButcherTableau);
  return irkModel;
}


// ///////////////////////
// Definition


// Constructors/initializers/accessors


template<class Scalar>
ImplicitRKModelEvaluator<Scalar>::ImplicitRKModelEvaluator()
{
  setTimeStepPointCalled_ = false;
  isInitialized_ = false;
}


// Overridden from ImplicitRKModelEvaluatorBase


template<class Scalar>
void ImplicitRKModelEvaluator<Scalar>::initializeIRKModel(
  const RCP<const Thyra::ModelEvaluator<Scalar> >& daeModel,
  const Thyra::ModelEvaluatorBase::InArgs<Scalar>& basePoint,
  const RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> >& irk_W_factory,
  const RCP<const RKButcherTableauBase<Scalar> >& irkButcherTableau
  )
{
  // ToDo: Assert input arguments!
  // How do I verify the basePoint is an authentic InArgs from daeModel?
  TEST_FOR_EXCEPTION( 
      is_null(basePoint.get_x()), 
      std::logic_error,
      "Error!  The basepoint x vector is null!"
      );
  TEST_FOR_EXCEPTION( 
      is_null(daeModel), 
      std::logic_error,
      "Error!  The model evaluator pointer is null!"
      );
  TEST_FOR_EXCEPTION( 
      !daeModel->get_x_space()->isCompatible(*(basePoint.get_x()->space())), 
      std::logic_error,
      "Error!  The basepoint input arguments are incompatible with the model evaluator vector space!"
      );
  TEST_FOR_EXCEPT(is_null(irk_W_factory));

  daeModel_ = daeModel;
  basePoint_ = basePoint;
  irk_W_factory_ = irk_W_factory;
  irkButcherTableau_ = irkButcherTableau;

  const int numStages = irkButcherTableau_->numStages();

  x_bar_space_ = productVectorSpace(daeModel_->get_x_space(),numStages);
  f_bar_space_ = productVectorSpace(daeModel_->get_f_space(),numStages);

  // HACK! Remove the preconditioner factory for now!
  if (irk_W_factory_->acceptsPreconditionerFactory())
    irk_W_factory_->unsetPreconditionerFactory();

  // ToDo: create the block diagonal preconditioner factory and set this on
  // irk_W_factory_!
  
  // Set up prototypical InArgs
  {
    typedef Thyra::ModelEvaluatorBase MEB;
    MEB::InArgsSetup<Scalar> inArgs;
    inArgs.setModelEvalDescription(this->description());
    inArgs.setSupports(MEB::IN_ARG_x);
    inArgs_ = inArgs;
  }
  // Set up prototypical OutArgs
  {
    typedef Thyra::ModelEvaluatorBase MEB;
    MEB::OutArgsSetup<Scalar> outArgs;
    outArgs.setModelEvalDescription(this->description());
    outArgs.setSupports(MEB::OUT_ARG_f);
    outArgs.setSupports(MEB::OUT_ARG_W_op);
    outArgs_ = outArgs;
  }
  // Set up nominal values
  nominalValues_ = inArgs_;

  isInitialized_ = true;
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
  // Verify x_old is compatible with the vector space in the DAE Model.
  TEST_FOR_EXCEPTION(!daeModel_->get_x_space()->isCompatible(*(x_old->space())), std::logic_error,
      "Error!  The incoming VectorBase object is not compatible with the DAE model that was provided at initialization!");
  x_old_ = x_old;
  t_old_ = t_old;
  delta_t_ = delta_t;
  setTimeStepPointCalled_ = true;
}


// Overridden from ModelEvaluator


template<class Scalar>
RCP<const Thyra::VectorSpaceBase<Scalar> >
ImplicitRKModelEvaluator<Scalar>::get_x_space() const
{
  TEST_FOR_EXCEPTION( !isInitialized_, std::logic_error,
      "Error, initializeIRKModel must be called first!\n"
      );
  return x_bar_space_;
}


template<class Scalar>
RCP<const Thyra::VectorSpaceBase<Scalar> >
ImplicitRKModelEvaluator<Scalar>::get_f_space() const
{
  TEST_FOR_EXCEPTION( !isInitialized_, std::logic_error,
      "Error, initializeIRKModel must be called first!\n"
      );
  return f_bar_space_;
}


template<class Scalar>
RCP<Thyra::LinearOpBase<Scalar> >
ImplicitRKModelEvaluator<Scalar>::create_W_op() const
{
  TEST_FOR_EXCEPTION( !isInitialized_, std::logic_error,
      "Error, initializeIRKModel must be called first!\n"
      );
  // Create the block structure for W_op_bar right away!
  const int numStages = irkButcherTableau_->numStages();
  RCP<Thyra::PhysicallyBlockedLinearOpBase<Scalar> >
    W_op_bar = Thyra::defaultBlockedLinearOp<Scalar>();
  W_op_bar->beginBlockFill( f_bar_space_, x_bar_space_ );
  for ( int i = 0; i < numStages; ++i )
    for ( int j = 0; j < numStages; ++j )
      W_op_bar->setNonconstBlock( i, j, daeModel_->create_W_op() );
  W_op_bar->endBlockFill();
  return W_op_bar;
}


template<class Scalar>
RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> >
ImplicitRKModelEvaluator<Scalar>::get_W_factory() const
{
  TEST_FOR_EXCEPTION( !isInitialized_, std::logic_error,
      "Error, initializeIRKModel must be called first!\n"
      );
  return irk_W_factory_;
}  


template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
ImplicitRKModelEvaluator<Scalar>::getNominalValues() const
{
  TEST_FOR_EXCEPTION( !isInitialized_, std::logic_error,
      "Error, initializeIRKModel must be called first!\n"
      );
  return nominalValues_;
}


template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
ImplicitRKModelEvaluator<Scalar>::createInArgs() const
{
  TEST_FOR_EXCEPTION( !isInitialized_, std::logic_error,
      "Error, initializeIRKModel must be called first!\n"
      );
  return inArgs_;
}


// Private functions overridden from ModelEvaluatorDefaultBase


template<class Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
ImplicitRKModelEvaluator<Scalar>::createOutArgsImpl() const
{
  TEST_FOR_EXCEPTION( !isInitialized_, std::logic_error,
      "Error, initializeIRKModel must be called first!\n"
      );
  return outArgs_;
}


template<class Scalar>
void ImplicitRKModelEvaluator<Scalar>::evalModelImpl(
  const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs_bar,
  const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs_bar
  ) const
{

  using Teuchos::rcp_dynamic_cast;
  typedef ScalarTraits<Scalar> ST;
  typedef Thyra::ModelEvaluatorBase MEB;
  typedef Thyra::VectorBase<Scalar> VB;
  typedef Thyra::ProductVectorBase<Scalar> PVB;
  typedef Thyra::BlockedLinearOpBase<Scalar> BLWB;

  TEST_FOR_EXCEPTION( !isInitialized_, std::logic_error,
      "Error!  initializeIRKModel must be called before evalModel\n"
      );

  TEST_FOR_EXCEPTION( !setTimeStepPointCalled_, std::logic_error,
      "Error!  setTimeStepPoint must be called before evalModel"
      );

  THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_GEN_BEGIN(
    "Rythmos::ImplicitRKModelEvaluator",inArgs_bar,outArgs_bar,daeModel_
    );

  //
  // A) Unwrap the inArgs and outArgs to get at product vectors and block op
  //

  const RCP<const PVB> x_bar = rcp_dynamic_cast<const PVB>(inArgs_bar.get_x(), true);
  const RCP<PVB> f_bar = rcp_dynamic_cast<PVB>(outArgs_bar.get_f(), true);
  const RCP<BLWB> W_op_bar = rcp_dynamic_cast<BLWB>(outArgs_bar.get_W_op(), true);

  //
  // B) Assemble f_bar and W_op_bar by looping over stages
  //

  MEB::InArgs<Scalar> daeInArgs = daeModel_->createInArgs();
  MEB::OutArgs<Scalar> daeOutArgs = daeModel_->createOutArgs();
  const RCP<VB> x_i = createMember(daeModel_->get_x_space());
  daeInArgs.setArgs(basePoint_);
  
  const int numStages = irkButcherTableau_->numStages();

  for ( int i = 0; i < numStages; ++i ) {

    // B.1) Setup the DAE's inArgs for stage f(i) ...
    assembleIRKState( i, irkButcherTableau_->A(), delta_t_, *x_old_, *x_bar, outArg(*x_i) );
    daeInArgs.set_x( x_i );
    daeInArgs.set_x_dot( x_bar->getVectorBlock(i) );
    daeInArgs.set_t( t_old_ + irkButcherTableau_->c()(i) * delta_t_ );
    Scalar alpha = ST::zero();
    if (i == 0) {
      alpha = ST::one();
    } else {
      alpha = ST::zero();
    }
    Scalar beta = delta_t_ * irkButcherTableau_->A()(i,0);
    daeInArgs.set_alpha( alpha );
    daeInArgs.set_beta( beta );

    // B.2) Setup the DAE's outArgs for stage f(i) ...
    if (!is_null(f_bar))
      daeOutArgs.set_f( f_bar->getNonconstVectorBlock(i) );
    if (!is_null(W_op_bar)) {
      daeOutArgs.set_W_op(W_op_bar->getNonconstBlock(i,0));
    }

    // B.3) Compute f_bar(i) and/or W_op_bar(i,0) ...
    daeModel_->evalModel( daeInArgs, daeOutArgs );
    daeOutArgs.set_f(Teuchos::null);
    daeOutArgs.set_W_op(Teuchos::null);
    
    // B.4) Evaluate the rest of the W_op_bar(i,j=1...numStages-1) ...
    if (!is_null(W_op_bar)) {
      for ( int j = 1; j < numStages; ++j ) {
        alpha = ST::zero();
        if (i == j) {
          alpha = ST::one();
        } else {
          alpha = ST::zero();
        }
        beta = delta_t_ * irkButcherTableau_->A()(i,j);
        daeInArgs.set_alpha( alpha );
        daeInArgs.set_beta( beta );
        daeOutArgs.set_W_op(W_op_bar->getNonconstBlock(i,j));
        daeModel_->evalModel( daeInArgs, daeOutArgs );
        daeOutArgs.set_W_op(Teuchos::null);
      }
    }

  }
  
  THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_END();
  
}


} // namespace Rythmos


#endif // RYTHMOS_IMPLICITRK_MODEL_EVALUATOR_HPP
