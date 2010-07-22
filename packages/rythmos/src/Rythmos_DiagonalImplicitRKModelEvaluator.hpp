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


#ifndef RYTHMOS_DIAGONALIMPLICITRK_MODEL_EVALUATOR_HPP
#define RYTHMOS_DIAGONALIMPLICITRK_MODEL_EVALUATOR_HPP


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
class DiagonalImplicitRKModelEvaluator : virtual public Thyra::StateFuncModelEvaluatorBase<Scalar>
{
public:

  /** \name Constructors/initializers/accessors */
  //@{

  /** \brief . */
  DiagonalImplicitRKModelEvaluator();

  /** \brief . */
  void initializeDIRKModel(
    const RCP<const Thyra::ModelEvaluator<Scalar> >& daeModel,
    const Thyra::ModelEvaluatorBase::InArgs<Scalar>& basePoint,
    const RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> >& dirk_W_factory,
    const RCP<const RKButcherTableauBase<Scalar> >& irkButcherTableau
    );

  /** \brief . */
  void setTimeStepPoint(
    const RCP<const Thyra::VectorBase<Scalar> > &x_old,
    const Scalar &t_old,
    const Scalar &delta_t
    );

  /** \brief . */
  void setStageSolution(
      int stageNumber,
      const Thyra::VectorBase<Scalar>& stage_solution
      );

  /** \brief . */
  void setCurrentStage(
      int currentStage
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
  RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> > dirk_W_factory_;
  RCP<const RKButcherTableauBase<Scalar> > dirkButcherTableau_;

  Thyra::ModelEvaluatorBase::InArgs<Scalar> inArgs_;
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> outArgs_;
  Thyra::ModelEvaluatorBase::InArgs<Scalar> nominalValues_;

  RCP<const Thyra::ProductVectorSpaceBase<Scalar> > stage_space_;
  RCP<Thyra::ProductVectorBase<Scalar> > stage_derivatives_;

  bool setTimeStepPointCalled_;
  RCP<const Thyra::VectorBase<Scalar> > x_old_;
  Scalar t_old_;
  Scalar delta_t_;
  int currentStage_;

  bool isInitialized_;

};


/** \brief Nonmember constructor function.
 *
 * \relates DiagonalImplicitRKModelEvaluator
 */
template<class Scalar>
RCP<DiagonalImplicitRKModelEvaluator<Scalar> >
diagonalImplicitRKModelEvaluator(
  const RCP<const Thyra::ModelEvaluator<Scalar> >& daeModel,
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &basePoint,
  const RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> > &dirk_W_factory,
  const RCP<const RKButcherTableauBase<Scalar> > &irkButcherTableau
  )
{
  RCP<DiagonalImplicitRKModelEvaluator<Scalar> >
    dirkModel = rcp(new DiagonalImplicitRKModelEvaluator<Scalar>());
  dirkModel->initializeDIRKModel(daeModel,basePoint,dirk_W_factory,irkButcherTableau);
  return dirkModel;
}


// ///////////////////////
// Definition


// Constructors/initializers/accessors


template<class Scalar>
DiagonalImplicitRKModelEvaluator<Scalar>::DiagonalImplicitRKModelEvaluator()
{
  setTimeStepPointCalled_ = false;
  isInitialized_ = false;
  currentStage_ = -1;
}


// Overridden from ImplicitRKModelEvaluatorBase


template<class Scalar>
void DiagonalImplicitRKModelEvaluator<Scalar>::initializeDIRKModel(
  const RCP<const Thyra::ModelEvaluator<Scalar> >& daeModel,
  const Thyra::ModelEvaluatorBase::InArgs<Scalar>& basePoint,
  const RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> >& dirk_W_factory,
  const RCP<const RKButcherTableauBase<Scalar> >& irkButcherTableau
  )
{
  typedef ScalarTraits<Scalar> ST;
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
  //TEST_FOR_EXCEPT(is_null(dirk_W_factory));

  daeModel_ = daeModel;
  basePoint_ = basePoint;
  dirk_W_factory_ = dirk_W_factory;
  dirkButcherTableau_ = irkButcherTableau;

  const int numStages = dirkButcherTableau_->numStages();

  using Teuchos::rcp_dynamic_cast;
  stage_space_ = productVectorSpace(daeModel_->get_x_space(),numStages);
  RCP<const Thyra::VectorSpaceBase<Scalar> > vs = rcp_dynamic_cast<const Thyra::VectorSpaceBase<Scalar> >(stage_space_,true);
  stage_derivatives_ = rcp_dynamic_cast<Thyra::ProductVectorBase<Scalar> >(createMember(vs),true);
  Thyra::V_S( rcp_dynamic_cast<Thyra::VectorBase<Scalar> >(stage_derivatives_).ptr(),ST::zero());

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
void DiagonalImplicitRKModelEvaluator<Scalar>::setTimeStepPoint(
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

template<class Scalar>
void DiagonalImplicitRKModelEvaluator<Scalar>::setStageSolution(
      int stageNumber,
      const Thyra::VectorBase<Scalar>& stage_solution
      )
{
  TEST_FOR_EXCEPT(stageNumber < 0);
  TEST_FOR_EXCEPT(stageNumber >= dirkButcherTableau_->numStages());
  Thyra::V_V(stage_derivatives_->getNonconstVectorBlock(stageNumber).ptr(),stage_solution);
}

template<class Scalar>
void DiagonalImplicitRKModelEvaluator<Scalar>::setCurrentStage(
      int currentStage
      )
{
  TEST_FOR_EXCEPT(currentStage < 0);
  TEST_FOR_EXCEPT(currentStage >= dirkButcherTableau_->numStages());
  currentStage_ = currentStage;
}



// Overridden from ModelEvaluator


template<class Scalar>
RCP<const Thyra::VectorSpaceBase<Scalar> >
DiagonalImplicitRKModelEvaluator<Scalar>::get_x_space() const
{
  TEST_FOR_EXCEPTION( !isInitialized_, std::logic_error,
      "Error, initializeDIRKModel must be called first!\n"
      );
  return daeModel_->get_x_space();
}


template<class Scalar>
RCP<const Thyra::VectorSpaceBase<Scalar> >
DiagonalImplicitRKModelEvaluator<Scalar>::get_f_space() const
{
  TEST_FOR_EXCEPTION( !isInitialized_, std::logic_error,
      "Error, initializeDIRKModel must be called first!\n"
      );
  return daeModel_->get_f_space();
}


template<class Scalar>
RCP<Thyra::LinearOpBase<Scalar> >
DiagonalImplicitRKModelEvaluator<Scalar>::create_W_op() const
{
  TEST_FOR_EXCEPTION( !isInitialized_, std::logic_error,
      "Error, initializeDIRKModel must be called first!\n"
      );
  return daeModel_->create_W_op();
}


template<class Scalar>
RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> >
DiagonalImplicitRKModelEvaluator<Scalar>::get_W_factory() const
{
  TEST_FOR_EXCEPTION( !isInitialized_, std::logic_error,
      "Error, initializeDIRKModel must be called first!\n"
      );
  return daeModel_->get_W_factory();
}  


template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
DiagonalImplicitRKModelEvaluator<Scalar>::getNominalValues() const
{
  TEST_FOR_EXCEPTION( !isInitialized_, std::logic_error,
      "Error, initializeDIRKModel must be called first!\n"
      );
  return nominalValues_;
}


template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
DiagonalImplicitRKModelEvaluator<Scalar>::createInArgs() const
{
  TEST_FOR_EXCEPTION( !isInitialized_, std::logic_error,
      "Error, initializeDIRKModel must be called first!\n"
      );
  return inArgs_;
}


// Private functions overridden from ModelEvaluatorDefaultBase


template<class Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
DiagonalImplicitRKModelEvaluator<Scalar>::createOutArgsImpl() const
{
  TEST_FOR_EXCEPTION( !isInitialized_, std::logic_error,
      "Error, initializeDIRKModel must be called first!\n"
      );
  return outArgs_;
}


template<class Scalar>
void DiagonalImplicitRKModelEvaluator<Scalar>::evalModelImpl(
  const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs_stage,
  const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs_stage
  ) const
{

  typedef ScalarTraits<Scalar> ST;
  typedef Thyra::ModelEvaluatorBase MEB;

  TEST_FOR_EXCEPTION( !isInitialized_, std::logic_error,
      "Error!  initializeDIRKModel must be called before evalModel\n"
      );

  TEST_FOR_EXCEPTION( !setTimeStepPointCalled_, std::logic_error,
      "Error!  setTimeStepPoint must be called before evalModel"
      );

  TEST_FOR_EXCEPTION( currentStage_ == -1, std::logic_error,
      "Error!  setCurrentStage must be called before evalModel"
      );

  THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_GEN_BEGIN(
    "Rythmos::DiagonalImplicitRKModelEvaluator",inArgs_stage,outArgs_stage,daeModel_
    );

  //
  // A) Unwrap the inArgs and outArgs 
  //

  const RCP<const Thyra::VectorBase<Scalar> > x_in = inArgs_stage.get_x();
  const RCP<Thyra::VectorBase<Scalar> > f_out = outArgs_stage.get_f();
  const RCP<Thyra::LinearOpBase<Scalar> > W_op_out = outArgs_stage.get_W_op();

  //
  // B) Assemble f_out and W_op_out for given stage
  //

  MEB::InArgs<Scalar> daeInArgs = daeModel_->createInArgs();
  MEB::OutArgs<Scalar> daeOutArgs = daeModel_->createOutArgs();
  const RCP<Thyra::VectorBase<Scalar> > x_i = createMember(daeModel_->get_x_space());
  daeInArgs.setArgs(basePoint_);
  
  // B.1) Setup the DAE's inArgs for stage f(currentStage_) ...
  V_V(stage_derivatives_->getNonconstVectorBlock(currentStage_).ptr(),*x_in);
  assembleIRKState( currentStage_, dirkButcherTableau_->A(), delta_t_, *x_old_, *stage_derivatives_, outArg(*x_i) );
  daeInArgs.set_x( x_i );
  daeInArgs.set_x_dot( x_in );
  daeInArgs.set_t( t_old_ + dirkButcherTableau_->c()(currentStage_) * delta_t_ );
  daeInArgs.set_alpha(ST::one());
  daeInArgs.set_beta( delta_t_ * dirkButcherTableau_->A()(currentStage_,currentStage_) );

  // B.2) Setup the DAE's outArgs for stage f(i) ...
  if (!is_null(f_out))
    daeOutArgs.set_f( f_out );
  if (!is_null(W_op_out))
    daeOutArgs.set_W_op(W_op_out);

  // B.3) Compute f_out(i) and/or W_op_out ...
  daeModel_->evalModel( daeInArgs, daeOutArgs );
  daeOutArgs.set_f(Teuchos::null);
  daeOutArgs.set_W_op(Teuchos::null);
  
  THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_END();
  
}


} // namespace Rythmos


#endif // RYTHMOS_DIAGONALIMPLICITRK_MODEL_EVALUATOR_HPP
