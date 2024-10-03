// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_01.hpp
    \brief Header file to define the ModelEvaluator used in the Tempus demonstration.
*/

#include "Thyra_ModelEvaluator.hpp" // Interface
#include "Thyra_StateFuncModelEvaluatorBase.hpp" // Implementation

#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_DetachedMultiVectorView.hpp"
#include "Thyra_DefaultSerialDenseLinearOpWithSolveFactory.hpp"
#include "Thyra_DefaultMultiVectorLinearOpWithSolve.hpp"
#include "Thyra_DefaultLinearOpSource.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_DefaultMultiVectorProductVector.hpp"

/*****************
** DECLARATION. **
*****************/

template<class Scalar>
class SinCosModelEvaluator : public Thyra::StateFuncModelEvaluatorBase<Scalar> {
public:

  SinCosModelEvaluator();
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar>> get_x_space() const;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar>> get_f_space() const;
  Thyra::ModelEvaluatorBase ::InArgs<Scalar>         getNominalValues() const;
  Thyra::ModelEvaluatorBase::InArgs<Scalar>          createInArgs() const;

private:

  void setupInOutArgs_() const;

  // methods overridden from base class
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> createOutArgsImpl() const;
  void evalModelImpl(const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
                     const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const;

private: // SinCos data members

  int dim_;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar>> x_space_;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar>> f_space_;
  mutable bool isInitialized_;
  mutable Thyra::ModelEvaluatorBase::InArgs<Scalar>  inArgs_;
  mutable Thyra::ModelEvaluatorBase::OutArgs<Scalar> outArgs_;
  mutable Thyra::ModelEvaluatorBase::InArgs<Scalar>  nominalValues_;

};


/********************
** IMPLEMENTATION. **
********************/

template<class Scalar>
SinCosModelEvaluator<Scalar>::SinCosModelEvaluator() {
  isInitialized_ = false;
  dim_ = 2;
  // Create x_space and f_space
  x_space_ = Thyra::defaultSpmdVectorSpace<Scalar>(dim_);
  f_space_ = Thyra::defaultSpmdVectorSpace<Scalar>(dim_);
}


template<class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar>> SinCosModelEvaluator<Scalar>::get_x_space() const {
  return x_space_;
}


template<class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar>> SinCosModelEvaluator<Scalar>::get_f_space() const {
  return f_space_;
}


template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar> SinCosModelEvaluator<Scalar>::getNominalValues() const {
  TEUCHOS_TEST_FOR_EXCEPTION( !isInitialized_, std::logic_error,"Error, setupInOutArgs_ must be called first!\n");
  return nominalValues_;
}


template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar> SinCosModelEvaluator<Scalar>::createInArgs() const {
  setupInOutArgs_();
  return inArgs_;
}


template<class Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar> SinCosModelEvaluator<Scalar>::createOutArgsImpl() const {
  setupInOutArgs_();
  return outArgs_;
}


template<class Scalar>
void SinCosModelEvaluator<Scalar>::evalModelImpl(
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
  const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const {

  using Teuchos::RCP;
  TEUCHOS_TEST_FOR_EXCEPTION( !isInitialized_, std::logic_error,"Error, setupInOutArgs_ must be called first!\n");

  const Thyra::ConstDetachedVectorView<Scalar> x_in_view(inArgs.get_x());
  const RCP<Thyra::VectorBase<Scalar> > f_out = outArgs.get_f();

  // Evaluate the Explicit ODE f(x,t) [= 0]
  if (!is_null(f_out)) {
    Thyra::DetachedVectorView<Scalar> f_out_view( *f_out );
    f_out_view[0] =  x_in_view[1];
    f_out_view[1] = -x_in_view[0];
  }

}


template<class Scalar>
void SinCosModelEvaluator<Scalar>::setupInOutArgs_() const {
  if (isInitialized_) {
    return;
  }

  { // Set up prototypical InArgs.
    Thyra::ModelEvaluatorBase::InArgsSetup<Scalar> inArgs;
    inArgs.setModelEvalDescription(this->description());
    inArgs.setSupports( Thyra::ModelEvaluatorBase::IN_ARG_t );
    inArgs.setSupports( Thyra::ModelEvaluatorBase::IN_ARG_x );
    inArgs_ = inArgs;
  }

  { // Set up prototypical OutArgs,
    Thyra::ModelEvaluatorBase::OutArgsSetup<Scalar> outArgs;
    outArgs.setModelEvalDescription(this->description());
    outArgs.setSupports( Thyra::ModelEvaluatorBase::OUT_ARG_f );
    outArgs_ = outArgs;
  }

  // Set up nominal values.
  nominalValues_ = inArgs_;
  const Teuchos::RCP<Thyra::VectorBase<Scalar>> x_ic = createMember(x_space_);
  { // scope to delete DetachedVectorView
    Thyra::DetachedVectorView<Scalar> x_ic_view( *x_ic );
    x_ic_view[0] = 1.0;
    x_ic_view[1] = 1.0;
  }
  nominalValues_.set_x(x_ic);

  isInitialized_ = true;
}
