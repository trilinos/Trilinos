//@HEADER

// ***********************************************************************
//
//                     Rythmos Package
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

#include "SinCosModel.hpp"

#ifdef SINCOSMODEL_DEBUG
#include <iostream>
#endif

namespace Rythmos {

// non-member Constructor
RCP<SinCosModel> sinCosModel(bool implicit) 
{
  RCP<SinCosModel> model = rcp(new SinCosModel);
  if (implicit) {
    model->setImplicitFlag(implicit);
  }
  return(model);
}

// Constructor
SinCosModel::SinCosModel()
{
  dim_ = 2;
  isImplicit_ = false;
  isInitialized_ = false;
}

void SinCosModel::setImplicitFlag(bool implicit) 
{
  if (isImplicit_ != implicit) {
    isImplicit_ = implicit;
    isInitialized_ = false;
  }
}

ModelEvaluatorBase::InArgs<double> SinCosModel::getExactSolution(double t) const
{
  initialize_();
  ModelEvaluatorBase::InArgs<double> inArgs = inArgs_;
  double exact_t = t;
  inArgs.set_t(exact_t);
  RCP<VectorBase<double> > exact_x = createMember(x_space_);
  // get access to underlying double* for exact_x
  // Fill as follows:
  // exact_x[0] = sin(t)
  // exact_x[1] = cos(t)
  inArgs.set_x(exact_x);
  if (isImplicit_) {
    RCP<VectorBase<double> > exact_x_dot = createMember(x_space_);
    // get access to underlying double* for exact_x_dot
    // Fill as follows:
    // exact_x_dot[0] = cos(t)
    // exact_x_dot[1] = -sin(t)
    inArgs.set_x_dot(exact_x_dot);
  }
  return(inArgs);
}

RCP<const VectorSpaceBase<double> >
SinCosModel::get_x_space() const
{
  initialize_();
  return x_space_;
}


RCP<const VectorSpaceBase<double> >
SinCosModel::get_f_space() const
{
  initialize_();
  return f_space_;
}


ModelEvaluatorBase::InArgs<double>
SinCosModel::getNominalValues() const
{
  initialize_();
  return nominalValues_;
}


RCP<Thyra::LinearOpWithSolveBase<double> >
SinCosModel::create_W() const
{
  initialize_();

  RCP<Thyra::MultiVectorBase<double> > matrix = Thyra::createMembers(x_space_, dim_);
  RCP<Thyra::LinearOpWithSolveBase<double> > W = 
    Thyra::linearOpWithSolve<double>(
      *Thyra::defaultSerialDenseLinearOpWithSolveFactory<double>(),
      matrix
      );
  return W;
}


RCP<Thyra::LinearOpBase<double> >
SinCosModel::create_W_op() const
{
  initialize_();
  return(create_W());
}


ModelEvaluatorBase::InArgs<double>
SinCosModel::createInArgs() const
{
  initialize_();
  return inArgs_;
}


// Private functions overridden from ModelEvaulatorDefaultBase


ModelEvaluatorBase::OutArgs<double>
SinCosModel::createOutArgsImpl() const
{
  initialize_();
  return outArgs_;
}


void SinCosModel::evalModelImpl(
  const ModelEvaluatorBase::InArgs<double> &inArgs,
  const ModelEvaluatorBase::OutArgs<double> &outArgs
  ) const
{
  const RCP<const VectorBase<double> > x_in = inArgs.get_x().assert_not_null();
  double t = inArgs.get_t();

  const RCP<const VectorBase<double> > x_dot_in;
  double beta = inArgs.get_beta();
  double alpha;
  if (isImplicit_) {
    x_dot_in = inArgs.get_x_dot().assert_not_null();
    alpha = inArgs.get_alpha();
  }

  const RCP<VectorBase<double> > f_out = outArgs.get_f();
  // TODO:  Should this be get_W?
  const RCP<Thyra::LinearOpBase<double> > W_op_out = outArgs.get_W_op();

  // Get direct access to underlying double* for x_in
  // Get direct access to underlying double* for f_out
  if (!isImplicit_) { // isImplicit_ == false
    if (!is_null(f_out)) {
      // Load the following into f_out:
      // f_out[0] = x_in[1]
      // f_out[1] = -x_in[0]
    }
    if (!is_null(W_op_out)) {
      // Load the following into W_op_out:
      // W_op_out[0][0] = 0.0
      // W_op_out[0][1] = +beta
      // W_op_out[1][0] = -beta
      // W_op_out[1][1] = 0.0
    }
  } else { // isImplicit_ == true
    if (!is_null(f_out)) {
      // Get direct access to underlying double* for x_dot_in
      // Get direct access to underlying double** for W_op_out
      // Load the following into f_out:
      // f_out[0] = x_dot_in[0] - x_in[1]
      // f_out[1] = x_dot_in[1] + x_int[0]
    }
    if (!is_null(W_op_out)) {
      // Load the following into W_op_out:
      // W_op_out[0][0] = alpha
      // W_op_out[0][1] = -beta
      // W_op_out[1][0] = +beta
      // W_op_out[1][1] = alpha
    }
  }

}

// private

void SinCosModel::initialize_() const
{
  if (!isInitialized_) {
    
    // Create x_space and f_space
    x_space_ = Thyra::defaultSpmdVectorSpaceBase<double>(dim_);
    f_space_ = Thyra::defaultSpmdVectorSpaceBase<double>(dim_);

    {
      // Set up prototypical InArgs
      ModelEvaluatorBase::InArgsSetup<double> inArgs;
      inArgs.setModelEvalDescription(this->description());
      inArgs.setSupports( ModelEvaluatorBase::IN_ARG_t );
      inArgs.setSupports( ModelEvaluatorBase::IN_ARG_x );
      inArgs.setSupports( ModelEvaluatorBase::IN_ARG_beta );
      if (isImplicit_) {
        inArgs.setSupports( ModelEvaluatorBase::IN_ARG_x_dot );
        inArgs.setSupports( ModelEvaluatorBase::IN_ARG_alpha );
      }
      inArgs_ = inArgs;
    }
    {
      // Set up prototypical OutArgs
      ModelEvaluatorBase::OutArgsSetup<double> outArgs;
      outArgs.setModelEvalDescription(this->description());
      outArgs.setSupports( ModelEvaluatorBase::OUT_ARG_f );
      outArgs.setSupports( ModelEvaluatorBase::OUT_ARG_W_op );
      outArgs_ = outArgs;
    }

    // Set up nominal values 
    nominalValues_ = inArgs_;
    const RCP<VectorBase<double> > x_ic = createMember(x_space_);
    // access underlying double *
    // Fill with x_ic[0] = 0.0
    //           x_ic[1] = 1.0
    nominalValues_.set_x(x_ic);
    double t_ic = 0.0;
    nominalValues_.set_t(t_ic);
    if (isImplicit_) {
      const RCP<VectorBase<double> > x_dot_ic = createMember(x_space_);
      // access underlying double *
      // Fill with x_dot_ic[0] = 1.0
      //           x_dot_ic[1] = 0.0
      nominalValues_.set_x_dot(x_dot_ic);
    }
    isInitialized_ = true;

  }

}

} // namespace Rythmos

