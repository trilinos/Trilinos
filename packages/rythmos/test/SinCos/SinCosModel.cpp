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

#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_DetachedMultiVectorView.hpp"
#include "Thyra_DefaultSerialDenseLinearOpWithSolveFactory.hpp"
#include "Thyra_DefaultMultiVectorLinearOpWithSolve.hpp"
#include "Thyra_DefaultLinearOpSource.hpp"

#ifdef SINCOSMODEL_DEBUG
#include <iostream>
#endif // SINCOSMODEL_DEBUG

namespace Rythmos {

// non-member Constructor
RCP<SinCosModel> sinCosModel(bool implicit) 
{
  RCP<SinCosModel> model = rcp(new SinCosModel);
  model->setImplicitFlag(implicit);
  return(model);
}

// Constructor
SinCosModel::SinCosModel()
{
  dim_ = 2;
  isImplicit_ = false;
  isInitialized_ = false;

  // Create x_space and f_space
  x_space_ = Thyra::defaultSpmdVectorSpace<double>(dim_);
  f_space_ = Thyra::defaultSpmdVectorSpace<double>(dim_);
}

void SinCosModel::setImplicitFlag(bool implicit) 
{
  if (isImplicit_ != implicit) {
    isInitialized_ = false;
  }
  isImplicit_ = implicit;
  initialize_();
}


ModelEvaluatorBase::InArgs<double> SinCosModel::getExactSolution(double t) const
{
  TEST_FOR_EXCEPTION( !isInitialized_, std::logic_error,
      "Error, setImplicitFlag must be called first!\n"
      );
  ModelEvaluatorBase::InArgs<double> inArgs = inArgs_;
  double exact_t = t;
  inArgs.set_t(exact_t);
  RCP<VectorBase<double> > exact_x = createMember(x_space_);
  { // scope to delete DetachedVectorView
    Thyra::DetachedVectorView<double> exact_x_view(*exact_x);
    exact_x_view[0] = sin(t);
    exact_x_view[1] = cos(t);
  }
  inArgs.set_x(exact_x);
  if (isImplicit_) {
    RCP<VectorBase<double> > exact_x_dot = createMember(x_space_);
    { // scope to delete DetachedVectorView
      Thyra::DetachedVectorView<double> exact_x_dot_view(*exact_x_dot);
      exact_x_dot_view[0] = cos(t);
      exact_x_dot_view[1] = -sin(t);
    }
    inArgs.set_x_dot(exact_x_dot);
  }
  return(inArgs);
}

RCP<const Thyra::VectorSpaceBase<double> >
SinCosModel::get_x_space() const
{
  return x_space_;
}


RCP<const Thyra::VectorSpaceBase<double> >
SinCosModel::get_f_space() const
{
  return f_space_;
}


ModelEvaluatorBase::InArgs<double>
SinCosModel::getNominalValues() const
{
  TEST_FOR_EXCEPTION( !isInitialized_, std::logic_error,
      "Error, setImplicitFlag must be called first!\n"
      );
  return nominalValues_;
}




RCP<Thyra::LinearOpWithSolveBase<double> >
SinCosModel::create_W() const
{
  RCP<const Thyra::LinearOpWithSolveFactoryBase<double> > W_factory = this->get_W_factory();
  RCP<Thyra::LinearOpBase<double> > matrix = this->create_W_op();
  RCP<Thyra::LinearOpWithSolveBase<double> > W = 
    Thyra::linearOpWithSolve<double>(
      *W_factory,
      matrix
      );
  return W;
}
//RCP<Thyra::LinearOpWithSolveBase<double> >
//SinCosModel::create_W() const
//{
//  return Thyra::multiVectorLinearOpWithSolve<double>();
//}





RCP<Thyra::LinearOpBase<double> >
SinCosModel::create_W_op() const
{
  RCP<Thyra::MultiVectorBase<double> > matrix = Thyra::createMembers(x_space_, dim_);
  return(matrix);
}




RCP<const Thyra::LinearOpWithSolveFactoryBase<double> > 
SinCosModel::get_W_factory() const
{
  RCP<Thyra::LinearOpWithSolveFactoryBase<double> > W_factory = 
    Thyra::defaultSerialDenseLinearOpWithSolveFactory<double>();
  return W_factory;
}
//  RCP<const Thyra::LinearOpBase<double> > fwdOp = this->create_W_op();
//  RCP<Thyra::LinearOpWithSolveBase<double> > W = 
//    Thyra::linearOpWithSolve<double>(
//      *W_factory,
//      fwdOp
//      );
//  W_factory->initializeOp(
//      Thyra::defaultLinearOpSource<double>(fwdOp),
//      &*W,
//      Thyra::SUPPORT_SOLVE_UNSPECIFIED
//      );



ModelEvaluatorBase::InArgs<double>
SinCosModel::createInArgs() const
{
  TEST_FOR_EXCEPTION( !isInitialized_, std::logic_error,
      "Error, setImplicitFlag must be called first!\n"
      );
  return inArgs_;
}


// Private functions overridden from ModelEvaulatorDefaultBase


ModelEvaluatorBase::OutArgs<double>
SinCosModel::createOutArgsImpl() const
{
  TEST_FOR_EXCEPTION( !isInitialized_, std::logic_error,
      "Error, setImplicitFlag must be called first!\n"
      );
  return outArgs_;
}


void SinCosModel::evalModelImpl(
  const ModelEvaluatorBase::InArgs<double> &inArgs,
  const ModelEvaluatorBase::OutArgs<double> &outArgs
  ) const
{
  TEST_FOR_EXCEPTION( !isInitialized_, std::logic_error,
      "Error, setImplicitFlag must be called first!\n"
      );

  const RCP<const VectorBase<double> > x_in = inArgs.get_x().assert_not_null();
  Thyra::ConstDetachedVectorView<double> x_in_view( *x_in ); 

  //double t = inArgs.get_t();

  RCP<const VectorBase<double> > x_dot_in;
  double beta = inArgs.get_beta();
  double alpha = -1.0;
  if (isImplicit_) {
    x_dot_in = inArgs.get_x_dot().assert_not_null();
    alpha = inArgs.get_alpha();
  }

  const RCP<VectorBase<double> > f_out = outArgs.get_f();
  const RCP<Thyra::LinearOpBase<double> > W_out = outArgs.get_W_op();


  if (!isImplicit_) { // isImplicit_ == false
    if (!is_null(f_out)) {
      Thyra::DetachedVectorView<double> f_out_view( *f_out ); 
      f_out_view[0] = x_in_view[1];
      f_out_view[1] = -x_in_view[0];
    }
    if (!is_null(W_out)) {
      RCP<Thyra::MultiVectorBase<double> > matrix = Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<double> >(W_out,true);
      Thyra::DetachedMultiVectorView<double> matrix_view( *matrix );
      matrix_view(0,0) = 0.0;
      matrix_view(0,1) = +beta;
      matrix_view(1,0) = -beta;
      matrix_view(1,1) = 0.0;
    }
  } else { // isImplicit_ == true
    if (!is_null(f_out)) {
      Thyra::DetachedVectorView<double> f_out_view( *f_out ); 
      Thyra::ConstDetachedVectorView<double> x_dot_in_view( *x_dot_in );
      f_out_view[0] = x_dot_in_view[0] - x_in_view[1];
      f_out_view[1] = x_dot_in_view[1] + x_in_view[0];
    }
    if (!is_null(W_out)) {
      RCP<Thyra::MultiVectorBase<double> > matrix = Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<double> >(W_out,true);
      Thyra::DetachedMultiVectorView<double> matrix_view( *matrix );
      matrix_view(0,0) = alpha;
      matrix_view(0,1) = -beta;
      matrix_view(1,0) = +beta;
      matrix_view(1,1) = alpha;
    }
  }
}

// private

void SinCosModel::initialize_() 
{
  if (!isInitialized_) {
    
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
    { // scope to delete DetachedVectorView
      Thyra::DetachedVectorView<double> x_ic_view( *x_ic );
      x_ic_view[0] = 0.0;
      x_ic_view[1] = 1.0;
    }
    nominalValues_.set_x(x_ic);
    double t_ic = 0.0;
    nominalValues_.set_t(t_ic);
    if (isImplicit_) {
      const RCP<VectorBase<double> > x_dot_ic = createMember(x_space_);
      { // scope to delete DetachedVectorView
        Thyra::DetachedVectorView<double> x_dot_ic_view( *x_dot_ic );
        x_dot_ic_view[0] = 1.0;
        x_dot_ic_view[1] = 0.0;
      }
      nominalValues_.set_x_dot(x_dot_ic);
    }
    isInitialized_ = true;

  }

}

} // namespace Rythmos

