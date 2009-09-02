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

#include "PolynomialModel.hpp"

#include "Teuchos_Polynomial.hpp"

#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_DetachedMultiVectorView.hpp"

#ifdef SINCOSMODEL_DEBUG
#include <iostream>
#endif

namespace Rythmos {

// non-member Constructor
RCP<PolynomialModel> polynomialModel(const RCP<const Teuchos::Polynomial<double> >& poly) 
{
  RCP<PolynomialModel> model = rcp(new PolynomialModel);
  model->setPolynomial(poly);
  return(model);
}

// Constructor
PolynomialModel::PolynomialModel()
{
  dim_ = 1;
  isInitialized_ = false;

  // Create x_space and f_space
  x_space_ = Thyra::defaultSpmdVectorSpace<double>(dim_);
  f_space_ = Thyra::defaultSpmdVectorSpace<double>(dim_);
}

void PolynomialModel::setPolynomial( const RCP<const Teuchos::Polynomial<double> >& poly )
{
  poly_ = poly;
  initialize_();
}


ModelEvaluatorBase::InArgs<double> PolynomialModel::getExactSolution(double t) const
{
  TEST_FOR_EXCEPTION( !isInitialized_, std::logic_error,
      "Error, setPolynomial must be called first!\n"
      );
  ModelEvaluatorBase::InArgs<double> inArgs = inArgs_;
  double exact_t = t;
  inArgs.set_t(exact_t);
  RCP<VectorBase<double> > exact_x = createMember(x_space_);
  { // scope to delete DetachedVectorView
    Thyra::DetachedVectorView<double> exact_x_view(*exact_x);
    exact_x_view[0] = 0.0; // TODO
  }
  inArgs.set_x(exact_x);
  return(inArgs);
}

RCP<const Thyra::VectorSpaceBase<double> >
PolynomialModel::get_x_space() const
{
  return x_space_;
}


RCP<const Thyra::VectorSpaceBase<double> >
PolynomialModel::get_f_space() const
{
  return f_space_;
}


ModelEvaluatorBase::InArgs<double>
PolynomialModel::getNominalValues() const
{
  TEST_FOR_EXCEPTION( !isInitialized_, std::logic_error,
      "Error, setPolynomial must be called first!\n"
      );
  return nominalValues_;
}



ModelEvaluatorBase::InArgs<double>
PolynomialModel::createInArgs() const
{
  TEST_FOR_EXCEPTION( !isInitialized_, std::logic_error,
      "Error, setPolynomial must be called first!\n"
      );
  return inArgs_;
}


// Private functions overridden from ModelEvaulatorDefaultBase


ModelEvaluatorBase::OutArgs<double>
PolynomialModel::createOutArgsImpl() const
{
  TEST_FOR_EXCEPTION( !isInitialized_, std::logic_error,
      "Error, setPolynomial must be called first!\n"
      );
  return outArgs_;
}


void PolynomialModel::evalModelImpl(
  const ModelEvaluatorBase::InArgs<double> &inArgs,
  const ModelEvaluatorBase::OutArgs<double> &outArgs
  ) const
{
  TEST_FOR_EXCEPTION( !isInitialized_, std::logic_error,
      "Error, setPolynomial must be called first!\n"
      );

  //const RCP<const VectorBase<double> > x_in = inArgs.get_x().assert_not_null();
  //Thyra::ConstDetachedVectorView<double> x_in_view( *x_in ); 

  double t = inArgs.get_t();

  double p;
  poly_->evaluate(t,&p);

  const RCP<VectorBase<double> > f_out = outArgs.get_f();

  if (!is_null(f_out)) {
      Thyra::DetachedVectorView<double> f_out_view( *f_out ); 
      f_out_view[0] = p;
  }
}

// private

void PolynomialModel::initialize_() 
{
  if (!isInitialized_) {
    
    {
      // Set up prototypical InArgs
      ModelEvaluatorBase::InArgsSetup<double> inArgs;
      inArgs.setModelEvalDescription(this->description());
      inArgs.setSupports( ModelEvaluatorBase::IN_ARG_t );
      inArgs_ = inArgs;
    }
    {
      // Set up prototypical OutArgs
      ModelEvaluatorBase::OutArgsSetup<double> outArgs;
      outArgs.setModelEvalDescription(this->description());
      outArgs.setSupports( ModelEvaluatorBase::OUT_ARG_f );
      outArgs_ = outArgs;
    }

    // Set up nominal values 
    nominalValues_ = inArgs_;
    double t_ic = 0.0;
    nominalValues_.set_t(t_ic);
    isInitialized_ = true;

  }

}

} // namespace Rythmos

