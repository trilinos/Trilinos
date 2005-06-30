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

#include "Thyra_EpetraModelEvaluator.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "EpetraExt_ModelEvaluator.hpp"

namespace Thyra {

// Constructors/initializers/accessors.

EpetraModelEvaluator::EpetraModelEvaluator()
{}

EpetraModelEvaluator::EpetraModelEvaluator( const Teuchos::RefCountPtr<const EpetraExt::ModelEvaluator>& epetraModel )
{
  initialize(epetraModel);
}

void EpetraModelEvaluator::initialize( const Teuchos::RefCountPtr<const EpetraExt::ModelEvaluator>& epetraModel )
{
  epetraModel_ = epetraModel;
  space_x_ = create_MPIVectorSpaceBase(epetraModel_->get_x_map());
  space_f_ = create_MPIVectorSpaceBase(epetraModel_->get_f_map());
}

Teuchos::RefCountPtr<const EpetraExt::ModelEvaluator> EpetraModelEvaluator::getEpetraModel() const
{
  return epetraModel_;
}

void EpetraModelEvaluator::uninitialize( Teuchos::RefCountPtr<const EpetraExt::ModelEvaluator>* epetraModel )
{
  if(epetraModel) *epetraModel = epetraModel_;
  epetraModel_ = Teuchos::null;
}

// Overridden from ModelEvaulator.

Teuchos::RefCountPtr<const VectorSpaceBase<double> >
EpetraModelEvaluator::get_x_space() const
{
  return space_x_;
}

Teuchos::RefCountPtr<const VectorSpaceBase<double> >
EpetraModelEvaluator::get_f_space() const
{
  return space_f_;
}

Teuchos::RefCountPtr<const VectorBase<double> >
EpetraModelEvaluator::get_x_init() const
{
  return create_MPIVectorBase( epetraModel_->get_x_init(), space_x_ );
}

double EpetraModelEvaluator::get_t_init() const
{
  return epetraModel_->get_t_init();
}

EpetraModelEvaluator::InArgs<double> EpetraModelEvaluator::createInArgs() const
{
  InArgsSetup<double> inArgs;
  typedef EpetraExt::ModelEvaluator EME;
  EME::InArgs
    epetraInArgs = epetraModel_->createInArgs();
  inArgs.setSupports(IN_ARG_x,epetraInArgs.supports(EME::IN_ARG_x));
  inArgs.setSupports(IN_ARG_t,epetraInArgs.supports(EME::IN_ARG_t));
  return inArgs;
}

EpetraModelEvaluator::OutArgs<double> EpetraModelEvaluator::createOutArgs() const
{
  OutArgsSetup<double> outArgs;
  typedef EpetraExt::ModelEvaluator EME;
  EME::OutArgs
    epetraOutArgs = epetraModel_->createOutArgs();
  outArgs.setSupports(OUT_ARG_f,epetraOutArgs.supports(EME::OUT_ARG_f));
  return outArgs;
}

void EpetraModelEvaluator::evalModel( const InArgs<double>& inArgs, const OutArgs<double>& outArgs ) const
{
  // ToDo: Fill in!
  TEST_FOR_EXCEPT(true);
}

} // namespace Thyra
