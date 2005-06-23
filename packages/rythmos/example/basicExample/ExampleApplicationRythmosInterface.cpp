//
// @HEADER
// ***********************************************************************
// 
//                           Rythmos Package
//                 Copyright (2005) Sandia Corporation
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

#include "ExampleApplicationRythmosInterface.hpp"

ExampleApplicationRythmosInterface::ExampleApplicationRythmosInterface(Teuchos::ParameterList &params)
{
  problem_ = Teuchos::rcp(new ExampleApplication(params));
  epetra_map_ = problem_->get_epetra_map();
//  std::cout << "ExampleApplicationRythmosInterface::ExampleApplicationRythmosInterface()" << std::endl;
//  std::cout << "epetra_map_.get() = " << std::endl;
//  std::cout << epetra_map_.get() << std::endl;
  thyra_vs_ = Thyra::create_MPIVectorSpaceBase(epetra_map_);
}

ExampleApplicationRythmosInterface::ExampleApplicationRythmosInterface()
{
}

ExampleApplicationRythmosInterface::~ExampleApplicationRythmosInterface()
{
}

int ExampleApplicationRythmosInterface::evalModel(const Rythmos::InArgs<double> &inargs, const Rythmos::OutArgs<double> &outargs) const
{
  // input arguments:
  Teuchos::RefCountPtr<Thyra::VectorBase<double> > x = inargs.get_x();
  double t = inargs.get_t();

  // output arguments:
  Teuchos::RefCountPtr<Thyra::VectorBase<double> > F = outargs.get_F();

  (*problem_).evalResidual(&*(Thyra::get_Epetra_Vector(*epetra_map_,F)),*(Thyra::get_Epetra_Vector(*epetra_map_,x)),t);
  return 0;
}

Teuchos::RefCountPtr<Thyra::VectorBase<double> > ExampleApplicationRythmosInterface::get_vector() const
{
  return(Thyra::create_MPIVectorBase(problem_->get_x0(),thyra_vs_));
}

Teuchos::RefCountPtr<const Epetra_Map> ExampleApplicationRythmosInterface::get_Epetra_Map() const
{ 
  return(epetra_map_); 
}
