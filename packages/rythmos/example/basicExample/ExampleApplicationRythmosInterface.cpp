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

//-----------------------------------------------------------------------------
// Function      : ExampleApplicationRythmosInterface::ExampleApplicationRythmosInterface
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 05/17/05
//-----------------------------------------------------------------------------
ExampleApplicationRythmosInterface::ExampleApplicationRythmosInterface()
{
  // 05/26/05 tscoffe:  This is where a parameter list could be passed in and
  // used in constructing the problem.
  double lambda = -0.5;
  int numelements = 1;
  problem_ = Teuchos::rcp(new ExampleApplication(lambda,numelements));
  epetra_map_ = problem_->get_epetra_map();
  thyra_vs_ = Thyra::create_MPIVectorSpaceBase(epetra_map_);
}

//-----------------------------------------------------------------------------
// Function      : ExampleApplication::~ExampleApplication
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 05/05/05
//-----------------------------------------------------------------------------
ExampleApplicationRythmosInterface::~ExampleApplicationRythmosInterface()
{
}

//-----------------------------------------------------------------------------
// Function      : ExampleApplication::evalModel
// Purpose       : Evaluate residual
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 05/17/05
//-----------------------------------------------------------------------------
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

//-----------------------------------------------------------------------------
// Function      : ExampleApplication::get_vector
// Purpose       : Get nominal vector
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 05/26/05
//-----------------------------------------------------------------------------
Teuchos::RefCountPtr<Thyra::VectorBase<double> > ExampleApplicationRythmosInterface::get_vector() const
{
  return(Thyra::create_MPIVectorBase(problem_->get_x0(),thyra_vs_));
}

//-----------------------------------------------------------------------------
// Function      : ExampleApplication::get_Epetra_Map
// Purpose       : Get Epetra Map
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, SNL
// Creation Date : 06/07/05
//-----------------------------------------------------------------------------
Teuchos::RefCountPtr<const Epetra_Map> ExampleApplicationRythmosInterface::get_Epetra_Map() const
{ 
  return(epetra_map_); 
}

