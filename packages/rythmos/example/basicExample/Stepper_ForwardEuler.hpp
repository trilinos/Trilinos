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

#ifndef Rythmos_STEPPER_FORWARDEULER_H
#define Rythmos_STEPPER_FORWARDEULER_H

#include "Stepper.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Thyra_VectorBase.hpp"
#include "NonlinearModel.hpp"

namespace Rythmos {

//-----------------------------------------------------------------------------
// Class         : FowardEuler
// Purpose       : Base class for defining stepper functionality
// Special Notes :
// Creator       : Todd Coffey, SNL
// Creation Date : 05/26/05
//-----------------------------------------------------------------------------
template<class Scalar>
class ForwardEuler : public virtual Stepper<Scalar>
{
  public:
    
    // Constructor
    ForwardEuler();
    ForwardEuler(const Teuchos::RefCountPtr<NonlinearModel<Scalar> > &model_);
    
    // Destructor
    ~ForwardEuler();

    // Take a step _no larger_ than dt 
    Scalar TakeStep(Scalar dt);
   
    // Take a step 
    Scalar TakeStep();

    // Get solution vector
    const Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > &get_solution();

    // Get residual vector
    const Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > &get_residual();

  protected:

    Teuchos::RefCountPtr<NonlinearModel<Scalar> > model_;
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > solution_vector_;
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > residual_vector_;
    Scalar t_;

};

} // namespace Rythmos

#endif //Rythmos_STEPPER_FORWARDEULER_H
