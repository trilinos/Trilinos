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

#ifndef Rythmos_STEPPER_H
#define Rythmos_STEPPER_H


#include "Teuchos_RefCountPtr.hpp"
#include "Thyra_VectorBase.hpp"
#include "NonlinearModel.hpp"

namespace Rythmos {

//-----------------------------------------------------------------------------
// Class         : Stepper
// Purpose       : Base class for defining stepper functionality
// Special Notes :
// Creator       : Todd Coffey, SNL
// Creation Date : 05/26/05
//-----------------------------------------------------------------------------
template<class Scalar> 
class Stepper
{
  public:
    
    // Destructor
    ~Stepper();

    // Cosntructor
    // 05/26/05 tscoffe:  I think I'll pass in a parameter list at construction
    // time also to specify what stepper options I want.
    Stepper(Teuchos::RefCountPtr<NonlinearModel<Scalar> > &model)
      { model_ = model; };

    // Take a step _no larger_ than dt 
    virtual Scalar TakeStep(Scalar dt)=0;
   
    // Take a step 
    virtual Scalar TakeStep()=0;

    // Get solution vector
    Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > &get_solution()
      { return(solution_vector_); };
    
    // Get residual vector
    Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > &get_residual()
    { return(residual_vector_); };

  protected:

    Teuchos::RefCountPtr<NonlinearModel<Scalar> > model_;
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > solution_vector_;
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > residual_vector_;
    Scalar t_;

};

} // namespace Rythmos

#endif //Rythmos_STEPPER_H
