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

#ifndef Rythmos_STEPPER_H
#define Rythmos_STEPPER_H


#include "Teuchos_RefCountPtr.hpp"
#include "Thyra_VectorBase.hpp"
#include "Teuchos_Describable.hpp"
#include "Rythmos_InterpolationBufferBase.hpp"

namespace Rythmos {

/** \brief Base class for defining stepper functionality. */
template<class Scalar> 
class Stepper : virtual public InterpolationBufferBase<Scalar>
{
  public:
    
    /// Destructor
    virtual ~Stepper() {};

    /// Take a step _no larger_ than dt 
    virtual Scalar TakeStep(Scalar dt)=0;
   
    /// Take a step 
    virtual Scalar TakeStep()=0;

    /// Get solution vector
    virtual Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > get_solution() const = 0;

};

} // namespace Rythmos

#endif //Rythmos_STEPPER_H
