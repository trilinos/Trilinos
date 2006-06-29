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

#ifndef Rythmos_INTERPOLATION_BUFFER_H
#define Rythmos_INTERPOLATION_BUFFER_H

#include "Thyra_VectorBase.hpp"
#include "Teuchos_Describable.hpp"

namespace Rythmos {

/** \brief Base class for defining interpolation buffer functionality. */
template<class Scalar> 
class InterpolationBuffer : virtual public Teuchos::Describable
{
  public:
    
    /// Destructor
    virtual ~InterpolationBuffer() {};

    /// Add point to buffer
    virtual bool SetPoint(const ScalarMag& time, const Thyra::VectorBase<Scalar>& x, const Thyra::VectorBase<Scalar>& xdot);
    // 06/29/06 tscoffe:  This could be limiting, especially for extra
    // parameters passed through InArgs.  Also, in the ODE case, xdot = f(x),
    // so there's the question of including additional functions from OutArgs
    // in the interpolation buffer also.  One fix would be to use
    // std::vector<Thyra::VectorBase<Scalar> > for (x,xdot).  Would this handle
    // parameters in InArgs correctly?   If we use std::vector though, we don't
    // know that the first two elements represent an x, xdot relationship which
    // is advantagous for many interpolation algorithms.  An alternative would
    // be to use (x, xdot, std::vector<Thyra::VectorBase<Scalar> > u) as input,
    // with u optional.  And this adds further complications for SetRange when
    // the other interpolation buffer uses optional inputs and this one doesn't
    // or vice versa.  Another strategy would use InArgs directly as input and
    // process interpolation over all the valid data present.  The
    // InArgs/OutArgs approach implies that InArgs are the independent
    // variables and OutArgs are the dependent variables.  In the basic case of
    // storing the solution to a DE, the independent variable is time and the
    // dependent variable is everything else in the InArgs object.  This
    // splitting of the object is awkward.  Fundamentally, our independent
    // variable is a ScalarMag labeled "time" and everything else is part of
    // the dependent variable that will be interpolated between values of
    // "time" provided.

    /// Get value from buffer
    virtual bool GetPoint(const ScalarMag& time, Thyra::VectorBase<Scalar>* x, Thyra::VectorBase<Scalar>* xdot) const;

    /// Fill data in from another interpolation buffer
    virtual bool SetRange(const ScalarMagRange& time_range, const InterpolationBuffer<Scalar>& IB);

    /// Get interpolation nodes
    virtual bool GetNodes(std::vector<ScalarMag>* time_list) const;
};

} // namespace Rythmos

#endif //Rythmos_INTERPOLATION_BUFFER_H
