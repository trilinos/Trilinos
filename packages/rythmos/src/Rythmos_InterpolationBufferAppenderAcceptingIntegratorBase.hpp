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


#ifndef RYTHMOS_INTERPOLATION_BUFFER_APPENDER_ACCEPTING_INTEGRATOR_BASE_HPP
#define RYTHMOS_INTERPOLATION_BUFFER_APPENDER_ACCEPTING_INTEGRATOR_BASE_HPP


#include "Rythmos_Types.hpp"
#include "Rythmos_IntegratorBase.hpp"
#include "Rythmos_InterpolationBufferAppenderBase.hpp"

namespace Rythmos {

/** \brief Mix-in interface for integrator objects that accept an
 * interpolationBufferAppender object to be used for appending to the trailing
 * interpolation buffer
 *
 * ToDo: Finish documentation!
 */

template<class Scalar>
class InterpolationBufferAppenderAcceptingIntegratorBase : virtual public IntegratorBase<Scalar>
{
public:

  /** \brief . */
  virtual void setInterpolationBufferAppender(
    const RCP<InterpolationBufferAppenderBase<Scalar> > &interpBufferAppender
    ) = 0;

  /** \brief . */
  virtual RCP<const InterpolationBufferAppenderBase<Scalar> >
    getInterpolationBufferAppender() = 0;

  /** \brief . */
  virtual RCP<InterpolationBufferAppenderBase<Scalar> >
    getNonconstInterpolationBufferAppender() = 0;

  /** \brief . */
  virtual RCP<InterpolationBufferAppenderBase<Scalar> >
    unSetInterpolationBufferAppender() = 0;

};


} // namespace Rythmos


#endif // RYTHMOS_INTERPOLATION_BUFFER_APPENDER_ACCEPTING_INTEGRATOR_BASE_HPP
