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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER


#ifndef RYTHMOS_INTERPOLATOR_ACCEPTING_OBJECT_BASE_HPP
#define RYTHMOS_INTERPOLATOR_ACCEPTING_OBJECT_BASE_HPP


#include "Rythmos_Types.hpp"
#include "Rythmos_InterpolatorBase.hpp"

namespace Rythmos {

/** \brief Mix-in interface for objects that accept an interpolator object.
 *
 * ToDo: Finish documentation!
 */

template<class Scalar>
class InterpolatorAcceptingObjectBase 
{
public:

  /** \brief . */
  virtual ~InterpolatorAcceptingObjectBase() {}

  /** \brief . */
  virtual void setInterpolator(
    const RCP<InterpolatorBase<Scalar> > &interpolator
    ) = 0;

  /** \brief . */
  virtual RCP<InterpolatorBase<Scalar> >
    getNonconstInterpolator() = 0;

  /** \brief . */
  virtual RCP<const InterpolatorBase<Scalar> >
    getInterpolator() const = 0;

  /** \brief . */
  virtual RCP<InterpolatorBase<Scalar> >
    unSetInterpolator() = 0;

};


} // namespace Rythmos


#endif // RYTHMOS_INTERPOLATOR_ACCEPTING_OBJECT_BASE_HPP
