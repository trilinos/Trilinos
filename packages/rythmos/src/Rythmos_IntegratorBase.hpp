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

#ifndef Rythmos_INTEGRATOR_BASE_H
#define Rythmos_INTEGRATOR_BASE_H

#include "Rythmos_InterpolationBufferBase.hpp"
#include "Rythmos_StepperBase.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_as.hpp"


namespace Rythmos {


/** \brief Abstract interface for time integrators.
 *
 * A time integrator accepts a fully initialized stepper object and then
 * carries out the time integrator in some fasion.  The client drives the
 * integrator by requesting value of the state at different points in time.
 */
template<class Scalar> 
class IntegratorBase
  : virtual public Rythmos::InterpolationBufferBase<Scalar>
{
public:
  
  /** \brief . */
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;

  /** \name Constructors, Initializers, Misc */
  //@{
  
  /** \brief . */
  virtual void setInterpolationBuffer(
    const RCP<Rythmos::InterpolationBufferBase<Scalar> > &trailingInterpBuffer
    ) =0;

  /** \brief . */
  virtual void setStepper(
    const RCP<Rythmos::StepperBase<Scalar> > &stepper,
    const Scalar &finalTime
    ) =0;

  /** \brief This is a non-const version of getPoints which allows the
   * integrator class to step forward to get the points asked for.
   *
   * The time values are assumed to be sorted on input.
   * */
  virtual bool getFwdPoints(
    const Array<Scalar>& time_vec,
    Array<RCP<const Thyra::VectorBase<Scalar> > >* x_vec,
    Array<RCP<const Thyra::VectorBase<Scalar> > >* xdot_vec,
    Array<ScalarMag>* accuracy_vec
    ) =0;

  /** \brief 
   * This function returns the valid range of points that the integrator can integrate over.
   * This is different from getTimeRange because that only returns the time
   * range over which the integrator already has data.
   * */
  virtual TimeRange<Scalar> getFwdTimeRange() const =0;

  //@}

};
 
}

#endif // Rythmos_INTEGRATOR_BASE_H

