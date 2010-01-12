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
#include "Teuchos_as.hpp"


namespace Rythmos {


namespace Exceptions {


/** brief Thrown if <tt>Rythmos::IntegratorBase::getFwdPoints()</tt> fails to
 * return given time points.
 */
class GetFwdPointsFailed : public ::Rythmos::Exceptions::ExceptionBase
{public: GetFwdPointsFailed(const std::string &my_what):ExceptionBase(my_what) {}};


} // namespace Exceptions


/** \brief Abstract interface for time integrators.
 *
 * A time integrator accepts a fully initialized stepper object (and a final
 * time) and then carries out the time integration in some fasion.  The client
 * drives the integrator by requesting value of the state at different points
 * in time.  If possible, the client should request time points only forward
 * in time if possible since.
 */
template<class Scalar> 
class IntegratorBase : virtual public InterpolationBufferBase<Scalar>
{
public:
  
  /** \brief . */
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;

  /** \brief . */
  virtual RCP<IntegratorBase<Scalar> > cloneIntegrator() const
    { return Teuchos::null; }

  /** \brief Specify the stepper to use for integration which effectively
   * reinitializes the intergrator.
   *
   * \param stepper [inout,persisting] Gives the stepper that will be used to
   * advance the time solution.  Note that it is expected that the stepper
   * will be fully initialized and ready to start taking time steps.  This
   * includes having an initial condition which defines the initial time.
   *
   * \param finalTime [in] Gives the final time that the integrator will allow
   * itself to integrate too.
   *
   * \param landOnFinalTime [in] If <tt>true</tt>, then the integrator should
   * stop exactly (within roundoff) on <tt>finalTime</tt>.  If <tt>false</tt>,
   * then the integrator can step over the final time if the stepper desires
   * it.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>!is_null(stepper)</tt>
   * <li><tt>stepper->getTimeRange().size() >= 0.0</tt>
   * <li><tt>compareTimeValues(finalTime,stepper->getTimeRange().upper()) >= 0</tt>
   * <li>If setStepper() and setInterpolationBuffer() have already been called then
   *     <tt>compareTimeValues(getTrailingInterpBuffer()->getTimeRange().upper(),stepper->getTimeRange().lower())==0</tt>.
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>this->getStepper() == stepper</tt>
   * <li><tt>compareTimeValues(this->getFwdTimeRange().lower(),stepper->getTimeRange().lower())==0</tt>
   * <li><tt>compareTimeValues(this->getFwdTimeRange().upper(),finalTime)==0</tt>
   * </ul>
   */
  virtual void setStepper(
    const RCP<StepperBase<Scalar> > &stepper,
    const Scalar &finalTime,
    const bool landOnFinalTime = true
    ) =0;

  /** \brief Get the current stepper that is set.
   *
   * \returns This function can return <tt>returnVal==null</tt> which case
   * <tt>*this</tt> is in an uninitialized state.
   */
  virtual Teuchos::RCP<const StepperBase<Scalar> > getStepper() const =0;

  /** \brief Get the current stepper that is set.
   *
   * \returns This function can return <tt>returnVal==null</tt> which case
   * <tt>*this</tt> is in an uninitialized state.
   */
  virtual Teuchos::RCP<StepperBase<Scalar> > getNonconstStepper() const =0;

  /** \brief Remove the stepper and set <tt>*this</tt> to an unitilaized
   * state.
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>is_null(this->getStepper()) == true</tt>
   * <li><tt>this->getTimeRange().isValid() == false</tt>
   * </ul>
   */
  virtual RCP<StepperBase<Scalar> > unSetStepper() =0;

  /** \brief Get values at time points both inside and outside (forward) of
   * current TimeRange.
   *
   * \param time_vec [in] Array (length <tt>n</tt>) of time points to get.
   *
   * \param x_vec [out] On output, if <tt>x_vec != 0</tt>, <tt>*x_vec</tt>
   * will be resized to <tt>n = time_vec.size()</tt> and <tt>(*x_vec)[i]</tt>
   * will be the state vector at time <tt>time_vec[i]</tt>, for
   * <tt>i=0...n-1</tt>.  This argument can be left NULL in which case it will
   * not be filled.
   *
   * \param xdot_vec [out] On output, if <tt>xdot_vec != 0</tt>,
   * <tt>*xdot_vec</tt> will be resized to <tt>n = time_vec.size()</tt> and
   * <tt>(*xdot_vec)[i]</tt> will be the state derivative vector at time
   * <tt>time_vec[i]</tt>, for <tt>i=0...n-1</tt>.  This argument can be left
   * NULL in which case it will not be filled.
   *
   * \param accuracy_vec [out] This contains an estimate of the accuracy of
   * the interpolation.  This argument can be left NULL in which case it will
   * not be filled.  If you asked for a node, this should be zero.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>range.lower() <= time_vec[i] <= range.upper()</tt>, for
   *     <tt>i=0...n-1</tt>, where <tt>range = this->getFwdTimeRange()</tt>.
   * <li><tt>time_vec</tt> must have unique and sorted values in ascending order
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li> Returns all of requested time points if no exception is thrown.
   * <li> Note that <tt>this->getTimeRange().lower()</tt> may be greater after
   *      this function returns than before this function was called!  That is
   *      why this is a non-const function!
   * </ul>
   *
   * \exception Throwns <tt>Exceptions::GetFwdPointsFailed</tt> if all of the
   * time points could not be reached for some reason (e.g. the max number of
   * time-step iterations was exceeded).
   *
   * This is a non-const version of the const function <tt>getPoints()</tt>
   * which allows the integrator class to step forward to get the points asked
   * for.
   */
  virtual void getFwdPoints(
    const Array<Scalar>& time_vec,
    Array<RCP<const Thyra::VectorBase<Scalar> > >* x_vec,
    Array<RCP<const Thyra::VectorBase<Scalar> > >* xdot_vec,
    Array<ScalarMag>* accuracy_vec
    ) =0;

  /** \brief Return the valid range of points that the integrator can
   * integrate over.
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>this->getFwdTimeRange().lower() == this->getTimeRange().lower()</tt>
   * <li><tt>this->getFwdTimeRange().upper() >= this->getTimeRange().upper()</tt>
   * </ul>
   */
  virtual TimeRange<Scalar> getFwdTimeRange() const =0;

};


// 2007/09/14: rabartl: ToDo: Below, Move these functions into a file
// Rythmos_IntegratorBaseHelpers.hpp.


/** \brief Nonmember helper function to get x at a (forward) time t.
 *
 * \relates IntegratorBase
 */
template<class Scalar> 
RCP<const Thyra::VectorBase<Scalar> >
get_fwd_x( IntegratorBase<Scalar>& integrator, const Scalar t )
{
  Array<Scalar> time_vec;
  time_vec.push_back(t);
  Array<RCP<const Thyra::VectorBase<Scalar> > > x_vec;
  integrator.getFwdPoints(time_vec,&x_vec,0,0);
  return x_vec[0];
}


/** \brief Nonmember helper function to get x and/or x_dot at s (forward)
 * time t.
 *
 * \param integrator [modified] This is the integrator that the x and x_dot
 * will be gotten from.  This is declared non-const since the state of the
 * integrator may change when getting forward points in time.
 *
 * \param t [in] The time point to get x and x_dot at.
 *
 * \param x [out] RCP to state vector to get if <tt>x!=0</tt>.
 *
 * \param x_dot [out] RCP to state derivative vector to get if
 * <tt>x_dot!=0</tt>.
 *
 * Note that <tt>*x</tt> (if <tt>x!=0</tt>) and <tt>*x_dot</tt> (if
 * <tt>x_dot!=0</tt>) should be null on input since the RCP will get set
 * internally.  Note that on output <tt>*x</tt> and <tt>*x_dot</tt> is a const
 * vector so the client can not change them.  The client should destroy (or
 * set to null) the RCP objects as soon as it is finished reading the vectors.
 * The input values of <tt>x==0</tt> and/or <tt>x_dot==0</tt> are valid which
 * means that these vectors will not be returned.
 *
 * \relates IntegratorBase
 */
template<class Scalar> 
void get_fwd_x_and_x_dot(
  IntegratorBase<Scalar>& integrator,
  const Scalar t,
  const Ptr<RCP<const Thyra::VectorBase<Scalar> > > &x,
  const Ptr<RCP<const Thyra::VectorBase<Scalar> > > &x_dot
  )
{
  Array<Scalar> time_vec;
  time_vec.push_back(t);
  Array<RCP<const Thyra::VectorBase<Scalar> > > x_vec;
  Array<RCP<const Thyra::VectorBase<Scalar> > > x_dot_vec;
  integrator.getFwdPoints(
    time_vec,
    nonnull(x) ? &x_vec : 0,
    nonnull(x_dot) ? &x_dot_vec : 0,
    0
    );
  if (nonnull(x))
    *x = x_vec[0];
  if (nonnull(x_dot))
    *x_dot = x_dot_vec[0];
}


/** \brief Depricated */
template<class Scalar> 
void get_fwd_x_and_x_dot(
  IntegratorBase<Scalar>& integrator,
  const Scalar t,
  RCP<const Thyra::VectorBase<Scalar> > *x,
  RCP<const Thyra::VectorBase<Scalar> > *x_dot
  )
{
  get_fwd_x_and_x_dot(integrator, t, ptr(x), ptr(x_dot));
}

 
} // namespace Rythmos


#endif // Rythmos_INTEGRATOR_BASE_H

