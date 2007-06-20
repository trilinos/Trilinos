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

#ifndef Rythmos_INTERPOLATION_BUFFER_BASE_H
#define Rythmos_INTERPOLATION_BUFFER_BASE_H

#include "Rythmos_Types.hpp"
#include "Thyra_VectorBase.hpp"
#include "Teuchos_Describable.hpp"
#include "Teuchos_ParameterListAcceptor.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_implicit_cast.hpp"
#include "Teuchos_Assert.hpp"


namespace Rythmos {


/** \brief Represent a time range.
 *
 * The compiler-generated default constructor, copy constructor, and
 * assignment operators are allowed and perform correctly.
 *
 * ToDo: Put in checks for the range if needed.
 */
template<class Time>
class TimeRange {
public:
  /** \brief Construct an invalid range. */
  TimeRange()
    : lower_(0.0), upper_(-1.0)
    {}
  /** \brief Construct a valid range. */
  TimeRange( const Time lower, const Time upper )
    : lower_(lower), upper_(upper)
    {
      TEUCHOS_ASSERT( lower_ <= upper_ ); 
    }
  /** \brief . */
  bool isValid() const { return (lower_ <= upper_); }
  /** \brief . */
  Time lower() const { return lower_; }
  /** \brief . */
  Time upper() const { return upper_; }
private:
  Time lower_;
  Time upper_;
};


/** \brief Nonmember constructor.
 *
 * \relates TimeRange
 */
template<class Time>
inline
TimeRange<Time> timeRange(const Time lower, const Time upper)
{
  return TimeRange<Time>(lower,upper);
}


/** \brief Nonmember constructor.
 *
 * \relates TimeRange
 */
template<class Time>
inline
TimeRange<Time> invalidTimeRange()
{
  return TimeRange<Time>();
}


/** \brief Base class for an interpolation buffer.
 *
 * An interpolation buffer represents the data and the ability to represent
 * and interpolate the values of some transient solution <tt>x</tt> and
 * <tt>x_dot</tt> over a contiguous range of time.
 *
 * 2007/05/21: rabartl ToDo: We really need to discuss the concepts that are
 * represented in this interface and we especially need to define what
 * "accuracy" is and how it is related to what points are added and retrieved.
 *
 * \section Rythmos_InterpolationBufferBase_Definitions_sec Definitions
 *
 * <ul>
 *
 * <li><b>Accuracy</b>: The approximate maximum error (in what Norm???) 
 * between any two node points for the function being approximated.  The
 * accuracy of an interpolation buffer will often come from the local
 * trucation error estimate from an accuracy-controlling time step algorithm.
 *
 * <li><b>Order</b>: The degree of polynomical that can be represented exactly
 * by the buffer interface.  For example, a second-order
 * (i.e. <tt>order==2</tt>) interpolation buffer can exactly represent any
 * polynomical up to degree two.
 *
 * </ul>
 *
 * ToDo: Finish documentation!
 */
template<class Scalar> 
class InterpolationBufferBase 
  : virtual public Teuchos::Describable
  , virtual public Teuchos::ParameterListAcceptor
  , virtual public Teuchos::VerboseObject<InterpolationBufferBase<Scalar> >
{
public:

  /** \brief . */
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;

  /** \brief Return the space for <tt>x</tt> and <tt>x_dot</tt>.
   *
   * This space can be used to create vectors for calling <tt>setPoints()</tt>
   * for instance and is also useful for writing unit testing software.
   *
   * Also note that this space may not be the same as the space returned from
   * <tt>StepperBase::getModel()->get_x_sapce()</tt> in some concrete
   * <tt>StepperBase</tt> subclasses.
   */
  virtual RCP<const Thyra::VectorSpaceBase<Scalar> >
  get_x_space() const
    {
      return Teuchos::null;
    }
  // 2007/05/17: rabartl: ToDo: Remove this default implementation and
  // propogate these changes down to the subclasses.  Every interpolation
  // buffer object should be able to implement this function.  I just don't
  // have time to propogate this yet but I will soon.

  /** \brief Add points to the buffer.
   *
   * \param time_vec
   *          [in] Array (length n) of time points.
   * \param x_vec
   *          [in] Array (length n) of state vectors at each time point in
   *          <tt>time_vec</tt>.  Specifically, <tt>x_vec[i]</tt> is the state
   *          vector at time <tt>time_vec[i]</tt>, for <tt>i=0...n-1</tt>.
   *          The RCP for the vectors may or may not be copied by
   *          <tt>*this</tt>.  An implementation is not required to copy the
   *          RCP's to the vector objects but instead might just use the
   *          vectors or do a deep copy.
   * \param xdot_vec
   *          [in] Array (length n) of state time differential vector at each
   *          time point in <tt>time_vec</tt>.  Specifically,
   *          <tt>xdot_vec[i]</tt> is the state differential vector at time
   *          <tt>time_vec[i]</tt>, for <tt>i=0...n-1</tt>.  The RCP for the
   *          vectors may or may not be copied by <tt>*this</tt>.  An
   *          implementation is not required to copy the RCP's to the vector
   *          objects but instead might just use the vectors or do a deep
   *          copy.
   * \param accuracy_vec
   *          [in] Array (length n) that gives the approximate accuracy of
   *          each of the input values for <tt>x</tt> and <tt>x_dot</tt>.  See the
   *          definition of accruacy given above.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>time_vec.size()!=0</tt>
   * <li><tt>time_vec.size()==x_vec.size()</tt>
   * <li><tt>time_vec.size()==xdot_vec.size()</tt>
   * <li><tt>time_vec.size()==accuracy_vec.size()</tt>
   * </ul>
   *
   * ToDo: Verify that the above are the correct preconditions?  For example,
   * is it okay for <tt>xdot_vec[i].get() == 0 </tt> for any <tt>i</tt>?
   *
   * <b>Postconditions:</b><ul>
   * <li>What are the post conditions?
   * </ul>
   *
   * ToDo: What is the role of the return value?  Is an implementation allowed
   * to simply ignore some time points or all of the time points?  How is this
   * helpful to the client?  How can we tighten up the specification of the
   * behavior so that incorrect usage results in an exception (in at least
   * debug mode)?
   */
  virtual bool setPoints(
    const Array<Scalar>& time_vec,
    const Array<RCP<const Thyra::VectorBase<Scalar> > >& x_vec,
    const Array<RCP<const Thyra::VectorBase<Scalar> > >& xdot_vec,
    const Array<ScalarMag> & accuracy_vec
    ) = 0;

  /** \brief Fill data in from another interpolation buffer.
   *
   * \param  range
   *           [in] The time range in <tt>interpBuffer</tt> that will be
   *           extracted into <tt>*this</tt> interpolation buffer.
   * \param  interpBuffer
   *           [in] The interpolation buffer that will be queried to get
   *           interpolated values to put into <tt>*this</tt> interpolation
   *           buffer.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>(range.lower() == this->getTimeRange().upper()) ||
   *         (range.upper() == this->getTimeRange().lower())</tt>
   * <li><tt>interpBuffer.getTimeRange().lower() <= range.lower()</tt> 
   * <li><tt>range.upper() <= interpBuffer.getTimeRange().upper()</tt> 
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li>What are the post conditions?
   * </ul>
   *
   * The time points in the interval [range.lower(),range.upper()] will be
   * inserted into the interpolation buffer.
   *
   * ToDo: What is the role of the return value?  Is an implementation allowed
   * to simply ignore parts of the time range?  How is this helpful to the
   * client?  How can we tighten up the specification of the behavior so that
   * incorrect usage results in an exception (in at least debug mode)?
   */
  virtual bool setRange(
    const TimeRange<Scalar>& range,
    const InterpolationBufferBase<Scalar>& interpBuffer
    ) = 0;

  /** \brief Return the range of time values where interpolation calls can be
   * performed.
   *
   * A return value of <tt>returnVal.isValid()==false</tt> means that there is
   * no time range for which interpolation can be performed.  Otherwise,
   * <tt>returnVal</tt> gives the time range for which state information can
   * be returned.
   */
  virtual TimeRange<Scalar> getTimeRange() const = 0;

  /** \brief Get values from the buffer at different time points.
   *
   * \param time_vec
   *          [in] Array (length n) of time points to get.
   * \param x_vec
   *          [out] On output, if <tt>x_vec != 0</tt>, <tt>*x_vec</tt> will be
   *          resized to <tt>n = time_vec.size()</tt> and <tt>(*x_vec)[i]</tt>
   *          will be the state vector at time <tt>time_vec[i]</tt>, for
   *          <tt>i=0...n-1</tt>.  This argument can be left NULL in which
   *          case it will not be filled.
   * \param xdot_vec
   *          [out] On output, if <tt>xdot_vec != 0</tt>, <tt>*xdot_vec</tt>
   *          will be resized to <tt>n = time_vec.size()</tt> and
   *          <tt>(*xdot_vec)[i]</tt> will be the state derivative vector at
   *          time <tt>time_vec[i]</tt>, for <tt>i=0...n-1</tt>.  This
   *          argument can be left NULL in which case it will not be filled.
   * \param accuracy_vec
   *          [out] On output, if <tt>accuracy_vec != 0</tt>,
   *          <tt>*accuracy_vec</tt> will be resized to <tt>n =
   *          time_vec.size()</tt> and <tt>(*accuracy_vec)[i]</tt> will give
   *          the accuracy of the retrieved interpolants in <tt>x_vec</tt> and
   *          <tt>xdot_vec</tt>.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>range.lower() <= time_vec[i] <= range.upper()</tt>, for
   *     <tt>i=0...n-1</tt>, where <tt>range = this->getTimeRange()</tt>.
   * </ul>
   *
   * 2007/06/08: rabartl: ToDo: Should we require that the user sort the time
   * point?  this will make the implementation by the subclasses much easier.
   * I think this is a reasonable thing to request.  If we want to allow the
   * user to pass in an unsorted time array, then we can provide an
   * unsortedGetPoints(...) helper functions that will take an unsorted
   * time_vec array, sort it, call getPoints(...) and then put the returned
   * values back in the right place.
   *
   * 2007/06/08: rabartl: ToDo: Perhaps we need to make getPoints(...) a
   * non-const member function since it might actually chagne the state of
   * <tt>*this</tt>!  In the case of the integrator class, calling
   * getPoints(...)  can actually make it so that the next call to
   * getPoints(...) might fail if old buffer memory is lost.  Do we need both
   * a const version (i.e. only return what I already can give you) and a
   * nonconst version (i.e. I can get the points you need)?  We need to
   * rethink the concept of a time range to differentiate between what *this
   * can give you without changing its state vs. a time range of what the
   * buffer could give you if you asked for it.  Right now this is all very
   * squishy.
   *
   * rabartl: ToDo: What is the role of the return value?  Is it allowed for
   * the implementation to refuse to return values for only certain times or
   * for all times?  If the time values fall in the allowed range, should not
   * an implementation be forced to return the interpolated values?  How can
   * we tighten up the specification of this function?
   */
  virtual bool getPoints(
    const Array<Scalar>& time_vec,
    Array<RCP<const Thyra::VectorBase<Scalar> > >* x_vec,
    Array<RCP<const Thyra::VectorBase<Scalar> > >* xdot_vec,
    Array<ScalarMag>* accuracy_vec
    ) const = 0;

  /** \brief Get interpolation nodes.
   *
   * This function will return the time points where actual data is stored.
   * This information can be used to get the actual nodal values themselves
   * using the <tt>getPoints()</tt> function.
   *
   * 2007/05/21: rabartl: ToDo: What is the role of the return value?  Do some
   * steppers simply not have any nodes?  If so, perhaps we need a query
   * function like <tt>hasInterpolationNodes()</tt> that will tell the client
   * this.
   */
  virtual bool getNodes(Array<Scalar>* time_vec) const = 0;

  /** \brief Remove interpolation nodes.
   *
   * ToDo: What is the role of the return value?  Do some steppers simply not
   * have any nodes?  If so, then we should simply thrown an exception if a
   * client tries to remove nodes that don't exist.
   */
  virtual bool removeNodes(Array<Scalar>& time_vec) = 0;

  /** \brief Get order of interpolation.
   *
   * See the above section \ref Rythmos_InterpolationBufferBase_Definitions_sec
   * for the definition of "order".
   */
  virtual int getOrder() const = 0;

};


/** \brief Get a single point <tt>x(t)</tt> from an interpolation buffer.
 *
 * \relates InterpolationBufferBase
 */
template<class Scalar>
RCP<const Thyra::VectorBase<Scalar> >
get_x( const InterpolationBufferBase<Scalar> &interpBuffer, const Scalar &t )
{
  using Teuchos::implicit_cast;
  Array<Scalar> time_vec;
  time_vec.push_back(t);
  Array<RCP<const Thyra::VectorBase<Scalar> > > x_vec;
  interpBuffer.getPoints(time_vec,&x_vec,0,0);
  TEUCHOS_ASSERT( 1 == implicit_cast<int>(x_vec.size()) );
  return x_vec[0];
}


} // namespace Rythmos


#endif //Rythmos_INTERPOLATION_BUFFER_BASE_H
