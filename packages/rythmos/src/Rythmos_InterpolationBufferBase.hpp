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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef Rythmos_INTERPOLATION_BUFFER_BASE_H
#define Rythmos_INTERPOLATION_BUFFER_BASE_H

#include "Rythmos_Types.hpp"
#include "Rythmos_TimeRange.hpp"

#include "Thyra_VectorBase.hpp"

#include "Teuchos_Describable.hpp"
#include "Teuchos_ParameterListAcceptor.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_implicit_cast.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_as.hpp"


namespace Rythmos {


/** \brief Base class for an interpolation buffer.
 *
 * An interpolation buffer represents the data and the ability to represent
 * and interpolate the values of some transient solution <tt>x</tt> and
 * <tt>x_dot</tt> over a contiguous range of time.
 *
 * \section Rythmos_InterpolationBufferBase_Definitions_sec Definitions
 *
 * <ul>
 *
 * <li><b>Order</b>: The degree of polynomial that can be represented exactly
 * by the buffer interface.  For example, a second-order
 * (i.e. <tt>order==2</tt>) interpolation buffer can exactly represent any
 * polynomial up to degree two.
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
   * This space can be used to create vectors for calling <tt>addPoints()</tt>
   * for instance and is also useful for writing unit testing software.
   *
   * Also note that this space may not be the same as the space returned from
   * <tt>StepperBase::getModel()->get_x_space()</tt> in some concrete
   * <tt>StepperBase</tt> subclasses.
   */
  virtual RCP<const Thyra::VectorSpaceBase<Scalar> >
  get_x_space() const =0;

  /** \brief Add points to the buffer.
   *
   * \param time_vec [in] Array (length n) of time points.
   *
   * \param x_vec [in] Array (length n) of state vectors at each time point in
   * <tt>time_vec</tt>.  Specifically, <tt>x_vec[i]</tt> is the state vector
   * at time <tt>time_vec[i]</tt>, for <tt>i=0...n-1</tt>.  The RCP for the
   * vectors may or may not be copied by <tt>*this</tt>.  An implementation is
   * not required to copy the RCP's to the vector objects but instead might
   * just use the vectors or do a deep copy.
   *
   * \param xdot_vec [in] Array (length n) of state time differential vector
   * at each time point in <tt>time_vec</tt>.  Specifically,
   * <tt>xdot_vec[i]</tt> is the state differential vector at time
   * <tt>time_vec[i]</tt>, for <tt>i=0...n-1</tt>.  The RCP for the vectors
   * may or may not be copied by <tt>*this</tt>.  An implementation is not
   * required to copy the RCP's to the vector objects.  The implementation may
   * copy the actually vectors or might just perform a shallow copy by copying
   * the RCP objects.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>time_vec.size()!=0</tt>
   * <li><tt>time_vec</tt> must have unique and sorted values in ascending order
   * <li><tt>time_vec.size()==x_vec.size()</tt>
   * <li><tt>time_vec.size()==xdot_vec.size()</tt>
   * <li><tt>x_vec[i] != null</tt> for <tt>i = 0...n-1</tt>
   * <li><tt>xdot_vec[i] != null</tt> for <tt>i = 0...n-1</tt>
   * <li><tt>getTimeRange().isInRange(time_vec[i]) == false</tt>, for <tt>i=0...n-1</tt>
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>getTimeRange().isInRange(time_vec[i]) == true</tt>, for <tt>i=0...n-1</tt>
   * <li><tt>getTimeRange().lower()</tt> may increase after the call.
   * </ul>
   */
  /* 11/24/08 tscoffe:  Proposed new interface for addPoints
   * virtual void addPoints(
   *   const ArrayView<const Scalar>& time_vec,
   *   const ArrayView<const Ptr<const Thyra::VectorBase<Scalar> > >& x_vec,
   *   const ArrayView<const Ptr<const Thyra::VectorBase<Scalar> > >& xdot_vec
    ) = 0;
   */

  virtual void addPoints(
    const Array<Scalar>& time_vec,
    const Array<RCP<const Thyra::VectorBase<Scalar> > >& x_vec,
    const Array<RCP<const Thyra::VectorBase<Scalar> > >& xdot_vec
    ) = 0;

  /*
  virtual void addPoints(
    const ArrayView<const Scalar>& time_vec,
    const ArrayView<const Ptr<const VectorBase<Scalar> > >& x_vec,
    const ArrayView<const Ptr<const VectorBase<Scalar> > >& xdot_vec
    ) =0;
  */


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
   * \param time_vec [in] Array (length n) of time points to get.
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
   * the interpolation.  If you asked for a node, this should be zero.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>range.isInRange(time_vec[i])</tt>, for
   *     <tt>i=0...n-1</tt>, where <tt>range = this->getTimeRange()</tt>.
   * <li><tt>time_vec</tt> must have unique and sorted values in ascending order
   * </ul>
   *
   */
  /* 11/24/08 tscoffe:  Proposed new interface for getPoints
   * virtual void getPoints(
   *   const ArrayView<const Scalar>& time_vec,
   *   const ArrayView<const Ptr<Thyra::VectorBase<Scalar> > >& x_vec,
   *   const ArrayView<const Ptr<Thyra::VectorBase<Scalar> > >& xdot_vec,
   *   const ArrayView<ScalarMag>& accuracy_vec
   *   ) const = 0;
   */
  virtual void getPoints(
    const Array<Scalar>& time_vec,
    Array<RCP<const Thyra::VectorBase<Scalar> > >* x_vec,
    Array<RCP<const Thyra::VectorBase<Scalar> > >* xdot_vec,
    Array<ScalarMag>* accuracy_vec
    ) const = 0;

  /*
  virtual void getPoints(
    const ArrayView<const Scalar>& time_vec,
    const ArrayView<const Ptr<VectorBase<Scalar> > >& x_vec,
    const ArrayView<const Ptr<VectorBase<Scalar> > >& xdot_vec,
    const ArrayView<ScalarMag>& accuracy_vec
    ) const = 0;
  */

  /** \brief Get interpolation nodes.
   *
   * This function will return the time points where actual data is stored.
   * This information can be used to get the actual nodal values themselves
   * using the <tt>getPoints()</tt> function.
   *
   * This function may return no nodes at all, and yet still have a valid
   * timeRange.
   */
  // 11/24/08 tscoffe:  Proposal:  get rid of "getNodes" in abstract base interface
  virtual void getNodes(Array<Scalar>* time_vec) const = 0;

  /** \brief Remove nodes from the interpolation buffer.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>time_vec</tt> must have unique and sorted values in ascending order
   * <li><tt>range.isInRange(time_vec[i])</tt>, for
   *     <tt>i=0..n-1</tt>, where <tt>range = this->getTimeRange()</tt>.
   * <li>Each <tt>time_vec[i]</tt> must equal a node in the interpolation buffer.
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li>Each time_vec[i] is no longer a node in the interpolation buffer.
   * </ul>
   */
  // 11/24/08 tscoffe:  Proposal:  get rid of "removeNodes" in abstract base interface
  virtual void removeNodes(Array<Scalar>& time_vec) =0;


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

/* 11/24/08 tscoffe:  Proposed new get_x function:
 * template<class Scalar>
 * void get_x( const InterpolationBufferBase<Scalar> &interpBuffer, const Scalar &t, const Ptr<Thyra::VectorBase<Scalar> >& x_out )
 * This will copy data into your vector and it won't be responsible for allocating new memory
 */


/** \brief Get a single point <tt>xdot(t)</tt> from an interpolation buffer.
 *
 * \relates InterpolationBufferBase
 */
template<class Scalar>
RCP<const Thyra::VectorBase<Scalar> >
get_xdot( const InterpolationBufferBase<Scalar> &interpBuffer, const Scalar &t )
{
  using Teuchos::implicit_cast;
  Array<Scalar> time_vec;
  time_vec.push_back(t);
  Array<RCP<const Thyra::VectorBase<Scalar> > > xdot_vec;
  interpBuffer.getPoints(time_vec,0,&xdot_vec,0);
  TEUCHOS_ASSERT( 1 == implicit_cast<int>(xdot_vec.size()) );
  return xdot_vec[0];
}

/* 11/24/08 tscoffe:  Proposed new get_xdot function
 * template<class Scalar>
 * void get_xdot( const InterpolationBufferBase<Scalar> &interpBuffer, const Scalar &t, const Ptr<Thyra::VectorBase<Scalar> >& xdot_out )
 * This will copy data into your vector and it won't be responsible for allocating new memory
 */

/** \brief Nonmember helper function to get x and x_dot at t.
 *
 * \relates InterpolationBufferBase
 */
template<class Scalar>
void get_x_and_x_dot(
  const InterpolationBufferBase<Scalar> &interpBuffer,
  const Scalar t,
  const Ptr<RCP<const Thyra::VectorBase<Scalar> > > &x,
  const Ptr<RCP<const Thyra::VectorBase<Scalar> > > &x_dot
  )
{
  Array<Scalar> time_vec;
  time_vec.push_back(t);
  Array<RCP<const Thyra::VectorBase<Scalar> > > x_vec;
  Array<RCP<const Thyra::VectorBase<Scalar> > > x_dot_vec;
  interpBuffer.getPoints(
    time_vec,
    nonnull(x) ? &x_vec : 0,
    nonnull(x_dot) ? &x_dot_vec : 0,
    0
    );
  if (nonnull(x)) *x = x_vec[0];
  if (nonnull(x_dot)) *x_dot = x_dot_vec[0];
}

/** \brief Deprecated. */
template<class Scalar>
void get_x_and_x_dot(
  const InterpolationBufferBase<Scalar> &interpBuffer,
  const Scalar t,
  RCP<const Thyra::VectorBase<Scalar> > *x,
  RCP<const Thyra::VectorBase<Scalar> > *x_dot
  )
{
  get_x_and_x_dot(interpBuffer, t, Teuchos::ptr(x), Teuchos::ptr(x_dot));
}

/* 11/24/08 tscoffe:  Proposed new get_x_and_xdot function:
 * template<class Scalar>
 * void get_x_and_xdot( const InterpolationBufferBase<Scalar> &interpBuffer,
 *   const Scalar &t,
 *   const Ptr<Thyra::VectorBase<Scalar> >& x_out,
 *   const Ptr<Thyra::VectorBase<Scalar> >& xdot_out )
 * This will copy data into your vector and it won't be responsible for allocating new memory
 */

} // namespace Rythmos


#endif //Rythmos_INTERPOLATION_BUFFER_BASE_H
