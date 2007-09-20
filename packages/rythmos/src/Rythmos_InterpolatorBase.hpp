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

#ifndef Rythmos_INTERPOLATOR_BASE_H
#define Rythmos_INTERPOLATOR_BASE_H

#include "Rythmos_DataStore.hpp"
#include "Rythmos_InterpolationBufferHelpers.hpp"


namespace Rythmos {


/** \brief Base strategy class for interpolation functionality.
 *
 * ToDo: Finish documentation!
 */
template<class Scalar> 
class InterpolatorBase
  : virtual public Teuchos::Describable
  , virtual public Teuchos::ParameterListAcceptor
  , virtual public Teuchos::VerboseObject<InterpolatorBase<Scalar> >
{
public:

  /** \brief Return if this interpolator supports cloning or not.
   *
   * If <tt>returnVal==true</tt>, then <tt>cloneInterpolator()</tt> will clone
   * <tt>*this</tt> object and return an non-null RCP.
   *
   * The default implementation of this function simply returns false.
   */
  virtual bool supportsCloning() const;

  /** \brief Clone the interpoloator object if supported.
   *
   * <b>Postconditions:</b><ul>
   * <tt>[<tt>supportsCloning()==true</tt>] <tt>returnVal != Teuchos::null</tt>
   * <tt>[<tt>supportsCloning()==false</tt>] <tt>returnVal == Teuchos::null</tt>
   * </ul>
   *
   * The default implementation returns <tt>Teuchos::null</tt> which is
   * consistent with the default implementation of <tt>supportsCloning()</tt>.
   * If this function is overridden in a base class to support cloning, then
   * <tt>supportsCloning()</tt> must be overridden to return <tt>true</tt>.
   */
  virtual Teuchos::RCP<InterpolatorBase<Scalar> > cloneInterpolator() const;

  /** \brief Perform an interpolation.
   *
   * This function must support passing node values back out directly,
   * handling when only one node value is passed in, and dealing with
   * <tt>xdot==Teuchos::null</tt>.  There is no guarantee at this time that
   * <tt>data_out</tt> will be in the same order as <tt>t_values</tt>.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>data_in</tt> must have unique time values and be sorted in ascending
   *      time order
   * <li><tt>data_in.size()>=1</tt>
   * <li>if <tt>data_in.size()==1</tt> then <tt>t_values[0] == data_in[0].time</tt>
   * <li><tt>t_values</tt> must have unique and sorted values in ascending order
   * <li><tt>data_in.front().time <= t_values[i] <= data_in.back().time</tt> for
   *     all <tt>t=0..time_values.size()-1</tt>
   * <li><tt>data_in[i].x != Teuchos::null</tt> for all <tt>i=0..data_in.size()-1</tt>
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>data_out</tt> will come out sorted in ascending time order
   * <li>if <tt>data_in[i].xdot == Teuchos::null</tt> then all t_values in the
   *    interval <tt>data_in[i-1]..data_in[i+1]</tt> will have <tt>xdot =
   *    Tuechos::null</tt>.
   * </ul>
  */
  virtual void interpolate(
    const typename DataStore<Scalar>::DataStoreVector_t &data_in,
    const Array<Scalar> &t_values,
    typename DataStore<Scalar>::DataStoreVector_t *data_out
    ) const =0;

  /** \brief Return the order of the interpolation.
   *
   * 2007/05/18: rabartl: ToDo: Define what "order" means.
   */
  virtual int order() const =0;

};


/** \relates InterplatorBase . */
template<class Scalar>
void assertBaseInterpolatePreconditions(
  const typename DataStore<Scalar>::DataStoreVector_t &data_in,
  const Array<Scalar> &t_values,
  typename DataStore<Scalar>::DataStoreVector_t *data_out
  );


// ///////////////////////////////
// Implementations


template<class Scalar>
bool InterpolatorBase<Scalar>::supportsCloning() const
{
  return false;
}


template<class Scalar>
Teuchos::RCP<InterpolatorBase<Scalar> >
InterpolatorBase<Scalar>::cloneInterpolator() const
{
  return Teuchos::null;
}


} // namespace Rythmos


template<class Scalar>
void Rythmos::assertBaseInterpolatePreconditions(
  const typename DataStore<Scalar>::DataStoreVector_t &data_in,
  const Array<Scalar> &t_values,
  typename DataStore<Scalar>::DataStoreVector_t *data_out
  )
{
  TEST_FOR_EXCEPTION(
      data_in.size()==0, std::logic_error,
      "Error, data_in.size() == 0!\n"
      );
  Array<Scalar> time_vec;
  dataStoreVectorToVector<Scalar>(data_in, &time_vec, 0, 0, 0);
  assertTimePointsAreSorted<Scalar>(time_vec);
  assertTimePointsAreSorted<Scalar>(t_values);
  if (data_in.size() == 1) {
    TEST_FOR_EXCEPTION(
      t_values.size()>1, std::logic_error,
      "Error, data_in.size() == 1, but t_values.size() > 1!\n"
      );
    TEST_FOR_EXCEPTION(
      t_values[0]!=data_in[0].time, std::logic_error,
      "Error, data_in.size) == 1, but t_values[0] = " << 
      t_values[0] << " != " << data_in[0].time << " = data_in[0].time!\n"
      );
  }
  TimeRange<Scalar> range(data_in.front().time,data_in.back().time);
  for (int i=0; i<Teuchos::as<int>(t_values.size()) ; ++i) {
    TEST_FOR_EXCEPTION(
      !range.isInRange(t_values[i]), std::logic_error,
      "Error, t_values[" << i << "] = " << t_values[i] << 
      " is not in range of data_in = " << range << "!\n"
      );
  }
  TEST_FOR_EXCEPTION(
    data_out == 0, std::logic_error,
    "Error, data_out = NULL!\n"
    );
  for (int i=0; i<Teuchos::as<int>(data_in.size()) ; ++i) {
    TEST_FOR_EXCEPTION(
      data_in[i].x == Teuchos::null, std::logic_error,
      "Error, data_in[" << i << "].x == Teuchos::null.\n"
      );
  }
}

// 2007/9/16: rabartl: ToDo: Move the above function to a new file
// Rythmos_InterpolatorBaseHelpers.hpp.


#endif //Rythmos_INTERPOLATOR_BASE_H
