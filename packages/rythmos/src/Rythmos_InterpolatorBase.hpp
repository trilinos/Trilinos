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

  /** \brief Clone the interpolator object if supported.
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

  /** \brief Store pointer to interpolation nodes
   * 
   * This function represent a persisting relationship between the
   * interpolation nodes and the interpolator.  For simple interpolators like
   * linear and Hermite, this is not needed, but for interpolators like cubic
   * splines where there is some computational work in assembling the
   * interpolant, this is important.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>nodes</tt> must have unique time values and be sorted in ascending time order
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li>If this function is called twice and <tt>nodes</tt> is a different
   * pointer than was previously called, then it is possible that the
   * interpolant will be recomputed when <tt>interpolate</tt> is next called.
   * </ul>
   */
  virtual void setNodes(
    const RCP<const typename DataStore<Scalar>::DataStoreVector_t> & nodes
    ) =0;

  /** \brief Perform an interpolation.
   *
   * This function must support passing node values back out directly,
   * handling when only one node value is passed in, and dealing with
   * <tt>xdot==Teuchos::null</tt>.  
   *
   * <b>Preconditions:</b><ul>
   * <li>if <tt>nodes->size()==1</tt> then <tt>t_values[0] == (*nodes)[0].time</tt>
   * <li><tt>t_values</tt> must have unique and sorted values in ascending order
   * <li><tt>nodes->front().time <= t_values[i] <= nodes->back().time</tt> for
   *     all <tt>t=0..time_values.size()-1</tt>
   * <li><tt>(*nodes)[i].x != Teuchos::null</tt> for all <tt>i=0..nodes->size()-1</tt>
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>data_out</tt> will come out sorted in ascending time order
   * <li>if <tt>(*nodes)[i].xdot == Teuchos::null</tt> then all t_values in the
   *    interval <tt>(*nodes)[i-1]..(*nodes)[i+1]</tt> will have <tt>xdot =
   *    Tuechos::null</tt>.
   * </ul>
  */
  virtual void interpolate(
    const Array<Scalar> &t_values,
    typename DataStore<Scalar>::DataStoreVector_t *data_out
    ) const =0;

  /** \brief Return the order of the interpolation.
   *
   * 2007/05/18: rabartl: ToDo: Define what "order" means.
   */
  virtual int order() const =0;

};



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


#endif //Rythmos_INTERPOLATOR_BASE_H
