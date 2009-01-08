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

#ifndef RYTHMOS_INTERPOLATION_BUFFER_APPENDER_BASE_HPP
#define RYTHMOS_INTERPOLATION_BUFFER_APPENDER_BASE_HPP

#include "Rythmos_InterpolationBufferBase.hpp"
#include "Rythmos_Types.hpp"
#include "Thyra_VectorBase.hpp"
#include "Teuchos_Describable.hpp"
#include "Teuchos_ParameterListAcceptor.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_implicit_cast.hpp"
#include "Teuchos_Assert.hpp"


namespace Rythmos {


/** \brief Base class for strategy objects that append data from one
 * InterplationBufferBase object to another.
 */
template<class Scalar>
class InterpolationBufferAppenderBase
  : virtual public Teuchos::Describable
  , virtual public Teuchos::ParameterListAcceptor
  , virtual public Teuchos::VerboseObject<InterpolationBufferAppenderBase<Scalar> >
{
public:

  /** \brief Append or Prepend data from one interpolation buffer into another.
   *
   * \param interpBuffSink [in/out] The interpolation buffer that will recieve the
   * data from <tt>interpBuffSource</tt> interpolation buffer.
   *
   * \param interpBuffSource [in] The interpolation buffer that will be queried to get
   * interpolated values to put into <tt>interpBuffSink</tt> interoplation buffer.
   *
   * \param range [in] The time range in <tt>interpBuffSource</tt> that will be
   * converted into <tt>interpBuffSink</tt> interpolation buffer.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>(range.lower() == interpBuffSink->getTimeRange().upper()) ||
   *         (range.upper() == interpBuffSink->getTimeRange().lower())</tt>
   * <li><tt>interpBuffSource.getTimeRange().lower() <= range.lower()</tt>
   * <li><tt>range.upper() <= interpBuffSource.getTimeRange().upper()</tt>
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li>if prior to the call <tt>range.lower() == interpBuffSink->getTimeRange().upper()</tt> 
   *     then after the call <tt>interpBuffSink->getTimeRange().upper() == range.upper()</tt>
   * <li>if prior to the call <tt>range.upper() == interpBuffSink->getTimeRange().lower()</tt>
   *     then after the call <tt>interpBuffSink->getTimeRange().lower() == range.lower()</tt>
   * </ul>
   */
  virtual void append(
    const InterpolationBufferBase<Scalar>& interpBuffSource,
    const TimeRange<Scalar>& range,
    const Ptr<InterpolationBufferBase<Scalar> > &interpBuffSink
    ) =0;

protected:

  /** \brief . */
  void assertAppendPreconditions(
    const InterpolationBufferBase<Scalar>& interpBuffSource,
    const TimeRange<Scalar>& range,
    const InterpolationBufferBase<Scalar>& interpBuffSink
    ) const;

};


template<class Scalar>
void InterpolationBufferAppenderBase<Scalar>::assertAppendPreconditions(
  const InterpolationBufferBase<Scalar>& interpBuffSource, 
  const TimeRange<Scalar>& range,
  const InterpolationBufferBase<Scalar>& interpBuffSink
  ) const
{
  // If the time range of interpBuffSink is invalid, then its just empty
  if (interpBuffSink.getTimeRange().isValid()) {
    TEST_FOR_EXCEPTION(
      ( compareTimeValues(range.lower(),interpBuffSink.getTimeRange().upper()) != 0 &&
        compareTimeValues(range.upper(),interpBuffSink.getTimeRange().lower()) != 0    ),
      std::logic_error, 
      "Error, import range = [" << range.lower() << "," << range.upper() << "] is not an append nor a prepend "
      "of the base range = [" << interpBuffSink.getTimeRange().lower() << "," << interpBuffSink.getTimeRange().upper() << "] "
      "interpolation buffer.\n"
      );
  }
  TEST_FOR_EXCEPTION(
    compareTimeValues(range.lower(),interpBuffSource.getTimeRange().lower())<0,
    std::logic_error,
    "Error, append range's lower bound = " << range.lower() << " does not sit inside incoming"
    " interpolation buffer's time range = "
    "[" << interpBuffSource.getTimeRange().lower() << "," << interpBuffSource.getTimeRange().upper() << "].\n"
    );
  TEST_FOR_EXCEPTION(
    compareTimeValues(interpBuffSource.getTimeRange().upper(),range.upper())<0,
    std::logic_error,
    "Error, append range's upper bound = " << range.upper() << "does not sit inside incoming"
    " interpolation buffer's time range = "
    "[" << interpBuffSource.getTimeRange().lower() << "," << interpBuffSource.getTimeRange().upper() << "].\n"
    );
}


} // namespace Rythmos


#endif //RYTHMOS_INTERPOLATION_BUFFER_APPENDER_BASE_HPP
