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

#ifndef RYTHMOS_SMART_INTERPOLATION_BUFFER_APPENDER_HPP
#define RYTHMOS_SMART_INTERPOLATION_BUFFER_APPENDER_HPP

#include "Rythmos_InterpolationBufferAppenderBase.hpp"
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"


namespace Rythmos {


/** \brief Smart interplation buffer class. */
template<class Scalar>
class SmartInterpolationBufferAppender
  : virtual public InterpolationBufferAppenderBase<Scalar>,
    virtual public Teuchos::ParameterListAcceptorDefaultBase
{
  public:
    /** \brief Concrete implementation that attempts to use the order of
     * interpolation between the two interpolation buffers to be a bit smarter
     * about copying data between them.
     */
    void append(
        const InterpolationBufferBase<Scalar>& interpBuffSource, 
        const TimeRange<Scalar>& range,
        const Ptr<InterpolationBufferBase<Scalar> > &interpBuffSink 
        );

    /** \name Overridden from Teuchos::ParameterListAcceptorDefaultBase */
    //@{

    /** \brief . */
    void setParameterList(RCP<Teuchos::ParameterList> const& paramList);

    /** \brief . */
    RCP<const Teuchos::ParameterList> getValidParameters() const;

    //@}
};


//
// Implementations
//


template<class Scalar>
void SmartInterpolationBufferAppender<Scalar>::append(
  const InterpolationBufferBase<Scalar>& interpBuffSource, 
  const TimeRange<Scalar>& range,
  const Ptr<InterpolationBufferBase<Scalar> > &interpBuffSink
    ) 
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    true, std::logic_error,
    "This class has never been tested before and should not be used\n"
    "until it is"
    );
  // 2007/12/05: rabartl: This code has not been tested so don't use this
  // until a test has been writen for this!
#ifdef HAVE_RYTHMOS_DEBUG
  this->assertAppendPreconditions(interpBuffSource,range,*interpBuffSink);
#endif // HAVE_RYTHMOS_DEBUG
  if (interpBuffSink->getOrder() >= interpBuffSource.getOrder()) {
    // The incoming interpolation buffer's order of interpolation is lower than
    // the base interpolation buffer's order of interpolation.  In this case,
    // we just copy the data over.
    PointwiseInterpolationBufferAppender<Scalar> defaultAppender;
    defaultAppender.append(interpBuffSink,interpBuffSource,range);
  } else {
    // In this case, the incoming interpolation buffer's order of interpolation
    // is higher than the base interpolation buffer's, so we'll ask it to
    // interpolate points before inserting into the base interpolation buffer.
    TEUCHOS_TEST_FOR_EXCEPTION(
        true,std::logic_error,
        "Error, the smart interpolation buffer appender is not implemented\n"
        "for appending interpolation buffers with higher order interpolation\n"
        "into interpolation buffers with a lower order of interpolation yet!"
        );
    // 08/09/07 tscoffe:  Note:  you can't use selectPointsInTimeRange
    // immediately because you may need to interpolate points before the first
    // node inside the range.  I.e. interpBuffSource.getNodes = [... , range.lower(), ... , range.upper(), ... ]
  }
}

template<class Scalar>
void SmartInterpolationBufferAppender<Scalar>::setParameterList(RCP<Teuchos::ParameterList> const& paramList)
{
  TEUCHOS_TEST_FOR_EXCEPT( is_null(paramList) );
  paramList->validateParameters(*this->getValidParameters());
  setMyParamList(paramList);
  Teuchos::readVerboseObjectSublist(&*paramList,this);
}

template<class Scalar>
RCP<const Teuchos::ParameterList> SmartInterpolationBufferAppender<Scalar>::getValidParameters() const
{
  static RCP<Teuchos::ParameterList> validPL;
  if (is_null(validPL)) {
    RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    Teuchos::setupVerboseObjectSublist(&*pl);
    validPL = pl;
  }
  return (validPL);
}

} // namespace Rythmos


#endif //RYTHMOS_SMART_INTERPOLATION_BUFFER_APPENDER_HPP
