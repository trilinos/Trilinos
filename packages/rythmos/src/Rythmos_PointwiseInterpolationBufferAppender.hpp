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

#ifndef RYTHMOS_POINTWISE_INTERPOLATION_BUFFER_APPENDER_HPP
#define RYTHMOS_POINTWISE_INTERPOLATION_BUFFER_APPENDER_HPP

#include "Rythmos_InterpolationBufferAppenderBase.hpp"
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"


namespace Rythmos {


/** \brief Concrete InterplationBufferAppender subclass that just transfers
 * notes without any regard for accuracy or order.
 */
template<class Scalar>
class PointwiseInterpolationBufferAppender
  : virtual public InterpolationBufferAppenderBase<Scalar>,
    virtual public Teuchos::ParameterListAcceptorDefaultBase
{
public:

  /** \brief . */
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;

  /** \brief Concrete implementation that simply copies the nodal points
   * between the interpolation buffers.
   */
  void append(
    const InterpolationBufferBase<Scalar>& interpBuffSource, 
    const TimeRange<Scalar>& range,
    const Ptr<InterpolationBufferBase<Scalar> > &interpBuffSink 
    );

  /** \name Overridden from Teuchos::Describable */
  //@{

  /** \brief . */
  void describe(
    Teuchos::FancyOStream &out,
    const Teuchos::EVerbosityLevel verbLevel
    ) const;

  //@}

  /** \name Overridden from ParameterListAcceptorDefaultBase */
  //@{

  /** \brief . */
  void setParameterList(const RCP<ParameterList> &paramList);

  /** \brief . */
  RCP<const ParameterList> getValidParameters() const;

  //@}

};



/** \brief Nonmember constructor function.
 *
 * \relates PointwiseInterpolationBufferAppender
 */
template<class Scalar>
RCP<PointwiseInterpolationBufferAppender<Scalar> >
pointwiseInterpolationBufferAppender()
{
  return Teuchos::rcp(new PointwiseInterpolationBufferAppender<Scalar>);
}


//
// Implementations
//


template<class Scalar>
void PointwiseInterpolationBufferAppender<Scalar>::append(
  const InterpolationBufferBase<Scalar>& interpBuffSource, 
  const TimeRange<Scalar>& range,
  const Ptr<InterpolationBufferBase<Scalar> > &interpBuffSink 
  ) 
{
#ifdef RYTHMOS_DEBUG
  this->assertAppendPreconditions(interpBuffSource,range,*interpBuffSink);
#endif // RYTHMOS_DEBUG

  Array<Scalar> time_vec_in;
  interpBuffSource.getNodes(&time_vec_in);

  Array<Scalar> time_vec;
  selectPointsInTimeRange(&time_vec,time_vec_in,range);
  // 2007/12/05: rabrtl: ToDo: Make the output argument last!

  Array<RCP<const Thyra::VectorBase<Scalar> > > x_vec;
  Array<RCP<const Thyra::VectorBase<Scalar> > > xdot_vec;
  Array<ScalarMag> accuracy_vec;
  interpBuffSource.getPoints(time_vec, &x_vec, &xdot_vec, &accuracy_vec);

  interpBuffSink->addPoints(time_vec, x_vec, xdot_vec);

}


template<class Scalar>
void PointwiseInterpolationBufferAppender<Scalar>::describe(
  Teuchos::FancyOStream &out,
  const Teuchos::EVerbosityLevel verbLevel
  ) const
{
  using Teuchos::as;
  if (
    (as<int>(verbLevel) == as<int>(Teuchos::VERB_DEFAULT))
    || (as<int>(verbLevel) >= as<int>(Teuchos::VERB_LOW))
    )
  {
    out << this->description() << std::endl;
  }
}


template<class Scalar>
void PointwiseInterpolationBufferAppender<Scalar>::setParameterList(
  const RCP<ParameterList> &paramList
  )
{
  TEST_FOR_EXCEPT( is_null(paramList) );
  paramList->validateParameters(*this->getValidParameters());
  Teuchos::readVerboseObjectSublist(&*paramList,this);
  setMyParamList(paramList);
}


template<class Scalar>
RCP<const ParameterList>
PointwiseInterpolationBufferAppender<Scalar>::getValidParameters() const
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


#endif //RYTHMOS_POINTWISE_INTERPOLATION_BUFFER_APPENDER_HPP
