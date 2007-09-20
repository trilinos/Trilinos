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

#ifndef Rythmos_INTERPOLATION_BUFFER_APPENDER_H
#define Rythmos_INTERPOLATION_BUFFER_APPENDER_H

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

template<class Scalar>
class InterpolationBufferAppenderBase
  : virtual public Teuchos::Describable
  , virtual public Teuchos::ParameterListAcceptor
  , virtual public Teuchos::VerboseObject<InterpolationBufferAppenderBase<Scalar> >
{
  public:

    /** \brief Append or Prepend data from one interpolation buffer into another.
     *
     * \param IB_base
     *          [in/out] The interpolation buffer that will recieve the data from
     *          <tt>IB_in</tt> interpolation buffer.
     * \param IB_in
     *          [in] The interpolation buffer that will be queried to get
     *          interpolated values to put into <tt>IB_base</tt> interoplation
     *          buffer.
     * \param range
     *          [in] The time range in <tt>IB_in</tt> that will be converted
     *          into <tt>IB_base</tt> interpolation buffer.
     *
     * <b>Preconditions:</b><ul>
     * <li><tt>(range.lower() == IB_base->getTimeRange().upper()) ||
     *         (range.upper() == IB_base->getTimeRange().lower())</tt>
     * <li><tt>IB_in.getTimeRange().lower() <= range.lower()</tt>
     * <li><tt>range.upper() <= IB_in.getTimeRange().upper()</tt>
     * </ul>
     *
     * <b>Postconditions:</b><ul>
     * <li>if prior to the call <tt>range.lower() == IB_base->getTimeRange().upper()</tt> 
     *     then after the call <tt>IB_base->getTimeRange().upper() == range.upper()</tt>
     * <li>if prior to the call <tt>range.upper() == IB_base->getTimeRange().lower()</tt>
     *     then after the call <tt>IB_base->getTimeRange().lower() == range.lower()</tt>
     * </ul>
     */
    virtual void import(
        InterpolationBufferBase<Scalar>* IB_base, 
        const InterpolationBufferBase<Scalar>& IB_in, 
        const TimeRange<Scalar>& range
        ) =0;
  protected:
    void assertImportPreconditions(
        const InterpolationBufferBase<Scalar>& IB_base, 
        const InterpolationBufferBase<Scalar>& IB_in, 
        const TimeRange<Scalar>& range
        ) const;
};

template<class Scalar>
void InterpolationBufferAppenderBase<Scalar>::assertImportPreconditions(
        const InterpolationBufferBase<Scalar>& IB_base, 
        const InterpolationBufferBase<Scalar>& IB_in, 
        const TimeRange<Scalar>& range
        ) const
{
  TEST_FOR_EXCEPTION((range.lower() != IB_base.getTimeRange().upper()) || (range.upper() != IB_base.getTimeRange().lower()),
      std::logic_error, 
      "Error, import range is not an append nor a prepend of the base interpolation buffer.\n"
      );
  TEST_FOR_EXCEPTION(range.lower() < IB_in.getTimeRange().lower(),
      std::logic_error,
      "Error, import range's lower bound does not sit inside incoming interpolation buffer's time range.\n"
      );
  TEST_FOR_EXCEPTION(IB_in.getTimeRange().upper() < range.upper(),
      std::logic_error,
      "Error, import range's upper bound does not sit inside incoming interpolation buffer's time range.\n"
      );
}


template<class Scalar>
class InterpolationBufferAppenderDefault : virtual public InterpolationBufferAppenderBase<Scalar>
{
  public:

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;

    /** \brief Concrete implementation that simply copies the nodal points
     * between the interpolation buffers.
     */
    void import(
        InterpolationBufferBase<Scalar>* IB_base, 
        const InterpolationBufferBase<Scalar>& IB_in, 
        const TimeRange<Scalar>& range
        );
    
    /// Redefined from Teuchos::Describable
    /** \brief . */
    std::string description() const;

    /** \brief . */
    void describe(
      Teuchos::FancyOStream       &out
      ,const Teuchos::EVerbosityLevel      verbLevel
      ) const;

    /// Redefined from Teuchos::ParameterListAcceptor
    /** \brief . */
    void setParameterList(RCP<Teuchos::ParameterList> const& paramList);

    /** \brief . */
    RCP<Teuchos::ParameterList> getParameterList();

    /** \brief . */
    RCP<Teuchos::ParameterList> unsetParameterList();

  private:
    RCP<Teuchos::ParameterList> parameterList_;

};

template<class Scalar>
void InterpolationBufferAppenderDefault<Scalar>::import(
        InterpolationBufferBase<Scalar>* IB_base, 
        const InterpolationBufferBase<Scalar>& IB_in, 
        const TimeRange<Scalar>& range
    ) 
{
#ifdef TEUCHOS_DEBUG
  InterpolationBufferAppenderBase<Scalar>::assertImportPreconditions(*IB_base,IB_in,range);
#endif // TEUCHOS_DEBUG

  Array<Scalar> time_vec;
  Array<RCP<const Thyra::VectorBase<Scalar> > > x_vec;
  Array<RCP<const Thyra::VectorBase<Scalar> > > xdot_vec;
  Array<ScalarMag> accuracy_vec;

  Array<Scalar> time_vec_in;
  IB_in.getNodes(&time_vec_in);

  selectPointsInTimeRange(&time_vec,time_vec_in,range);

  IB_in.getPoints(time_vec, &x_vec, &xdot_vec, &accuracy_vec);
  IB_base->addPoints(time_vec, x_vec, xdot_vec);
}


template<class Scalar>
std::string InterpolationBufferAppenderDefault<Scalar>::description() const
{
  std::string name = "Rythmos::InterpolationBufferAppenderDefault";
  return(name);
}

template<class Scalar>
void InterpolationBufferAppenderDefault<Scalar>::describe(
      Teuchos::FancyOStream       &out
      ,const Teuchos::EVerbosityLevel      verbLevel
      ) const
{
  if ( (Teuchos::as<int>(verbLevel) == Teuchos::as<int>(Teuchos::VERB_DEFAULT) ) ||
       (Teuchos::as<int>(verbLevel) >= Teuchos::as<int>(Teuchos::VERB_LOW)     )
     ) {
    out << description() << "::describe" << std::endl;
  }
}

template<class Scalar>
void InterpolationBufferAppenderDefault<Scalar>::setParameterList(RCP<Teuchos::ParameterList> const& paramList)
{
  TEST_FOR_EXCEPTION(paramList==Teuchos::null,std::logic_error,"Error, paramList == Teuchos::null!\n");
  parameterList_ = paramList;
}

template<class Scalar>
RCP<Teuchos::ParameterList> InterpolationBufferAppenderDefault<Scalar>::getParameterList()
{
  return(parameterList_);
}

template<class Scalar>
RCP<Teuchos::ParameterList> InterpolationBufferAppenderDefault<Scalar>::unsetParameterList()
{
  RCP<Teuchos::ParameterList> temp_param_list = parameterList_;
  parameterList_ = Teuchos::null;
  return(temp_param_list);
}

template<class Scalar>
class InterpolationBufferAppenderSmart : virtual public InterpolationBufferAppenderBase<Scalar>
{
  public:
    /** \brief Concrete implementation that attempts to use the order of
     * interpolation between the two interpolation buffers to be a bit smarter
     * about copying data between them.
     */
    void import(
        InterpolationBufferBase<Scalar>* IB_base, 
        const InterpolationBufferBase<Scalar>& IB_in, 
        const TimeRange<Scalar>& range
        );
};

template<class Scalar>
void InterpolationBufferAppenderSmart<Scalar>::import(
        InterpolationBufferBase<Scalar>* IB_base, 
        const InterpolationBufferBase<Scalar>& IB_in, 
        const TimeRange<Scalar>& range
    ) 
{
#ifdef TEUCHOS_DEBUG
  InterpolationBufferAppenderBase<Scalar>::assertImportPreconditions(*IB_base,IB_in,range);
#endif // TEUCHOS_DEBUG
  if (IB_base->getOrder() >= IB_in.getOrder()) {
    // The incoming interpolation buffer's order of interpolation is lower than
    // the base interpolation buffer's order of interpolation.  In this case,
    // we just copy the data over.
    InterpolationBufferAppenderDefault<Scalar> defaultAppender;
    defaultAppender.import(IB_base,IB_in,range);
  } else {
    // In this case, the incoming interpolation buffer's order of interpolation
    // is higher than the base interpolation buffer's, so we'll ask it to
    // interpolate points before inserting into the base interpolation buffer.
    TEST_FOR_EXCEPTION(
        true,std::logic_error,
        "Error, the smart interpolation buffer appender is not implemented for importing interpolation buffers with higher order interpolation into interpolation buffers with a lower order of interpolation yet.\n"
        );
    // 08/09/07 tscoffe:  Note:  you can't use selectPointsInTimeRange
    // immediately because you may need to interpolate points before the first
    // node inside the range.  I.e. IB_in.getNodes = [... , range.lower(), ... , range.upper(), ... ]
  }
}

/* This is the function from InterpolationBuffer.hpp which includes some of the
 * import functionality I want to include above.
 
template<class Scalar>
bool InterpolationBuffer<Scalar>::setRange(
  const TimeRange<Scalar>& range,
  const InterpolationBufferBase<Scalar>& IB
  )
{
  const Scalar time_lower = range.lower();
  const Scalar time_upper = range.upper();
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"IB::setRange");
  if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) ) {
    *out << "time_lower = " << time_lower << std::endl;
    *out << "time_upper = " << time_upper << std::endl;
    *out << "IB = " << IB.description() << std::endl;
  }
  Array<Scalar> input_nodes;
  bool status = IB.getNodes(&input_nodes);
  if (!status) { 
    return(status);
  }
  std::sort(input_nodes.begin(),input_nodes.end());
  if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) ) {
    *out << "input_nodes after sorting = " << std::endl;
    for (unsigned int i=0 ; i<input_nodes.size() ; ++i) {
      *out << "input_nodes[" << i << "] = " << input_nodes[i] << std::endl;
    }
  }
  // Remove nodes outside the range [time_lower,time_upper]
  typename Array<Scalar>::iterator input_it_lower = input_nodes.begin();
  for (; input_it_lower != input_nodes.end() ; input_it_lower++) {
    if (*input_it_lower >= time_lower) {
      break;
    }
  }
  if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) ) {
    int n0 = 0;
    int n1 = input_it_lower - input_nodes.begin();
    *out << "Removing input_nodes before time_lower with indices: [" << n0 << "," << n1 << ")" << std::endl;
    for (int i=n0 ; i<n1; ++i) {
      *out << "  input_nodes[" << i << "] = " << input_nodes[i] << std::endl;
    }
  }
  // tscoffe 10/19/06 Note:  erase removes the range [it_begin,it_end)
  if (input_it_lower - input_nodes.begin() >= 0) {
    input_nodes.erase(input_nodes.begin(),input_it_lower);
  }
  typename Array<Scalar>::iterator input_it_upper = input_nodes.end();
  input_it_upper--;
  for (; input_it_upper != input_nodes.begin() ; input_it_upper--) {
    if (*input_it_upper <= time_upper) {
      input_it_upper++;
      break;
    }
  }
  // This is to handle the case of one element in the vector
  if (*input_it_upper <= time_upper) {
    input_it_upper++;
  }
  if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) ) {
    int n0 = input_it_upper - input_nodes.begin();
    int n1 = input_nodes.size();
    *out << "Removing input_nodes after time_upper with indices [" << n0 << "," << n1 << ")" << std::endl;
    for (int i=n0 ; i<n1; ++i) {
      *out << "  input_nodes[" << i << "] = " << input_nodes[i] << std::endl;
    }
  }
  if (static_cast<unsigned int>(input_it_upper - input_nodes.begin()) < input_nodes.size()) {
    input_nodes.erase(input_it_upper,input_nodes.end());
  }
  if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) ) {
    *out << "input_nodes remaining:" << std::endl;
    for (unsigned int i=0 ; i<input_nodes.size() ; ++i) {
      *out << "input_nodes[" << i << "] = " << input_nodes[i] << std::endl;
    }
  }

  // Ask IB to interpolate more points if IB's order is higher than ours
  typedef Teuchos::ScalarTraits<Scalar> ST;
  Scalar h_safety = Scalar(2*ST::one());
  int IBOrder = IB.getOrder();
  if (IBOrder >= interpolator->order()) {
    std::list<Scalar> add_nodes;
    for (unsigned int i=0 ; i<input_nodes.size()-1 ; ++i) {
      Scalar h_0 = input_nodes[i+1] - input_nodes[i];
      Scalar h = pow(h_0,(IBOrder/interpolator->order())/h_safety);
      if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) ) {
        *out << "i = " << i << std::endl;
        *out << "interpolator->order() = " << interpolator->order() << std::endl;
        *out << "IB.getOrder() = " << IB.getOrder() << std::endl;
        *out << "h = " << h << std::endl;
      }
      Scalar N = ceil(h_0/h);
      h = Scalar(h_0/N);
      if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) ) {
        *out << "h_0 = " << h_0 << std::endl;
        *out << "N = " << N << std::endl;
        *out << "h = " << h << std::endl;
        *out << "Inserting an additional " << N-1 << " points to be interpolated:" << std::endl;
      }
      for (int j=1 ; j<N ; ++j) {
        if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) ) {
          *out << input_nodes[i]+j*h << std::endl;
        }
        add_nodes.push_back(input_nodes[i]+j*h);
      }
    }
    input_nodes.insert(input_nodes.end(),add_nodes.begin(),add_nodes.end());
    std::sort(input_nodes.begin(),input_nodes.end());
  }
  // If IB's order is lower than ours, then simply grab the node values and continue.
  // If IB's order is higher than ours, then grab the node values and ask IB to
  // interpolate extra values so that our order of accuracy is approximately
  // the same as the other IB's.
  // One approach:
  // Lets say IB's order is p and our order is r (p>r).
  // Given a particular interval with spacing h_0, the order of accuracy of IB is h_0^p
  // We want to find a spacing h such that h^r = h_0^p.  Clearly, this is h = h_0^{p/r}.
  // Given this new spacing, divide up the interval h_0 into h_0/h subintervals
  // and ask for the IB to interpolate points.  This will match basic order of
  // accuracy estimates.  Its probably a good idea to include a fudge factor in
  // there too.  E.g. h = h_0^{p/r}/fudge.

  // Don't forget to check the interval [time_lower,time_upper].
  // Use setPoints and check return value to make sure we observe storage_limit.

  Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > > input_x;
  Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > > input_xdot;
  Array<ScalarMag> input_accuracy;
  status = IB.getPoints( input_nodes, &input_x, &input_xdot, &input_accuracy );
  if (!status) { 
    return(status);
  }
  // We could check that the accuracy meets our criteria here.
  status = setPoints( input_nodes, input_x, input_xdot, input_accuracy );
  return(status);
}

*/



} // namespace Rythmos


#endif //Rythmos_INTERPOLATION_BUFFER_APPENDER_H
