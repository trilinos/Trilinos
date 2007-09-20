//@HEADER
// ***********************************************************************
//
//                           Rythmos Package
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

#ifndef Rythmos_LINEAR_INTERPOLATOR_H
#define Rythmos_LINEAR_INTERPOLATOR_H

#include "Rythmos_InterpolatorBase.hpp"

namespace Rythmos {

template<class Scalar>
class LinearInterpolator : virtual public InterpolatorBase<Scalar>
{
public:

  /// Destructor
  ~LinearInterpolator() {};

  /// Constructor
  LinearInterpolator();

  /** \brief . */
  bool supportsCloning() const;

  /** \brief . */
  Teuchos::RCP<InterpolatorBase<Scalar> > cloneInterpolator() const;

  /** \brief . */
  void interpolate(
    const typename DataStore<Scalar>::DataStoreVector_t &data_in
    ,const Array<Scalar> &t_values
    ,typename DataStore<Scalar>::DataStoreVector_t *data_out
    ) const;

  /** \brief . */
  int order() const; 

  /** \brief . */
  std::string description() const;

  /** \brief . */
  void describe(
    Teuchos::FancyOStream &out,
    const Teuchos::EVerbosityLevel verbLevel
    ) const;

  /** \brief . */
  void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& paramList);

  /** \brief . */
  Teuchos::RCP<Teuchos::ParameterList> getParameterList();

  /** \brief . */
  Teuchos::RCP<Teuchos::ParameterList> unsetParameterList();

private:

  Teuchos::RCP<Teuchos::ParameterList> parameterList_;

};

template<class Scalar>
LinearInterpolator<Scalar>::LinearInterpolator()
{}

template<class Scalar>
bool LinearInterpolator<Scalar>::supportsCloning() const
{
  return true;
}

template<class Scalar>
Teuchos::RCP<InterpolatorBase<Scalar> >
LinearInterpolator<Scalar>::cloneInterpolator() const
{
  Teuchos::RCP<LinearInterpolator<Scalar> >
    interpolator = Teuchos::rcp(new LinearInterpolator<Scalar>);
  if (!is_null(parameterList_))
    interpolator->parameterList_ = parameterList(*parameterList_);
  return interpolator;
}

template<class Scalar>
void LinearInterpolator<Scalar>::interpolate(
    const typename DataStore<Scalar>::DataStoreVector_t &data_in
    ,const Array<Scalar> &t_values
    ,typename DataStore<Scalar>::DataStoreVector_t *data_out ) const
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
#ifdef TEUCHOS_DEBUG
  assertBaseInterpolatePreconditions(data_in,t_values,data_out);
#endif // TEUCHOS_DEBUG
  data_out->clear();
  if (t_values.size() == 0) {
    return;
  }
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"LI::interpolator");
  if ( Teuchos::as<int>(this->getVerbLevel()) >= Teuchos::as<int>(Teuchos::VERB_HIGH) ) {
    *out << "data_in:" << std::endl;
    for (unsigned int i=0 ; i<data_in.size() ; ++i) {
      *out << "data_in[" << i << "] = " << std::endl;
      data_in[i].describe(*out,Teuchos::VERB_EXTREME);
    }
    *out << "t_values = " << std::endl;
    for (unsigned int i=0 ; i<t_values.size() ; ++i) {
      *out << "t_values[" << i << "] = " << t_values[i] << std::endl;
    }
  }

  if (data_in.size() == 1) {
    // trivial case of one node
    // preconditions assert that t_values[0] == data_in[0].time so we can just pass it out
    DataStore<Scalar> DS(data_in[0]);
    data_out->push_back(DS);
  } else {
    // data_in.size() >= 2
    int n = 0; // index into t_values
    for (int i=0 ; i<Teuchos::as<int>(data_in.size())-1 ; ++i) {
      const Scalar& ti = data_in[i].time;
      const Scalar& tip1 = data_in[i+1].time;
      const Scalar  h = tip1-ti;
      while ((ti <= t_values[n]) && (t_values[n] <= tip1)) {
        // First we check for exact node matches:
        if (t_values[n] == ti) {
          DataStore<Scalar> DS(data_in[i]);
          data_out->push_back(DS);
        } else if (t_values[n] == tip1) {
          DataStore<Scalar> DS(data_in[i+1]);
          data_out->push_back(DS);
        } else {
          const Scalar& t = t_values[n];
          Teuchos::RCP<const Thyra::VectorBase<Scalar> > xi =      data_in[i  ].x;
          Teuchos::RCP<const Thyra::VectorBase<Scalar> > xip1 =    data_in[i+1].x;
          Teuchos::RCP<const Thyra::VectorBase<Scalar> > xdoti =   data_in[i  ].xdot;
          Teuchos::RCP<const Thyra::VectorBase<Scalar> > xdotip1 = data_in[i+1].xdot;
        
          // interpolate this point
          //
          // x(t) = (t-ti)/(tip1-ti) * xip1 + (1-(t-ti)/(tip1-ti)) * xi
          //
          // Above, it is easy to see that:
          //
          //    x(ti) = xi
          //    x(tip1) = xip1
          //
          DataStore<Scalar> DS;
          DS.time = t;
          const Scalar dt = t-ti;
          const Scalar dt_over_h = dt / h;
          const Scalar one_minus_dt_over_h = ST::one() - dt_over_h;
          // x
          RCP<Thyra::VectorBase<Scalar> > x = createMember(xi->space());
          Thyra::V_StVpStV(&*x,dt_over_h,*xip1,one_minus_dt_over_h,*xi);
          DS.x = x;
          // xdot
          RCP<Thyra::VectorBase<Scalar> > xdot;
          if ((xdoti != Teuchos::null) && (xdotip1 != Teuchos::null)) {
            xdot = createMember(xdoti->space());
            Thyra::V_StVpStV(&*xdot,dt_over_h,*xdotip1,one_minus_dt_over_h,*xdoti);
          }
          DS.xdot = xdot;
          // And finally we estimate our order of accuracy
          DS.accuracy = h;
          // Push DataStore object onto vector:
          data_out->push_back(DS);
          // Move to the next time to consider!
        }
        n++;
        if (n == Teuchos::as<int>(t_values.size())) {
          return;
        }
      }
    }
  } // data_in.size() == 1
}

template<class Scalar>
int LinearInterpolator<Scalar>::order() const
{
  return(1);
}

template<class Scalar>
std::string LinearInterpolator<Scalar>::description() const
{
  std::string name = "Rythmos::LinearInterpolator";
  return(name);
}

template<class Scalar>
void LinearInterpolator<Scalar>::describe(
      Teuchos::FancyOStream       &out
      ,const Teuchos::EVerbosityLevel      verbLevel
      ) const
{
  if ( (Teuchos::as<int>(verbLevel) == Teuchos::as<int>(Teuchos::VERB_DEFAULT) ) ||
       (Teuchos::as<int>(verbLevel) >= Teuchos::as<int>(Teuchos::VERB_LOW)     )
     ) {
    out << description() << "::describe" << std::endl;
  } else if (Teuchos::as<int>(verbLevel) >= Teuchos::as<int>(Teuchos::VERB_LOW)) {
  } else if (Teuchos::as<int>(verbLevel) >= Teuchos::as<int>(Teuchos::VERB_MEDIUM)) {
  } else if (Teuchos::as<int>(verbLevel) >= Teuchos::as<int>(Teuchos::VERB_HIGH)) {
  }
}

template <class Scalar>
void LinearInterpolator<Scalar>::setParameterList(
  Teuchos::RCP<Teuchos::ParameterList> const& paramList
  )
{
  parameterList_ = paramList;
  int outputLevel = parameterList_->get( "outputLevel", int(-1) );
  outputLevel = std::min(std::max(outputLevel,-1),4);
  this->setVerbLevel(static_cast<Teuchos::EVerbosityLevel>(outputLevel));
  // 2007/05/18: rabartl: ToDo: Replace with standard "Verbose Object"
  // sublist! and validate the sublist!

}

template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList> LinearInterpolator<Scalar>::getParameterList()
{
  return(parameterList_);
}

template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList> LinearInterpolator<Scalar>::unsetParameterList()
{
  Teuchos::RCP<Teuchos::ParameterList> temp_param_list = parameterList_;
  parameterList_ = Teuchos::null;
  return(temp_param_list);
}

} // namespace Rythmos

#endif // Rythmos_LINEAR_INTERPOLATOR_H
