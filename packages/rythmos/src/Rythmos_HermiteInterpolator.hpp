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

#ifndef Rythmos_HERMITE_INTERPOLATOR_H
#define Rythmos_HERMITE_INTERPOLATOR_H

#include "Rythmos_InterpolatorBase.hpp"

namespace Rythmos {

template<class Scalar>
class HermiteInterpolator : virtual public InterpolatorBase<Scalar>
{
  public:

    /// Destructor
    ~HermiteInterpolator() {};

    /// Constructor
    HermiteInterpolator();

    /// Interpolation:
    bool interpolate(
        const typename DataStore<Scalar>::DataStoreVector_t &data_in
        ,const std::vector<Scalar> &t_values
        ,typename DataStore<Scalar>::DataStoreVector_t *data_out
        ) const;

    /// Order of interpolation:
    int order() const; 

    /// Inherited from Teuchos::Describable
    /** \brief . */
    std::string description() const;

    /** \brief . */
    void describe(
      Teuchos::FancyOStream       &out
      ,const Teuchos::EVerbosityLevel      verbLevel
      ) const;

    /// Redefined from Teuchos::ParameterListAcceptor
    /** \brief . */
    void setParameterList(Teuchos::RefCountPtr<Teuchos::ParameterList> const& paramList);

    /** \brief . */
    Teuchos::RefCountPtr<Teuchos::ParameterList> getParameterList();

    /** \brief . */
    Teuchos::RefCountPtr<Teuchos::ParameterList> unsetParameterList();

  private:

    Teuchos::RefCountPtr<Teuchos::ParameterList> parameterList;

};

template<class Scalar>
HermiteInterpolator<Scalar>::HermiteInterpolator()
{
  Teuchos::RefCountPtr<Teuchos::FancyOStream> out = this->getOStream();
  out->precision(15);
  out->setMaxLenLinePrefix(28);
  out->pushLinePrefix("Rythmos::HermiteInterpolator");
  out->setShowLinePrefix(true);
  out->setTabIndentStr("    ");
}

template<class Scalar>
bool HermiteInterpolator<Scalar>::interpolate(
    const typename DataStore<Scalar>::DataStoreVector_t &data_in
    ,const std::vector<Scalar> &t_values
    ,typename DataStore<Scalar>::DataStoreVector_t *data_out ) const
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  Teuchos::RefCountPtr<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"HI::interpolator");
  if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
  {
    if (data_in.size() == 0)
      *out << "data_in = empty vector" << std::endl;
    else
    {
      *out << "data_in:" << std::endl;
      for (int i=0 ; i<data_in.size() ; ++i)
      {
        *out << "data_in[" << i << "] = " << std::endl;
        data_in[i].describe(*out,Teuchos::VERB_EXTREME);
      }
    }
    if (t_values.size() == 0)
      *out << "t_values = empty vector" << std::endl;
    else
    {
      *out << "t_values = " << std::endl;
      for (int i=0 ; i<t_values.size() ; ++i)
      {
        *out << "t_values[" << i << "] = " << t_values[i] << std::endl;
      }
    }
    if (data_out == NULL)
      *out << "data_out = NULL" << std::endl;
    else if (data_out->size() == 0)
      *out << "data_out = empty vector" << std::endl;
    else
    {
      for (int i=0; i<data_out->size() ; ++i)
      {
        *out << "data_out[" << i << "] = " << std::endl;
        (*data_out)[i].describe(*out,Teuchos::VERB_EXTREME);
      }
    }
  }
  // Sort data_in: 
  typename DataStore<Scalar>::DataStoreVector_t local_data_in;
  local_data_in.insert(local_data_in.end(),data_in.begin(),data_in.end());
  std::sort(local_data_in.begin(),local_data_in.end());
  // Clear the output:
  data_out->clear();
  // Sort incoming t values:
  std::list<Scalar> local_t_values;
  local_t_values.insert(local_t_values.end(),t_values.begin(),t_values.end());
  local_t_values.sort();
  // Check for node values to pass out directly:
  typename std::list<Scalar>::iterator time_it = local_t_values.begin();
  while (time_it != local_t_values.end())
  {
    typename DataStore<Scalar>::DataStoreVector_t::iterator data_in_it;
    data_in_it = std::find(local_data_in.begin(),local_data_in.end(),*time_it);
    if (data_in_it != local_data_in.end())
    {
      if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
        *out << "Passing out data (w/o interpolating) for t = " << *time_it << std::endl;
      data_out->push_back(*data_in_it);
      time_it = local_t_values.erase(time_it);
    }
    else
    {
      time_it++;
    }
  }
  if (local_t_values.size() == 0)
    return(true);
  // If there are fewer than 2 points in data_in, then return failure
  if (local_data_in.size() < 2)
    return(false);

  // If local_t_values are outside range of data_in, then return failure
  if ( ! ( (local_data_in.front() <= local_t_values.front()) && 
           (local_data_in.back()  >= local_t_values.back() ) )    )
    return(false);
  // Find t values on either side of time
  time_it = local_t_values.begin();
  for (int i=0 ; i<local_data_in.size()-1 ; ++i)
  {
    while ( (local_data_in[i] <= *time_it) && (local_data_in[i+1] >= *time_it) )
    {
      Scalar& t = *time_it;
      Scalar& ti = local_data_in[i].time;
      Scalar& tip1 = local_data_in[i+1].time;
      Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > xi = local_data_in[i].x;
      Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > xip1 = local_data_in[i+1].x;
      Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > xdoti = local_data_in[i].xdot;
      Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > xdotip1 = local_data_in[i+1].xdot;

      // 10/10/06 tscoffe:  this could be expensive:
      Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > tmp_vec = xi->clone_v(); 

      // interpolate this point
      DataStore<Scalar> DS;
      DS.time = t;
      Scalar h = tip1-ti;
      // First we work on x.
      V_StVpStV(&*tmp_vec,Scalar(ST::one()/h),*xi,Scalar(-ST::one()/h),*xip1);
      Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > x = xi->clone_v();
      V_StVpStV(&*x, ST::one(), *xi, t-ti, *tmp_vec);
      DS.x = x;
      // Check that xdot != Teuchos::null
      if ( (xdoti != Teuchos::null) && (xdotip1 != Teuchos::null) )
      {
        // Then we work on xdot.
        V_StVpStV(&*tmp_vec,Scalar(ST::one()/h),*xdoti,Scalar(-ST::one()/h),*xdotip1);
        Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > xdot = xdoti->clone_v();
        V_StVpStV(&*xdot, ST::one(), *xdoti, t-ti, *tmp_vec);
        DS.xdot = xdot;
      }
      else
        DS.xdot = xdoti;
      // And finally we estimate our order of accuracy
      DS.accuracy = h;
      // Push DataStore object onto vector:
      data_out->push_back(DS);
      // Increment time_it:
      time_it++;
    }
  }
  return(true);
}

template<class Scalar>
int HermiteInterpolator<Scalar>::order() const
{
  return(1);
}

template<class Scalar>
std::string HermiteInterpolator<Scalar>::description() const
{
  std::string name = "Rythmos::HermiteInterpolator";
  return(name);
}

template<class Scalar>
void HermiteInterpolator<Scalar>::describe(
      Teuchos::FancyOStream       &out
      ,const Teuchos::EVerbosityLevel      verbLevel
      ) const
{
  if ( (static_cast<int>(verbLevel) == static_cast<int>(Teuchos::VERB_DEFAULT) ) ||
       (static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW)     )
     )
  {
    out << description() << "::describe" << std::endl;
  }
  else if (static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
  {
  }
  else if (static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_MEDIUM))
  {
  }
  else if (static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_HIGH))
  {
  }
}

template <class Scalar>
void HermiteInterpolator<Scalar>::setParameterList(Teuchos::RefCountPtr<Teuchos::ParameterList> const& paramList)
{
  parameterList = paramList;
  int outputLevel = parameterList->get( "outputLevel", int(-1) );
  outputLevel = min(max(outputLevel,-1),4);
  this->setVerbLevel(static_cast<Teuchos::EVerbosityLevel>(outputLevel));
}

template <class Scalar>
Teuchos::RefCountPtr<Teuchos::ParameterList> HermiteInterpolator<Scalar>::getParameterList()
{
  return(parameterList);
}

template <class Scalar>
Teuchos::RefCountPtr<Teuchos::ParameterList> HermiteInterpolator<Scalar>::unsetParameterList()
{
  Teuchos::RefCountPtr<Teuchos::ParameterList> temp_param_list = parameterList;
  parameterList = Teuchos::null;
  return(temp_param_list);
}

} // namespace Rythmos

#endif // Rythmos_HERMITE_INTERPOLATOR_H
