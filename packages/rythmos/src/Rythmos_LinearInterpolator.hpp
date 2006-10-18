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

#include "Rythmos_Interpolator.hpp"

namespace Rythmos {

template<class Scalar>
class LinearInterpolator : virtual public Interpolator<Scalar>
{
  public:

    /// Destructor
    ~LinearInterpolator() {};

    /// Constructor
    LinearInterpolator();

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

  private:

#ifdef Rythmos_DEBUG
    int debugLevel;
    Teuchos::RefCountPtr<Teuchos::FancyOStream> debug_out;
#endif // Rythmos_DEBUG

};

template<class Scalar>
LinearInterpolator<Scalar>::LinearInterpolator()
{
#ifdef Rythmos_DEBUG
  debugLevel = 2;
  debug_out = Teuchos::VerboseObjectBase::getDefaultOStream();
  debug_out->precision(15);
  debug_out->setMaxLenLinePrefix(28);
  debug_out->pushLinePrefix("Rythmos::LinearInterpolator");
  debug_out->setShowLinePrefix(true);
  debug_out->setTabIndentStr("    ");
#endif // Rythmos_DEBUG
}

template<class Scalar>
bool LinearInterpolator<Scalar>::interpolate(
    const typename DataStore<Scalar>::DataStoreVector_t &data_in
    ,const std::vector<Scalar> &t_values
    ,typename DataStore<Scalar>::DataStoreVector_t *data_out ) const
{
#ifdef Rythmos_DEBUG
  Teuchos::OSTab ostab(debug_out,1,"interpolate");
  if (debugLevel > 1)
  {
    if (data_in.size() == 0)
      *debug_out << "data_in = empty vector" << std::endl;
    else
    {
      *debug_out << "data_in:" << std::endl;
      for (int i=0 ; i<data_in.size() ; ++i)
      {
        *debug_out << "data_in[" << i << "] = " << std::endl;
        data_in[i].describe(*debug_out,Teuchos::VERB_EXTREME);
      }
    }
    if (t_values.size() == 0)
      *debug_out << "t_values = empty vector" << std::endl;
    else
    {
      *debug_out << "t_values = " << std::endl;
      for (int i=0 ; i<t_values.size() ; ++i)
      {
        *debug_out << "t_values[" << i << "] = " << t_values[i] << std::endl;
      }
    }
    if (data_out == NULL)
      *debug_out << "data_out = NULL" << std::endl;
    else if (data_out->size() == 0)
      *debug_out << "data_out = empty vector" << std::endl;
    else
    {
      for (int i=0; i<data_out->size() ; ++i)
      {
        *debug_out << "data_out[" << i << "] = " << std::endl;
        (*data_out)[i].describe(*debug_out,Teuchos::VERB_EXTREME);
      }
    }
  }
#endif // Rythmos_DEBUG
  typedef Teuchos::ScalarTraits<Scalar> ST;
  // If there are fewer than 2 points in data_in, then return failure
  if (data_in.size() < 2)
    return(false);
  // Sort data_in: 
  typename DataStore<Scalar>::DataStoreVector_t local_data_in = data_in;
  std::sort(local_data_in.begin(),local_data_in.end());
  // Sort t_values:
  std::vector<Scalar> local_time_vec = t_values;
  std::sort(local_time_vec.begin(),local_time_vec.end());
  int N_t = local_time_vec.size();
  int N_D = local_data_in.size();
  // If time is outside range of t_values, then return failure
  if ( (local_time_vec[0] < local_data_in[0].time) || (local_time_vec[N_t] > local_data_in[N_D].time) )
    return(false);
  // Clear the output:
  data_out->clear();
  // Find t values on either side of time
  int j=0;
  for (int i=0 ; i<N_D-1 ; ++i)
  {
    while ((local_data_in[i].time <= local_time_vec[j]) && (local_time_vec[j] <= local_data_in[i+1].time ))
    {
      Scalar& t = local_time_vec[j];
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
      // Then we work on xdot.
      V_StVpStV(&*tmp_vec,Scalar(ST::one()/h),*xdoti,Scalar(-ST::one()/h),*xdotip1);
      Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > xdot = xdoti->clone_v();
      V_StVpStV(&*xdot, ST::one(), *xdoti, t-ti, *tmp_vec);
      DS.xdot = xdot;
      // And finally we estimate our order of accuracy
      DS.accuracy = h;
      // Push DataStore object onto vector:
      data_out->push_back(DS);
      // Increment local_time_vec time value
      j++;
    }
  }
  return(true);
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
  if (verbLevel == Teuchos::VERB_EXTREME)
  {
    out << description() << "::describe" << std::endl;
  }
}

} // namespace Rythmos

#endif // Rythmos_LINEAR_INTERPOLATOR_H
