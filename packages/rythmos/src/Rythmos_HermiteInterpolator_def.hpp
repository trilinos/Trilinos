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

#ifndef Rythmos_HERMITE_INTERPOLATOR_DEF_H
#define Rythmos_HERMITE_INTERPOLATOR_DEF_H

#include "Rythmos_HermiteInterpolator_decl.hpp"
#include "Rythmos_InterpolatorBaseHelpers.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"

namespace Rythmos {

template<class Scalar>
HermiteInterpolator<Scalar>::HermiteInterpolator()
{
  nodes_ = Teuchos::null;
  parameterList_ = Teuchos::null;
}

template<class Scalar>
void HermiteInterpolator<Scalar>::setNodes(
  const RCP<const typename DataStore<Scalar>::DataStoreVector_t> & nodes
  )
{
  nodes_ = nodes;
}

template<class Scalar>
void HermiteInterpolator<Scalar>::interpolate(
    const Array<Scalar> &t_values
    ,typename DataStore<Scalar>::DataStoreVector_t *data_out ) const
{

  //TEST_FOR_EXCEPT_MSG(true, "Error, ths function is not tested!" );

  typedef Teuchos::ScalarTraits<Scalar> ST;
#ifdef RYTHMOS_DEBUG
  assertInterpolatePreconditions((*nodes_),t_values,data_out);
#endif // RYTHMOS_DEBUG
  RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"HI::interpolator");
  if ( Teuchos::as<int>(this->getVerbLevel()) >= Teuchos::as<int>(Teuchos::VERB_HIGH) ) {
    *out << "(*nodes_):" << std::endl;
    for (Teuchos::Ordinal i=0 ; i<(*nodes_).size() ; ++i) {
      *out << "(*nodes_)[" << i << "] = " << std::endl;
      (*nodes_)[i].describe(*out,Teuchos::VERB_EXTREME);
    }
    *out << "t_values = " << std::endl;
    for (Teuchos::Ordinal i=0 ; i<t_values.size() ; ++i) {
      *out << "t_values[" << i << "] = " << t_values[i] << std::endl;
    }
    for (Teuchos::Ordinal i=0; i<data_out->size() ; ++i) {
      *out << "data_out[" << i << "] = " << std::endl;
      (*data_out)[i].describe(*out,Teuchos::VERB_EXTREME);
    }
  }
  data_out->clear();
  if (t_values.size() == 0) {
    return;
  }
  
  if ((*nodes_).size() == 1) {
    // trivial case of one node
    // preconditions assert that t_values[0] == (*nodes_)[0].time so we can just pass it out
    DataStore<Scalar> DS((*nodes_)[0]);
    data_out->push_back(DS);
  } else {
    // (*nodes_).size() >= 2
    int n = 0;
    for (int i=0 ; i<Teuchos::as<int>((*nodes_).size())-1 ; ++i) {
      const Scalar& t0 = (*nodes_)[i].time;
      const Scalar& t1 = (*nodes_)[i+1].time;
      while ((t0 <= t_values[n]) && (t_values[n] <= t1)) {
        const Scalar& t = t_values[n];
        // First we check for exact node matches:
        if (t == t0) {
          DataStore<Scalar> DS((*nodes_)[i]);
          data_out->push_back(DS);
        } else if (t == t1) {
          DataStore<Scalar> DS((*nodes_)[i+1]);
          data_out->push_back(DS);
        } else {
          RCP<const Thyra::VectorBase<Scalar> > x0    = (*nodes_)[i  ].x;
          RCP<const Thyra::VectorBase<Scalar> > x1    = (*nodes_)[i+1].x;
          RCP<const Thyra::VectorBase<Scalar> > xdot0 = (*nodes_)[i  ].xdot;
          RCP<const Thyra::VectorBase<Scalar> > xdot1 = (*nodes_)[i+1].xdot;
          
          // 10/10/06 tscoffe:  this could be expensive:
          RCP<Thyra::VectorBase<Scalar> > tmp_vec = x0->clone_v(); 
          RCP<Thyra::VectorBase<Scalar> > xdot_temp = x1->clone_v(); 
          Scalar dt = t1-t0;
          Scalar dt2 = dt*dt;
          Scalar t_t0 = t - t0;
          Scalar t_t1 = t - t1;
          Scalar tmp_t;

          // Compute numerical divided difference:
          Thyra::Vt_S(xdot_temp.ptr(),Scalar(ST::one()/dt));
          Thyra::Vp_StV(xdot_temp.ptr(),Scalar(-ST::one()/dt),*x0);

          // interpolate this point
          DataStore<Scalar> DS;
          DS.time = t;

          //  H_3(t) = x(t0) + xdot(t0)(t-t0) + ((x(t1)-x(t0))/(t1-t0) - xdot(t0))(t-t0)^2/(t1-t0)
          //           +(xdot(t1) - 2(x(t1)-x(t0))/(t1-t0) + xdot(t0))(t-t0)^2(t-t1)/(t1-t0)^2
          RCP<Thyra::VectorBase<Scalar> > x_vec = x0->clone_v(); 
          Thyra::Vp_StV(x_vec.ptr(),t_t0,*xdot0);
          tmp_t = t_t0*t_t0/dt;
          Thyra::V_StVpStV(tmp_vec.ptr(),tmp_t,*xdot_temp,Scalar(-ST::one()*tmp_t),*xdot0);
          Thyra::Vp_V(x_vec.ptr(),*tmp_vec);
          tmp_t = t_t0*t_t0*t_t1/dt2;
          Thyra::V_StVpStV(tmp_vec.ptr(),tmp_t,*xdot1,Scalar(-2*tmp_t),*xdot_temp);
          Thyra::Vp_StV(tmp_vec.ptr(),tmp_t,*xdot0);
          Thyra::Vp_V(x_vec.ptr(),*tmp_vec);
          DS.x = x_vec;

          //  H_3'(t) =        xdot(t0) + 2*((x(t1)-x(t0))/(t1-t0) - xdot(t0))(t-t0)/(t1-t0)
          //           +(xdot(t1) - 2(x(t1)-x(t0))/(t1-t0) + xdot(t0))[2*(t-t0)(t-t1) + (t-t0)^2]/(t1-t0)^2
          RCP<Thyra::VectorBase<Scalar> > xdot_vec = xdot0->clone_v(); 
          tmp_t = t_t0/dt;
          Thyra::Vp_StV(xdot_vec.ptr(),Scalar(2*tmp_t),*xdot_temp);
          Thyra::Vp_StV(xdot_vec.ptr(),Scalar(-2*tmp_t),*xdot0);
          tmp_t = Scalar((2*t_t0*t_t1+t_t0*t_t0)/dt2);
          Thyra::V_StVpStV(tmp_vec.ptr(),tmp_t,*xdot1,Scalar(-2*tmp_t),*xdot_temp);
          Thyra::Vp_StV(tmp_vec.ptr(),tmp_t,*xdot0);
          Thyra::Vp_V(xdot_vec.ptr(),*tmp_vec);
          DS.xdot = xdot_vec;
          
          // Accuracy:
          // f(x) - H_3(x) = (f^{(3)}(\xi(x))/(4!))(x-x0)^2(x-x1)^2
          DS.accuracy = (t_t0)*(t_t0)*(t_t1)*(t_t1);

          // Push DataStore object onto vector:
          data_out->push_back(DS);
        }
        n++;
        if (n == Teuchos::as<int>(t_values.size())) {
          return;
        }
      }
    }
  } // (*nodes_).size() == 1
}

// non-member constructor
template<class Scalar>
RCP<HermiteInterpolator<Scalar> > hermiteInterpolator()
{
  RCP<HermiteInterpolator<Scalar> > hi = rcp(new HermiteInterpolator<Scalar>() );
  return hi;
}

template<class Scalar>
int HermiteInterpolator<Scalar>::order() const
{
  return(3);
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
  if ( (Teuchos::as<int>(verbLevel) == Teuchos::as<int>(Teuchos::VERB_DEFAULT) ) ||
       (Teuchos::as<int>(verbLevel) >= Teuchos::as<int>(Teuchos::VERB_LOW)     )
     )
  {
    out << description() << "::describe" << std::endl;
  }
  else if (Teuchos::as<int>(verbLevel) >= Teuchos::as<int>(Teuchos::VERB_LOW))
  {}
  else if (Teuchos::as<int>(verbLevel) >= Teuchos::as<int>(Teuchos::VERB_MEDIUM))
  {}
  else if (Teuchos::as<int>(verbLevel) >= Teuchos::as<int>(Teuchos::VERB_HIGH))
  {}
}

template <class Scalar>
void HermiteInterpolator<Scalar>::setParameterList(RCP<ParameterList> const& paramList)
{
  TEST_FOR_EXCEPT(is_null(paramList));
  paramList->validateParametersAndSetDefaults(*this->getValidParameters());
  parameterList_ = paramList;
  Teuchos::readVerboseObjectSublist(&*parameterList_,this);
}

template <class Scalar>
RCP<ParameterList> HermiteInterpolator<Scalar>::getNonconstParameterList()
{
  return(parameterList_);
}

template <class Scalar>
RCP<ParameterList> HermiteInterpolator<Scalar>::unsetParameterList()
{
  RCP<ParameterList> temp_param_list = parameterList_;
  parameterList_ = Teuchos::null;
  return(temp_param_list);
}

template<class Scalar>
RCP<const Teuchos::ParameterList> HermiteInterpolator<Scalar>::getValidParameters() const
{
  static RCP<Teuchos::ParameterList> validPL;
  if (is_null(validPL)) {
    RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    Teuchos::setupVerboseObjectSublist(&*pl);
    validPL = pl;
  }
  return (validPL);
}

template<class Scalar>
void HermiteInterpolator<Scalar>::assertInterpolatePreconditions(
        const typename DataStore<Scalar>::DataStoreVector_t &data_in
        ,const Array<Scalar> &t_values
        ,typename DataStore<Scalar>::DataStoreVector_t *data_out
        ) const
{
  assertBaseInterpolatePreconditions(data_in,t_values,data_out);
  for (int i=0; i<Teuchos::as<int>(data_in.size()) ; ++i) {
    TEST_FOR_EXCEPTION(
        is_null(data_in[i].x), std::logic_error,
        "Error, data_in[" << i << "].x == Teuchos::null.\n"
        );
    TEST_FOR_EXCEPTION(
        is_null(data_in[i].xdot), std::logic_error,
        "Error, data_in[" << i << "].xdot == Teuchos::null.\n"
        );
  }
}

// 
// Explicit Instantiation macro
//
// Must be expanded from within the Rythmos namespace!
//

#define RYTHMOS_HERMITE_INTERPOLATOR_INSTANT(SCALAR) \
  \
  template class HermiteInterpolator< SCALAR >; \
  \
  template RCP<HermiteInterpolator< SCALAR > > hermiteInterpolator(); 


} // namespace Rythmos

#endif // Rythmos_HERMITE_INTERPOLATOR_DEF_H
