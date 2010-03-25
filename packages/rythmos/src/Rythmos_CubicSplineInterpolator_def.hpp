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

#ifndef Rythmos_CUBIC_SPLINE_INTERPOLATOR_DEF_H
#define Rythmos_CUBIC_SPLINE_INTERPOLATOR_DEF_H

#include "Rythmos_CubicSplineInterpolator_decl.hpp"
#include "Rythmos_InterpolatorBaseHelpers.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_VectorSpaceBase.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"

namespace Rythmos {

template<class Scalar>
Teuchos::RCP<Rythmos::CubicSplineInterpolator<Scalar> >
cubicSplineInterpolator()
{
  RCP<CubicSplineInterpolator<Scalar> > csi = Teuchos::rcp(new CubicSplineInterpolator<Scalar>() );
  return csi;
}

template<class Scalar>
void computeCubicSplineCoeff(
    const typename DataStore<Scalar>::DataStoreVector_t & data,
    const Ptr<CubicSplineCoeff<Scalar> > & coeffPtr
    )
{
  using Teuchos::outArg;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  using Thyra::createMember;
  TEST_FOR_EXCEPTION( 
      (data.size() < 2), std::logic_error,
      "Error!  A minimum of two data points is required for this cubic spline."
      );
  // time data in the DataStoreVector should be unique and sorted 
  Array<Scalar> t;
  Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > > x_vec, xdot_vec;
  dataStoreVectorToVector<Scalar>( data, &t, &x_vec, &xdot_vec, NULL );
#ifdef RYTHMOS_DEBUG
  assertTimePointsAreSorted<Scalar>( t );
#endif // RYTHMOS_DEBUG
  // 11/18/08 tscoffe:  Question:  Should I erase everything in coeffPtr or
  // re-use what I can?  For now, I'll erase and create new each time.
  CubicSplineCoeff<Scalar>& coeff = *coeffPtr;
  // If there are only two points, then we do something special and just create
  // a linear polynomial between the points and return.
  if (t.size() == 2) {
    coeff.t.clear(); 
    coeff.a.clear(); coeff.b.clear(); coeff.c.clear(); coeff.d.clear();
    coeff.t.reserve(2); 
    coeff.a.reserve(1); coeff.b.reserve(1); coeff.c.reserve(1); coeff.d.reserve(1);
    coeff.t.push_back(t[0]);
    coeff.t.push_back(t[1]);
    coeff.a.push_back(x_vec[0]->clone_v());
    coeff.b.push_back(createMember(x_vec[0]->space()));
    coeff.c.push_back(createMember(x_vec[0]->space()));
    coeff.d.push_back(createMember(x_vec[0]->space()));
    Scalar h = coeff.t[1] - coeff.t[0];
    V_StVpStV(outArg(*coeff.b[0]),ST::one()/h,*x_vec[1],-ST::one()/h,*x_vec[0]);
    V_S(outArg(*coeff.c[0]),ST::zero());
    V_S(outArg(*coeff.d[0]),ST::zero());
    return;
  }
  // Data objects we'll need:
  int n = t.length()-1; // Number of intervals
  coeff.t.clear(); coeff.t.reserve(n+1);
  coeff.a.clear(); coeff.a.reserve(n+1);
  coeff.b.clear(); coeff.b.reserve(n);
  coeff.c.clear(); coeff.c.reserve(n+1);
  coeff.d.clear(); coeff.d.reserve(n);

  Array<Scalar> h(n);
  Array<RCP<Thyra::VectorBase<Scalar> > > alpha(n), z(n+1);
  Array<Scalar> l(n+1), mu(n);
  for (int i=0 ; i<n ; ++i) {
    coeff.t.push_back(t[i]);
    coeff.a.push_back(x_vec[i]->clone_v());
    coeff.b.push_back(Thyra::createMember(x_vec[0]->space()));
    coeff.c.push_back(Thyra::createMember(x_vec[0]->space()));
    coeff.d.push_back(Thyra::createMember(x_vec[0]->space()));
    alpha[i] = Thyra::createMember(x_vec[0]->space());
    z[i] = Thyra::createMember(x_vec[0]->space());
  }
  coeff.a.push_back(x_vec[n]->clone_v());
  coeff.t.push_back(t[n]);
  coeff.c.push_back(Thyra::createMember(x_vec[0]->space()));
  z[n] = Thyra::createMember(x_vec[0]->space());
  Scalar zero = ST::zero();
  Scalar one = ST::one();
  Scalar two = Scalar(2*ST::one());
  Scalar three = Scalar(3*ST::one());

  // Algorithm starts here:
  for (int i=0 ; i<n ; ++i) {
    h[i] = coeff.t[i+1]-coeff.t[i];
  }
  for (int i=1 ; i<n ; ++i) {
    V_StVpStV(outArg(*(alpha[i])),three/h[i],*coeff.a[i+1],-3/h[i],*coeff.a[i]);
    Vp_StV(outArg(*(alpha[i])),-three/h[i-1],*coeff.a[i]);
    Vp_StV(outArg(*(alpha[i])),+three/h[i-1],*coeff.a[i-1]);
  }
  l[0] = one;
  mu[0] = zero;
  V_S(outArg(*(z[0])),zero);
  for (int i=1 ; i<n ; ++i) {
    l[i] = 2*(coeff.t[i+1]-coeff.t[i-1])-h[i-1]*mu[i-1];
    mu[i] = h[i]/l[i];
    V_StVpStV(outArg(*(z[i])),one/l[i],*alpha[i],-h[i-1]/l[i],*z[i-1]);
  }
  l[n] = one;
  V_S(outArg(*(z[n])),zero);
  V_S(outArg(*(coeff.c[n])),zero);
  for (int j=n-1 ; j >= 0 ; --j) {
    V_StVpStV(outArg(*(coeff.c[j])),one,*z[j],-mu[j],*coeff.c[j+1]);
    V_StVpStV(outArg(*(coeff.b[j])),one/h[j],*coeff.a[j+1],-one/h[j],*coeff.a[j]);
    Vp_StV(outArg(*(coeff.b[j])),-h[j]/three,*coeff.c[j+1]);
    Vp_StV(outArg(*(coeff.b[j])),-h[j]*two/three,*coeff.c[j]);
    V_StVpStV(outArg(*(coeff.d[j])),one/(three*h[j]),*coeff.c[j+1],-one/(three*h[j]),*coeff.c[j]);
  }
  // Pop the last entry off of a and c to make them the right size.
  coeff.a.erase(coeff.a.end()-1);
  coeff.c.erase(coeff.c.end()-1);
}


template<class Scalar>
void validateCubicSplineCoeff(const CubicSplineCoeff<Scalar>& coeff) 
{
  int t_n = coeff.t.size();
  int a_n = coeff.a.size();
  int b_n = coeff.b.size();
  int c_n = coeff.c.size();
  int d_n = coeff.d.size();
  TEST_FOR_EXCEPTION( 
      ((a_n != t_n-1) || (a_n != b_n) || (a_n != c_n) || (a_n != d_n)),
      std::logic_error,
      "Error!  The sizes of the data structures in the CubicSplineCoeff object do not match"
      );
}


template<class Scalar>
void evaluateCubicSpline(
    const CubicSplineCoeff<Scalar>& coeff,
    Teuchos::Ordinal j, 
    const Scalar& t,
    const Ptr<Thyra::VectorBase<Scalar> >& S,
    const Ptr<Thyra::VectorBase<Scalar> >& Sp, 
    const Ptr<Thyra::VectorBase<Scalar> >& Spp
    )
{
  using Teuchos::outArg;
  using Teuchos::as;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  // Assert preconditions:
  validateCubicSplineCoeff<Scalar>(coeff);
  TEST_FOR_EXCEPTION( as<Teuchos::Ordinal>(j) >= coeff.a.size(),
     std::out_of_range, "Error!, j is out of range" );

  Scalar dt = t-coeff.t[j];
  const Thyra::VectorBase<Scalar>& a = *(coeff.a[j]);
  const Thyra::VectorBase<Scalar>& b = *(coeff.b[j]);
  const Thyra::VectorBase<Scalar>& c = *(coeff.c[j]);
  const Thyra::VectorBase<Scalar>& d = *(coeff.d[j]);
  
  if (!Teuchos::is_null(S)) {
    // Evaluate S:
    //*S = (a) + (b)*dt + (c)*dt*dt + (d)*dt*dt*dt;
    V_StVpStV(outArg(*S),dt*dt*dt,d,dt*dt,c);
    Vp_StV(outArg(*S),dt,b);
    Vp_StV(outArg(*S),ST::one(),a);
  } 
  if (!Teuchos::is_null(Sp)) {
    // Evaluate S':
    //*Sp = (b) + (c)*2*dt + (d)*3*dt*dt;
    V_StVpStV(outArg(*Sp),Scalar(3*ST::one())*dt*dt,d,Scalar(2*ST::one())*dt,c);
    Vp_StV(outArg(*Sp),ST::one(),b);
  }
  if (!Teuchos::is_null(Spp)) {
    // Evaluate S'':
    //*Spp = (c)*2 + (d)*6*dt;
    V_StVpStV(outArg(*Spp),Scalar(6*ST::one())*dt,d,Scalar(2*ST::one()),c);
  }
}




template<class Scalar>
CubicSplineInterpolator<Scalar>::CubicSplineInterpolator()
{
  splineCoeffComputed_ = false;
  nodesSet_ = false;
}


template<class Scalar>
bool CubicSplineInterpolator<Scalar>::supportsCloning() const
{
  return true;
}


template<class Scalar>
RCP<InterpolatorBase<Scalar> >
CubicSplineInterpolator<Scalar>::cloneInterpolator() const
{
  RCP<CubicSplineInterpolator<Scalar> >
    interpolator = Teuchos::rcp(new CubicSplineInterpolator<Scalar>);
  if (!is_null(parameterList_))
    interpolator->parameterList_ = parameterList(*parameterList_);
  return interpolator;
}

template<class Scalar>
void CubicSplineInterpolator<Scalar>::setNodes(
    const RCP<const typename DataStore<Scalar>::DataStoreVector_t> & nodesPtr
    )
{
  nodes_ = nodesPtr;
  nodesSet_ = true;
  splineCoeffComputed_ = false;
#ifdef RYTHMOS_DEBUG
  const typename DataStore<Scalar>::DataStoreVector_t & nodes = *nodesPtr;
  // Copy nodes to internal data structure for verification upon calls to interpolate
  nodes_copy_ = Teuchos::rcp(new typename DataStore<Scalar>::DataStoreVector_t);
  nodes_copy_->reserve(nodes.size());
  for (int i=0 ; i<Teuchos::as<int>(nodes.size()) ; ++i) {
    nodes_copy_->push_back(*nodes[i].clone());
  }
#endif // RYTHMOS_DEBUG
}

template<class Scalar>
void CubicSplineInterpolator<Scalar>::interpolate(
  const Array<Scalar> &t_values,
  typename DataStore<Scalar>::DataStoreVector_t *data_out
  ) const
{
  using Teuchos::as;
  using Teuchos::outArg;
  typedef Teuchos::ScalarTraits<Scalar> ST;

  TEST_FOR_EXCEPTION( nodesSet_ == false, std::logic_error,
      "Error!, setNodes must be called before interpolate"
      );
#ifdef RYTHMOS_DEBUG
  // Check that our nodes_ have not changed between the call to setNodes and interpolate
  assertNodesUnChanged<Scalar>(*nodes_,*nodes_copy_);
  // Assert that the base interpolator preconditions are satisfied
  assertBaseInterpolatePreconditions(*nodes_,t_values,data_out);
#endif // RYTHMOS_DEBUG
  
  // Output info
  const RCP<FancyOStream> out = this->getOStream();
  const Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel(); 
  Teuchos::OSTab ostab(out,1,"CSI::interpolator");
  if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
    *out << "nodes_:" << std::endl;
    for (Teuchos::Ordinal i=0 ; i<(*nodes_).size() ; ++i) {
      *out << "nodes_[" << i << "] = " << std::endl;
      (*nodes_)[i].describe(*out,Teuchos::VERB_EXTREME);
    }
    *out << "t_values = " << std::endl;
    for (Teuchos::Ordinal i=0 ; i<t_values.size() ; ++i) {
      *out << "t_values[" << i << "] = " << t_values[i] << std::endl;
    }
  }

  data_out->clear();

  // Return immediately if no time points are requested ...
  if (t_values.size() == 0) {
    return;
  }

  if ((*nodes_).size() == 1) {
    // trivial case of one node.  Preconditions assert that t_values[0] ==
    // (*nodes_)[0].time so we can just pass it out
    DataStore<Scalar> DS((*nodes_)[0]);
    data_out->push_back(DS);
  }
  else { // (*nodes_).size() >= 2
    int n = 0; // index into t_values
    // Loop through all of the time interpolation points in the buffer and
    // satisfiy all of the requested time points that you find.  NOTE: The
    // loop will be existed once all of the time points are satisified (see
    // return below).
    for (Teuchos::Ordinal i=0 ; i < (*nodes_).size()-1; ++i) {
      const Scalar& ti = (*nodes_)[i].time;
      const Scalar& tip1 = (*nodes_)[i+1].time;
      const TimeRange<Scalar> range_i(ti,tip1);
      // For the interpolation range of [ti,tip1], satisify all of the
      // requested points in this range.
      while ( range_i.isInRange(t_values[n]) ) {
        // First we check for exact node matches:
        if (compareTimeValues(t_values[n],ti)==0) {
          DataStore<Scalar> DS((*nodes_)[i]);
          data_out->push_back(DS);
        }
        else if (compareTimeValues(t_values[n],tip1)==0) {
          DataStore<Scalar> DS((*nodes_)[i+1]);
          data_out->push_back(DS);
        } else {
          if (!splineCoeffComputed_) {
            computeCubicSplineCoeff<Scalar>(*nodes_,outArg(splineCoeff_));
            splineCoeffComputed_ = true;
          }
          DataStore<Scalar> DS;
          RCP<Thyra::VectorBase<Scalar> > x = createMember((*nodes_)[i].x->space());
          evaluateCubicSpline<Scalar>( splineCoeff_, i, t_values[n], outArg(*x) );
          DS.time = t_values[n];
          DS.x = x;
          DS.accuracy = ST::zero();
          data_out->push_back(DS);
        }
        // Move to the next user time point to consider!
        n++;
        if (n == as<int>(t_values.size())) {
          // WE ARE ALL DONE!  MOVE OUT!
          return;
        }
      }
      // Move on the the next interpolation time range
    }
  } // (*nodes_).size() == 1
}


template<class Scalar>
int CubicSplineInterpolator<Scalar>::order() const
{
  return(1);
}


template<class Scalar>
std::string CubicSplineInterpolator<Scalar>::description() const
{
  std::string name = "Rythmos::CubicSplineInterpolator";
  return(name);
}


template<class Scalar>
void CubicSplineInterpolator<Scalar>::describe(
  FancyOStream &out,
  const Teuchos::EVerbosityLevel verbLevel
  ) const
{
  using Teuchos::as;
  if ( (as<int>(verbLevel) == as<int>(Teuchos::VERB_DEFAULT) ) ||
       (as<int>(verbLevel) >= as<int>(Teuchos::VERB_LOW)     )
     )
  {
    out << description() << "::describe" << std::endl;
  }
  else if (as<int>(verbLevel) >= as<int>(Teuchos::VERB_LOW))
  {}
  else if (as<int>(verbLevel) >= as<int>(Teuchos::VERB_MEDIUM))
  {}
  else if (as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH))
  {}
}


template <class Scalar>
void CubicSplineInterpolator<Scalar>::setParameterList(
  RCP<ParameterList> const& paramList
  )
{
  TEST_FOR_EXCEPT(is_null(paramList));
  paramList->validateParametersAndSetDefaults(*this->getValidParameters());
  parameterList_ = paramList;
  Teuchos::readVerboseObjectSublist(&*parameterList_,this);
}


template <class Scalar>
RCP<ParameterList>
CubicSplineInterpolator<Scalar>::getNonconstParameterList()
{
  return(parameterList_);
}


template <class Scalar>
RCP<ParameterList>
CubicSplineInterpolator<Scalar>::unsetParameterList()
{
  RCP<ParameterList> temp_param_list;
  std::swap( temp_param_list, parameterList_ );
  return(temp_param_list);
}

template<class Scalar>
RCP<const Teuchos::ParameterList>
CubicSplineInterpolator<Scalar>::getValidParameters() const
{
  static RCP<Teuchos::ParameterList> validPL;
  if (is_null(validPL)) {
    RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    Teuchos::setupVerboseObjectSublist(&*pl);
    validPL = pl;
  }
  return (validPL);
}


// 
// Explicit Instantiation macro
//
// Must be expanded from within the Rythmos namespace!
//


#define RYTHMOS_CUBIC_SPLINE_INTERPOLATOR_INSTANT(SCALAR) \
  \
  template class CubicSplineInterpolator< SCALAR >; \
  \
  template class CubicSplineCoeff< SCALAR >; \
  template RCP<CubicSplineInterpolator< SCALAR > > cubicSplineInterpolator(); \
  template void computeCubicSplineCoeff( \
      const DataStore< SCALAR >::DataStoreVector_t & data, \
      const Ptr<CubicSplineCoeff< SCALAR > > & coeffPtr \
      ); \
  template void validateCubicSplineCoeff(const CubicSplineCoeff< SCALAR >& coeff); \
  template void evaluateCubicSpline( \
      const CubicSplineCoeff< SCALAR >& coeff, \
      Teuchos::Ordinal j,  \
      const  SCALAR & t, \
      const Ptr<Thyra::VectorBase< SCALAR > >& S, \
      const Ptr<Thyra::VectorBase< SCALAR > >& Sp,  \
      const Ptr<Thyra::VectorBase< SCALAR > >& Spp \
      ); 



} // namespace Rythmos


#endif // Rythmos_CUBIC_SPLINE_INTERPOLATOR_DEF_H
