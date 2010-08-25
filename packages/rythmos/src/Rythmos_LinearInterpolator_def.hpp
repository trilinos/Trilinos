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

#ifndef Rythmos_LINEAR_INTERPOLATOR_DEF_H
#define Rythmos_LINEAR_INTERPOLATOR_DEF_H

#include "Rythmos_LinearInterpolator_decl.hpp"
#include "Rythmos_InterpolatorBaseHelpers.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_VectorSpaceBase.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"


namespace Rythmos {


// non-member constructor
template<class Scalar>
RCP<LinearInterpolator<Scalar> > linearInterpolator()
{
  RCP<LinearInterpolator<Scalar> > li = rcp(new LinearInterpolator<Scalar>() );
  return li;
}

template<class Scalar>
LinearInterpolator<Scalar>::LinearInterpolator()
{
  nodes_ = Teuchos::null;
  parameterList_ = Teuchos::null;
}


template<class Scalar>
bool LinearInterpolator<Scalar>::supportsCloning() const
{
  return true;
}


template<class Scalar>
RCP<InterpolatorBase<Scalar> >
LinearInterpolator<Scalar>::cloneInterpolator() const
{
  RCP<LinearInterpolator<Scalar> >
    interpolator = Teuchos::rcp(new LinearInterpolator<Scalar>);
  if (!is_null(parameterList_))
    interpolator->parameterList_ = parameterList(*parameterList_);
  return interpolator;
}

template<class Scalar>
void LinearInterpolator<Scalar>::setNodes(
    const RCP<const typename DataStore<Scalar>::DataStoreVector_t> & nodes
    )
{
  nodes_ = nodes;
}


template<class Scalar>
void LinearInterpolator<Scalar>::interpolate(
  const Array<Scalar> &t_values,
  typename DataStore<Scalar>::DataStoreVector_t *data_out
  ) const
{

  using Teuchos::as;
  typedef Teuchos::ScalarTraits<Scalar> ST;

#ifdef RYTHMOS_DEBUG
  assertBaseInterpolatePreconditions(*nodes_,t_values,data_out);
#endif // RYTHMOS_DEBUG
  
  // Output info
  const RCP<FancyOStream> out = this->getOStream();
  const Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel(); 
  Teuchos::OSTab ostab(out,1,"LI::interpolator");
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

  // Return immediatly if not time points are requested ...
  if (t_values.size() == 0) {
    if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
      *out << "Returning because no time points were requested" << std::endl;
    }
    return;
  }

  if ((*nodes_).size() == 1) {
    // trivial case of one node.  Preconditions assert that t_values[0] ==
    // (*nodes_)[0].time so we can just pass it out
    DataStore<Scalar> DS((*nodes_)[0]);
    data_out->push_back(DS);
    if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
      *out << "Only a single node is in the buffer, so preconditions assert that this must be the point requested" << std::endl;
    }
  }
  else { // (*nodes_).size() >= 2
    if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
      *out << "More than two nodes, looping through the intervals looking for points to interpolate" << std::endl;
    }
    int n = 0; // index into t_values
    // Loop through all of the time interpolation points in the buffer and
    // satisfiy all of the requested time points that you find.  NOTE: The
    // loop will be existed once all of the time points are satisified (see
    // return below).
    for (int i=0 ; i < as<int>((*nodes_).size())-1; ++i) {
      if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
        *out << "Looking for interval containing: t_values["<<n<<"] = " << t_values[n] << std::endl;
      }
      const Scalar& ti = (*nodes_)[i].time;
      const Scalar& tip1 = (*nodes_)[i+1].time;
      const Scalar  h = tip1-ti;
      const TimeRange<Scalar> range_i(ti,tip1);
      // For the interpolation range of [ti,tip1], satisify all of the
      // requested points in this range.
      while ( range_i.isInRange(t_values[n]) ) {
        // First we check for exact node matches:
        if (compareTimeValues(t_values[n],ti)==0) {
          DataStore<Scalar> DS((*nodes_)[i]);
          data_out->push_back(DS);
          if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
            *out << "Found an exact node match (on left), shallow copying." << std::endl;
            *out << "Found t_values["<<n<<"] = " << t_values[n] << " on boundary of interval ["<<ti<<","<<tip1<<"]" << std::endl;
          }
        }
        else if (compareTimeValues(t_values[n],tip1)==0) {
          DataStore<Scalar> DS((*nodes_)[i+1]);
          data_out->push_back(DS);
          if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
            *out << "Found an exact node match (on right), shallow copying." << std::endl;
            *out << "Found t_values["<<n<<"] = " << t_values[n] << " on boundary of interval ["<<ti<<","<<tip1<<"]" << std::endl;
          }
        }
        else {
          if ( as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH) ) {
            *out << "Interpolating a point (creating a new vector)..." << std::endl;
            *out << "Found t_values["<<n<<"] = " << t_values[n] << " in interior of interval ["<<ti<<","<<tip1<<"]" << std::endl;
          }
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
          const Scalar& t = t_values[n];
          DS.time = t;
          // Get the time and interpolation node points
          RCP<const Thyra::VectorBase<Scalar> > xi = (*nodes_)[i].x;
          RCP<const Thyra::VectorBase<Scalar> > xip1 = (*nodes_)[i+1].x;
          RCP<const Thyra::VectorBase<Scalar> > xdoti = (*nodes_)[i].xdot;
          RCP<const Thyra::VectorBase<Scalar> > xdotip1 = (*nodes_)[i+1].xdot;
          // Get constants used in interplation
          const Scalar dt = t-ti;
          const Scalar dt_over_h = dt / h;
          const Scalar one_minus_dt_over_h = ST::one() - dt_over_h;
          // x = dt/h * xip1 + (1-dt/h) * xi
          RCP<Thyra::VectorBase<Scalar> > x;
          if (!is_null(xi) && !is_null(xip1)) {
            x = createMember(xi->space());
            Thyra::V_StVpStV(x.ptr(),dt_over_h,*xip1,one_minus_dt_over_h,*xi);
          }
          DS.x = x;
          // x = dt/h * xdotip1 + (1-dt/h) * xdoti
          RCP<Thyra::VectorBase<Scalar> > xdot;
          if (!is_null(xdoti) && !is_null(xdotip1)) {
            xdot = createMember(xdoti->space());
            Thyra::V_StVpStV(xdot.ptr(),dt_over_h,*xdotip1,one_minus_dt_over_h,*xdoti);
          }
          DS.xdot = xdot;
          // Estimate our accuracy ???
          DS.accuracy = h;
          // 2007/12/06: rabartl: Above, should the be a relative value of
          // some type.  What does this really mean?
          // Store this interplation
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
void LinearInterpolator<Scalar>::setParameterList(
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
LinearInterpolator<Scalar>::getNonconstParameterList()
{
  return(parameterList_);
}


template <class Scalar>
RCP<ParameterList>
LinearInterpolator<Scalar>::unsetParameterList()
{
  RCP<ParameterList> temp_param_list;
  std::swap( temp_param_list, parameterList_ );
  return(temp_param_list);
}

template<class Scalar>
RCP<const Teuchos::ParameterList> LinearInterpolator<Scalar>::getValidParameters() const
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

#define RYTHMOS_LINEAR_INTERPOLATOR_INSTANT(SCALAR) \
  \
  template class LinearInterpolator< SCALAR >; \
  \
  template RCP<LinearInterpolator< SCALAR > > linearInterpolator(); 

} // namespace Rythmos


#endif // Rythmos_LINEAR_INTERPOLATOR_DEF_H
