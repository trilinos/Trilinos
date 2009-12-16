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

#ifndef Rythmos_INTERPOLATION_BUFFER_DEF_H
#define Rythmos_INTERPOLATION_BUFFER_DEF_H

#include "Rythmos_InterpolationBuffer_decl.hpp"
#include "Rythmos_InterpolatorBaseHelpers.hpp"
#include "Rythmos_LinearInterpolator.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"

namespace {

  static std::string IBPolicyTypeInvalid_name = "Invalid Policy";
  static std::string IBPolicyTypeStatic_name = "Static Policy";
  static std::string IBPolicyTypeKeepNewest_name = "Keep Newest Policy";
  static std::string interpolationBufferPolicySelection_name = "InterpolationBufferPolicy";
  static std::string interpolationBufferPolicySelection_default = IBPolicyTypeKeepNewest_name;

  static std::string interpolationBufferStorageLimit_name = "StorageLimit";
  static int interpolationBufferStorageLimit_default = 0;

  Teuchos::Array<std::string>
    S_InterpolationBufferPolicyTypes = Teuchos::tuple<std::string>(
        IBPolicyTypeInvalid_name,
        IBPolicyTypeStatic_name,
        IBPolicyTypeKeepNewest_name
        );

  const Teuchos::RCP<Teuchos::StringToIntegralParameterEntryValidator<Rythmos::IBPolicy> >
    interpolationBufferPolicyValidator = Teuchos::rcp(
        new Teuchos::StringToIntegralParameterEntryValidator<Rythmos::IBPolicy>(
          S_InterpolationBufferPolicyTypes,
          Teuchos::tuple<Rythmos::IBPolicy>(
            Rythmos::BUFFER_POLICY_INVALID,
            Rythmos::BUFFER_POLICY_STATIC,
            Rythmos::BUFFER_POLICY_KEEP_NEWEST
            ),
          interpolationBufferPolicySelection_name
          )
        );

} // namespace


namespace Rythmos {

// ////////////////////////////
// Defintions


template<class Scalar>
InterpolationBuffer<Scalar>::InterpolationBuffer()
{
  this->defaultInitializeAll_();
  initialize(Teuchos::null,0);
}

template<class Scalar>
void InterpolationBuffer<Scalar>::defaultInitializeAll_()
{
  interpolator_ = Teuchos::null;
  storage_limit_ = -1;
  data_vec_ = Teuchos::null;
  paramList_ = Teuchos::null;
  policy_ = BUFFER_POLICY_INVALID;
}

template<class Scalar>
RCP<const Thyra::VectorSpaceBase<Scalar> >
InterpolationBuffer<Scalar>::get_x_space() const
{
  if (data_vec_->size() == 0) {
    RCP<const Thyra::VectorSpaceBase<Scalar> > space;
    return(space);
  } else {
    return((*data_vec_)[0].x->space());
  }
}


template<class Scalar>
void InterpolationBuffer<Scalar>::initialize( 
  const RCP<InterpolatorBase<Scalar> >& interpolator
  ,int storage
  )
{
  RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"IB::initialize");
  if ( Teuchos::as<int>(this->getVerbLevel()) >= Teuchos::as<int>(Teuchos::VERB_HIGH) ) {
    *out << "Initializing InterpolationBuffer" << std::endl;
    *out << "Calling setInterpolator..." << std::endl;
  }
  data_vec_ = rcp(new typename DataStore<Scalar>::DataStoreVector_t);
  setInterpolator(interpolator);
  if ( Teuchos::as<int>(this->getVerbLevel()) >= Teuchos::as<int>(Teuchos::VERB_HIGH) ) {
    *out << "Calling setStorage..." << std::endl;
  }
  setStorage(storage);
  policy_ = BUFFER_POLICY_KEEP_NEWEST;
}

template<class Scalar>
void InterpolationBuffer<Scalar>::setStorage( int storage )
{
  int storage_limit = std::max(2,storage); // Minimum of two points so interpolation is possible
  TEST_FOR_EXCEPTION(
    Teuchos::as<int>(data_vec_->size()) > storage_limit,
    std::logic_error,
    "Error, specified storage = " << storage_limit
    << " is below current number of vectors stored = " << data_vec_->size() << "!\n"
    );
  storage_limit_ = storage_limit;
  RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"IB::setStorage");
  if ( Teuchos::as<int>(this->getVerbLevel()) >= Teuchos::as<int>(Teuchos::VERB_HIGH) ) {
    *out << "storage_limit = " << storage_limit_ << std::endl;
  }
}


template<class Scalar>
int InterpolationBuffer<Scalar>::getStorage() const
{
  return(storage_limit_);
}


template<class Scalar>
void InterpolationBuffer<Scalar>::setInterpolator(
  const RCP<InterpolatorBase<Scalar> >& interpolator
  )
{
  if (interpolator == Teuchos::null) {
    interpolator_ = linearInterpolator<Scalar>();
  } else {
    interpolator_ = interpolator;
  }
  RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"IB::setInterpolator");
  if ( Teuchos::as<int>(this->getVerbLevel()) >= Teuchos::as<int>(Teuchos::VERB_HIGH) ) {
    *out << "interpolator = " << interpolator_->description() << std::endl;
  }
}

template<class Scalar>
RCP<InterpolatorBase<Scalar> >
  InterpolationBuffer<Scalar>::getNonconstInterpolator()
{
  return interpolator_;
}

template<class Scalar>
RCP<const InterpolatorBase<Scalar> >
  InterpolationBuffer<Scalar>::getInterpolator() const
{
  return interpolator_;
}

template<class Scalar>
RCP<InterpolatorBase<Scalar> > InterpolationBuffer<Scalar>::unSetInterpolator()
{
  RCP<InterpolatorBase<Scalar> > old_interpolator = interpolator_;
  interpolator_ = linearInterpolator<Scalar>();
  return old_interpolator;
}


template<class Scalar>
void InterpolationBuffer<Scalar>::addPoints( 
  const Array<Scalar>& time_vec
  ,const Array<RCP<const Thyra::VectorBase<Scalar> > >& x_vec
  ,const Array<RCP<const Thyra::VectorBase<Scalar> > >& xdot_vec 
  )
{
#ifdef RYTHMOS_DEBUG
  // Check preconditions
  assertTimePointsAreSorted(time_vec);
  int tsize = Teuchos::as<int>(time_vec.size());
  TEST_FOR_EXCEPTION(
    tsize == 0, std::logic_error,
    "Error, time_vec is empty!"
    );
  TEST_FOR_EXCEPTION(
    Teuchos::as<int>(x_vec.size()) != tsize, std::logic_error,
    "Error, size of x_vec = " << x_vec.size() << " != " << tsize << " = size of time_vec!\n"
    );
  TEST_FOR_EXCEPTION(
    Teuchos::as<int>(xdot_vec.size()) != tsize, std::logic_error,
    "Error, size of xdot_vec = " << x_vec.size() << " != " << tsize << " = size of time_vec!\n"
    );
  for (int i=0; i<tsize ; ++i) {
    TEST_FOR_EXCEPTION(
      x_vec[i] == Teuchos::null, std::logic_error,
      "Error, x_vec[" << i << "] == null!\n"
      );
//    TEST_FOR_EXCEPTION(
//      xdot_vec[i] == Teuchos::null, std::logic_error,
//      "Error, xdot_vec[" << i << "] == null!\n"
//      );
  }
  assertNoTimePointsInsideCurrentTimeRange(*this,time_vec);
#endif // RYTHMOS_DEBUG
  RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"IB::addPoints");
  if ( Teuchos::as<int>(this->getVerbLevel()) >= Teuchos::as<int>(Teuchos::VERB_HIGH) ) {
    *out << "time_vec = " << std::endl;
    for (Teuchos::Ordinal i=0 ; i<time_vec.size() ; ++i) {
      *out << "time_vec[" << i << "] = " << time_vec[i] << std::endl;
    }
    *out << "x_vec = " << std::endl;
    for (Teuchos::Ordinal i=0 ; i<x_vec.size() ; ++i) {
      *out << "x_vec[" << i << "] = " << std::endl;
      x_vec[i]->describe(*out,Teuchos::VERB_EXTREME);
    }
    *out << "xdot_vec = " << std::endl;
    for (Teuchos::Ordinal i=0 ; i<xdot_vec.size() ; ++i) {
      if (!is_null(xdot_vec[i])) {
        *out << "xdot_vec[" << i << "] = " << std::endl;
        xdot_vec[i]->describe(*out,Teuchos::VERB_EXTREME);
      }
    }
  }
  typename DataStore<Scalar>::DataStoreList_t input_data_list;
  vectorToDataStoreList<Scalar>(time_vec,x_vec,xdot_vec,&input_data_list);
  // Check that we're not going to exceed our storage limit:
  if (Teuchos::as<int>(data_vec_->size()+input_data_list.size()) > storage_limit_) { 
    if (policy_ == BUFFER_POLICY_STATIC) {
      TEST_FOR_EXCEPTION(
        true, std::logic_error,
        "Error, buffer would be over-full and buffer policy is BUFFER_POLICY_STATIC, these points can not be added\n"
        );
    } else if (policy_ == BUFFER_POLICY_KEEP_NEWEST) {
      if (input_data_list.front() > data_vec_->back()) {
        // Case:  all of new points are past end of existing points
        // Remove points from the beginning of data_vec, then add new points
        int num_extra_points = input_data_list.size()-(storage_limit_-data_vec_->size());
#ifdef RYTHMOS_DEBUG
        TEST_FOR_EXCEPTION( num_extra_points <= 0, std::logic_error, 
            "Error!  Buffer policy is keep newest and input data size = " << input_data_list.size() << ", storage limit  = " << storage_limit_ << ", and data_vec size = " << data_vec_->size() << ".  Somehow number of points to delete = " << num_extra_points << " <= 0!"
            );
#endif // RYTHMOS_DEBUG
        typename DataStore<Scalar>::DataStoreVector_t::iterator 
          data_it = data_vec_->begin();
        for (int i=0 ; i < num_extra_points ; ++i) {
          data_it++;
        }
        if ( Teuchos::as<int>(this->getVerbLevel()) >= Teuchos::as<int>(Teuchos::VERB_HIGH) ) {
          *out << "Removing " << num_extra_points 
               << " from beginning of data_vec to make room for new points." << std::endl;
        }
        data_vec_->erase(data_vec_->begin(),data_it);
      } else if (input_data_list.back() < data_vec_->front()) {
        // Case:  all of new points are before beginning of existing points
        // Remove points from end of data_vec, then add new points
        int num_extra_points = input_data_list.size()-(storage_limit_-data_vec_->size());
#ifdef RYTHMOS_DEBUG
        TEST_FOR_EXCEPTION( num_extra_points <= 0, std::logic_error, 
            "Error!  Buffer policy is keep newest and input data size = " << input_data_list.size() << ", storage limit  = " << storage_limit_ << ", and data_vec size = " << data_vec_->size() << ".  Somehow number of points to delete = " << num_extra_points << " <= 0!"
            );
#endif // RYTHMOS_DEBUG
        typename DataStore<Scalar>::DataStoreVector_t::iterator 
          data_it = data_vec_->end();
        for (int i=0 ; i < num_extra_points ; ++i) {
          data_it--;
        }
        if ( Teuchos::as<int>(this->getVerbLevel()) >= Teuchos::as<int>(Teuchos::VERB_HIGH) ) {
          *out << "Removing " << num_extra_points 
               << " from end of data_vec to make room for new points." << std::endl;
        }
        data_vec_->erase(data_it,data_vec_->end());
      } else {
        // Case:  Some points are before beginning of data_vec and some points are after end of data_vec
        TEST_FOR_EXCEPTION(
          true, std::logic_error,
          "Error, incoming points are on both sides of TimeRange, I don't know which points are newest in this case.\n"
          );
      }
    } else {
      // Unknown Buffer policy:
      TEST_FOR_EXCEPTION(
        true, std::logic_error,
        "Error, unknown buffer policy.\n"
        );
    }
  }
  // Clone the vectors in input_data_list
  std::list<DataStore<Scalar> > internal_input_data_list;
  typename DataStore<Scalar>::DataStoreList_t::iterator it_list;
  for (it_list = input_data_list.begin() ; it_list != input_data_list.end() ; it_list++) {
    RCP<DataStore<Scalar> > ds_clone = it_list->clone();
    internal_input_data_list.push_back(*ds_clone);
  }
  // Now add all the remaining points to data_vec
  data_vec_->insert(data_vec_->end(),internal_input_data_list.begin(),internal_input_data_list.end());
  // And sort data_vec:
  std::sort(data_vec_->begin(),data_vec_->end());
  if ( Teuchos::as<int>(this->getVerbLevel()) >= Teuchos::as<int>(Teuchos::VERB_HIGH) ) {
    *out << "data_vec at end of addPoints:" << std::endl;
    for (Teuchos::Ordinal i=0 ; i<data_vec_->size() ; ++i) {
      *out << "data_vec[" << i << "] = " << std::endl;
      (*data_vec_)[i].describe(*out,Teuchos::VERB_EXTREME);
    }
  }
}


template<class Scalar>
void InterpolationBuffer<Scalar>::getPoints(
  const Array<Scalar>& time_vec
  ,Array<RCP<const Thyra::VectorBase<Scalar> > >* x_vec
  ,Array<RCP<const Thyra::VectorBase<Scalar> > >* xdot_vec
  ,Array<ScalarMag>* accuracy_vec
  ) const
{
  RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"IB::getPoints");
  typename DataStore<Scalar>::DataStoreVector_t data_out;
  interpolate<Scalar>(*interpolator_, data_vec_, time_vec, &data_out);
  Array<Scalar> time_out_vec;
  dataStoreVectorToVector<Scalar>(data_out, &time_out_vec, x_vec, xdot_vec, accuracy_vec);
  TEST_FOR_EXCEPTION(
    (time_vec.size() != time_out_vec.size()), std::logic_error,
    "Error, number of time points returned from interpolator = " <<
    time_out_vec.size() << " != " << time_vec.size() << 
    " = number of time points requested\n"
    );
  if ( Teuchos::as<int>(this->getVerbLevel()) >= Teuchos::as<int>(Teuchos::VERB_HIGH) ) {
    *out << "Conversion of DataStoreVector to Vector successful" << std::endl;
  }
}


template<class Scalar>
TimeRange<Scalar> InterpolationBuffer<Scalar>::getTimeRange() const
{
  TimeRange<Scalar> timerange;
  if (data_vec_->size() > 0) {
    timerange = TimeRange<Scalar>(data_vec_->front().time,data_vec_->back().time);
  }
  return(timerange);
}


template<class Scalar>
void InterpolationBuffer<Scalar>::getNodes( Array<Scalar>* time_vec ) const
{
  int N = data_vec_->size();
  time_vec->clear();
  time_vec->reserve(N);
  for (int i=0 ; i<N ; ++i) {
    time_vec->push_back((*data_vec_)[i].time);
  }
  RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"IB::getNodes");
  if ( Teuchos::as<int>(this->getVerbLevel()) >= Teuchos::as<int>(Teuchos::VERB_HIGH) ) {
    *out << this->description() << std::endl;
    for (Teuchos::Ordinal i=0 ; i<time_vec->size() ; ++i) {
      *out << "time_vec[" << i << "] = " << (*time_vec)[i] << std::endl;
    }
  }
}


template<class Scalar>
void InterpolationBuffer<Scalar>::removeNodes( Array<Scalar>& time_vec ) 
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  int N = time_vec.size();
#ifdef RYTHMOS_DEBUG
  // Check preconditions:
  TimeRange<Scalar> range = this->getTimeRange();
  for (int i=0; i<N ; ++i) {
    TEST_FOR_EXCEPTION(
      ~(range.lower() <= time_vec[i]) && (time_vec[i] <= range.upper()),
      std::logic_error,
      "Error, time_vec[" << i << "] = " << time_vec[i] << 
      "is not in range of this interpolation buffer = [" << 
      range.lower() << "," << range.upper() << "]!\n"
      );
  }
#endif // RYTHMOS_DEBUG
  RCP<Thyra::VectorBase<Scalar> > vec_temp;
  ScalarMag z = ST::zero();
  for (int i=0; i<N ; ++i) {
    DataStore<Scalar> ds_temp(time_vec[i],vec_temp,vec_temp,z);
    typename DataStore<Scalar>::DataStoreVector_t::iterator 
      data_it = std::find(data_vec_->begin(),data_vec_->end(),ds_temp);
    TEST_FOR_EXCEPTION(
      data_it == data_vec_->end(), std::logic_error,
      "Error, time_vec[" << i << "] = " << time_vec[i] << "is not a node in the interpolation buffer!\n"
      );
    data_vec_->erase(data_it);
  }
}


template<class Scalar>
int InterpolationBuffer<Scalar>::getOrder() const
{
  return(interpolator_->order());
}


template<class Scalar>
std::string InterpolationBuffer<Scalar>::description() const
{
  std::string name = "Rythmos::InterpolationBuffer";
  return(name);
}


template<class Scalar>
void InterpolationBuffer<Scalar>::describe(
  Teuchos::FancyOStream                &out
  ,const Teuchos::EVerbosityLevel      verbLevel
  ) const
{
  if ( (Teuchos::as<int>(verbLevel) == Teuchos::as<int>(Teuchos::VERB_DEFAULT) ) ||
    (Teuchos::as<int>(verbLevel) >= Teuchos::as<int>(Teuchos::VERB_LOW)     )
    ) {
    out << description() << "::describe" << std::endl;
    out << "interpolator = " << interpolator_->description() << std::endl;
    out << "storage_limit = " << storage_limit_ << std::endl;
  } else if (Teuchos::as<int>(verbLevel) >= Teuchos::as<int>(Teuchos::VERB_LOW)) {
  } else if (Teuchos::as<int>(verbLevel) >= Teuchos::as<int>(Teuchos::VERB_MEDIUM)) {
  } else if (Teuchos::as<int>(verbLevel) >= Teuchos::as<int>(Teuchos::VERB_HIGH)) {
    out << "data_vec = " << std::endl;
    for (Teuchos::Ordinal i=0; i<data_vec_->size() ; ++i) {
      out << "data_vec[" << i << "] = " << std::endl;
      (*data_vec_)[i].describe(out,this->getVerbLevel());
    }
  }
}


template <class Scalar>
void InterpolationBuffer<Scalar>::setParameterList(RCP<Teuchos::ParameterList> const& paramList)
{
  TEST_FOR_EXCEPT( is_null(paramList) );
  paramList->validateParameters(*this->getValidParameters());
  paramList_ = paramList;

  Teuchos::readVerboseObjectSublist(&*paramList_,this);

  //int outputLevel = paramList_->get( "outputLevel", int(-1) );
  //outputLevel = std::min(std::max(outputLevel,-1),4);
  //this->setVerbLevel(static_cast<Teuchos::EVerbosityLevel>(outputLevel));
  IBPolicy policyLevel = interpolationBufferPolicyValidator->getIntegralValue(
      *paramList_, interpolationBufferPolicySelection_name, interpolationBufferPolicySelection_default
      );
  if (policyLevel != BUFFER_POLICY_INVALID) {
    policy_ = policyLevel;
  }
  int storage_limit = paramList_->get( interpolationBufferStorageLimit_name, interpolationBufferStorageLimit_default);
  setStorage(storage_limit);
}

template<class Scalar>
RCP<const Teuchos::ParameterList> InterpolationBuffer<Scalar>::getValidParameters() const
{
  static RCP<Teuchos::ParameterList> validPL;

  if (is_null(validPL)) {

    RCP<Teuchos::ParameterList>
      pl = Teuchos::parameterList();

    Teuchos::setupVerboseObjectSublist(&*pl);

    pl->set( 
        interpolationBufferPolicySelection_name, 
        interpolationBufferPolicySelection_default,
        "Interpolation Buffer Policy for when the maximum storage size is exceeded.  Static will throw an exception when the storage limit is exceeded.  Keep Newest will over-write the oldest data in the buffer when the storage limit is exceeded.",
        interpolationBufferPolicyValidator
        );

    pl->set(
        interpolationBufferStorageLimit_name,
        interpolationBufferStorageLimit_default,
        "Storage limit for the interpolation buffer."
        );

    validPL = pl;

  }
  return validPL;
}


template <class Scalar>
RCP<Teuchos::ParameterList>
InterpolationBuffer<Scalar>::getNonconstParameterList()
{
  return(paramList_);
}


template <class Scalar>
RCP<Teuchos::ParameterList> InterpolationBuffer<Scalar>::unsetParameterList()
{
  RCP<Teuchos::ParameterList> temp_param_list = paramList_;
  paramList_ = Teuchos::null;
  return(temp_param_list);
}

template <class Scalar>
IBPolicy InterpolationBuffer<Scalar>::getIBPolicy()
{
  return policy_;
}

// 
// Explicit Instantiation macro
//
// Must be expanded from within the Rythmos namespace!
//

#define RYTHMOS_INTERPOLATION_BUFFER_INSTANT(SCALAR) \
  \
  template class InterpolationBuffer< SCALAR >; \
  \
  template RCP<InterpolationBuffer< SCALAR > > interpolationBuffer(  \
    const RCP<InterpolatorBase< SCALAR > >& interpolator, \
    int storage  \
    ); 

} // namespace Rythmos


#endif // Rythmos_INTERPOLATION_BUFFER_DEF_H
