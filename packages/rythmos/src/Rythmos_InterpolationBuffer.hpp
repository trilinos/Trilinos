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

#ifndef Rythmos_INTERPOLATION_BUFFER_H
#define Rythmos_INTERPOLATION_BUFFER_H

#include "Rythmos_InterpolationBufferBase.hpp"
#include "Rythmos_InterpolatorBase.hpp"
#include "Rythmos_InterpolatorBaseHelpers.hpp"
#include "Rythmos_LinearInterpolator.hpp"
#include "Rythmos_DataStore.hpp"

#include "Thyra_VectorBase.hpp"
#include "Thyra_VectorStdOps.hpp"

namespace Rythmos {


/** \brief concrete class for interpolation buffer functionality. */
template<class Scalar> 
class InterpolationBuffer : virtual public InterpolationBufferBase<Scalar>
{
public:

  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;

  /// Redefined from Rythmos::InterpolationBufferBase
  /** \brief. */
  RCP<const Thyra::VectorSpaceBase<Scalar> > get_x_space() const;
    
  /** \brief. */
  InterpolationBuffer();
  InterpolationBuffer( const RCP<InterpolatorBase<Scalar> >& interpolator_, int storage_ );

  /// Initialize the buffer:
  void initialize( const RCP<InterpolatorBase<Scalar> >& interpolator_, int storage_ );

  /// Set the interpolator for this buffer
  void setInterpolator(const RCP<InterpolatorBase<Scalar> >& interpolator_);
    
  /// Unset the interpolator for this buffer
  RCP<InterpolatorBase<Scalar> >& unSetInterpolator();

  /// Set the maximum storage of this buffer
  void setStorage( int storage );
    
  /// Get the maximum storage of this buffer
  int getStorage() const;
        
  /// Destructor
  ~InterpolationBuffer() {};

  /// Add point to buffer
  void addPoints(
    const Array<Scalar>& time_vec
    ,const Array<RCP<const Thyra::VectorBase<Scalar> > >& x_vec
    ,const Array<RCP<const Thyra::VectorBase<Scalar> > >& xdot_vec);

  /// Get value from buffer
  void getPoints(
    const Array<Scalar>& time_vec
    ,Array<RCP<const Thyra::VectorBase<Scalar> > >* x_vec
    ,Array<RCP<const Thyra::VectorBase<Scalar> > >* xdot_vec
    ,Array<ScalarMag>* accuracy_vec
    ) const;

  /** \brief . */
  TimeRange<Scalar> getTimeRange() const;

  /// Get interpolation nodes
  void getNodes(Array<Scalar>* time_vec) const;

  /// Get order of interpolation
  int getOrder() const;
    
  /// Remove interpolation nodes
  void removeNodes(Array<Scalar>& time_vec);

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
  RCP<Teuchos::ParameterList> getNonconstParameterList();

  /** \brief . */
  RCP<Teuchos::ParameterList> unsetParameterList();

  enum IBPolicy {
    BUFFER_STATIC = 0
    ,BUFFER_KEEP_NEWEST = 1
  };
    
private:

  RCP<InterpolatorBase<Scalar> > interpolator_;
  int storage_limit_;
  RCP<typename DataStore<Scalar>::DataStoreVector_t> data_vec_;

  RCP<Teuchos::ParameterList> parameterList_;

  IBPolicy policy_;

};

// Nonmember constructor
template<class Scalar>
RCP<InterpolationBuffer<Scalar> > interpolationBuffer() 
{
  RCP<InterpolationBuffer<Scalar> > ib = rcp(new InterpolationBuffer<Scalar>());
  return ib;
}

// ////////////////////////////
// Defintions


template<class Scalar>
InterpolationBuffer<Scalar>::InterpolationBuffer()
{
  initialize(Teuchos::null,0);
}


template<class Scalar>
InterpolationBuffer<Scalar>::InterpolationBuffer( 
  const RCP<InterpolatorBase<Scalar> >& interpolator_
  ,int storage_ 
  )
{
  initialize(interpolator_,storage_);
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
  const RCP<InterpolatorBase<Scalar> >& interpolator_
  ,int storage_
  )
{
  RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"IB::initialize");
  if ( Teuchos::as<int>(this->getVerbLevel()) >= Teuchos::as<int>(Teuchos::VERB_HIGH) ) {
    *out << "Initializing InterpolationBuffer" << std::endl;
    *out << "Calling setInterpolator..." << std::endl;
  }
  data_vec_ = rcp(new typename DataStore<Scalar>::DataStoreVector_t);
  setInterpolator(interpolator_);
  if ( Teuchos::as<int>(this->getVerbLevel()) >= Teuchos::as<int>(Teuchos::VERB_HIGH) ) {
    *out << "Calling setStorage..." << std::endl;
  }
  setStorage(storage_);
  policy_ = BUFFER_KEEP_NEWEST;
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
  if (interpolator_ == Teuchos::null) {
    interpolator_ = Teuchos::rcp(new LinearInterpolator<Scalar>);
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
void InterpolationBuffer<Scalar>::addPoints( 
  const Array<Scalar>& time_vec
  ,const Array<RCP<const Thyra::VectorBase<Scalar> > >& x_vec
  ,const Array<RCP<const Thyra::VectorBase<Scalar> > >& xdot_vec 
  )
{
#ifdef TEUCHOS_DEBUG
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
    TEST_FOR_EXCEPTION(
      xdot_vec[i] == Teuchos::null, std::logic_error,
      "Error, xdot_vec[" << i << "] == null!\n"
      );
  }
  assertNoTimePointsInsideCurrentTimeRange(*this,time_vec);
#endif // TEUCHOS_DEBUG
  RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"IB::addPoints");
  if ( Teuchos::as<int>(this->getVerbLevel()) >= Teuchos::as<int>(Teuchos::VERB_HIGH) ) {
    *out << "time_vec = " << std::endl;
    for (unsigned int i=0 ; i<time_vec.size() ; ++i) {
      *out << "time_vec[" << i << "] = " << time_vec[i] << std::endl;
    }
    *out << "x_vec = " << std::endl;
    for (unsigned int i=0 ; i<x_vec.size() ; ++i) {
      *out << "x_vec[" << i << "] = " << std::endl;
      x_vec[i]->describe(*out,Teuchos::VERB_EXTREME);
    }
    *out << "xdot_vec = " << std::endl;
    for (unsigned int i=0 ; i<xdot_vec.size() ; ++i) {
      *out << "xdot_vec[" << i << "] = " << std::endl;
      xdot_vec[i]->describe(*out,Teuchos::VERB_EXTREME);
    }
  }
  typename DataStore<Scalar>::DataStoreList_t input_data_list;
  vectorToDataStoreList<Scalar>(time_vec,x_vec,xdot_vec,&input_data_list);
  // Check that we're not going to exceed our storage limit:
  if ((data_vec_->size()+input_data_list.size()) > Teuchos::as<unsigned int>(storage_limit_)) { 
    if (policy_ == BUFFER_STATIC) {
      TEST_FOR_EXCEPTION(
        true, std::logic_error,
        "Error, buffer is full and buffer policy is BUFFER_STATIC, no points can be added\n"
        );
    } else if (policy_ == BUFFER_KEEP_NEWEST) {
      if (input_data_list.front() > data_vec_->back()) {
        // Case:  all of new points are past end of existing points
        // Remove points from the beginning of data_vec, then add new points
        int num_extra_points = input_data_list.size();
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
        int num_extra_points = input_data_list.size();
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
          "Error, incoming points are on both sides of TimeRange, this feature not implemented yet.\n"
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
    for (unsigned int i=0 ; i<data_vec_->size() ; ++i) {
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
    for (unsigned int i=0 ; i<time_vec->size() ; ++i) {
      *out << "time_vec[" << i << "] = " << (*time_vec)[i] << std::endl;
    }
  }
}


template<class Scalar>
void InterpolationBuffer<Scalar>::removeNodes( Array<Scalar>& time_vec ) 
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  int N = time_vec.size();
#ifdef TEUCHOS_DEBUG
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
#endif // TEUCHOS_DEBUG
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
    for (unsigned int i=0; i<data_vec_->size() ; ++i) {
      out << "data_vec[" << i << "] = " << std::endl;
      (*data_vec_)[i].describe(out,this->getVerbLevel());
    }
  }
}


template <class Scalar>
void InterpolationBuffer<Scalar>::setParameterList(RCP<Teuchos::ParameterList> const& paramList)
{
  // 2007/12/05: rabartl: ToDo: Validate the parameter list!
  parameterList_ = paramList;
  int outputLevel = parameterList_->get( "outputLevel", int(-1) );
  outputLevel = std::min(std::max(outputLevel,-1),4);
  this->setVerbLevel(static_cast<Teuchos::EVerbosityLevel>(outputLevel));
  int policyLevel = parameterList_->get( "InterpolationBufferPolicy", int(1) );
  policyLevel = std::min(std::max(policyLevel,0),1);
  policy_ = static_cast<IBPolicy>(policyLevel);
  int storage_limit = parameterList_->get( "StorageLimit", storage_limit_ );
  if (storage_limit != storage_limit_) {
    this->setStorage(storage_limit);
  }
}


template <class Scalar>
RCP<Teuchos::ParameterList>
InterpolationBuffer<Scalar>::getNonconstParameterList()
{
  return(parameterList_);
}


template <class Scalar>
RCP<Teuchos::ParameterList> InterpolationBuffer<Scalar>::unsetParameterList()
{
  RCP<Teuchos::ParameterList> temp_param_list = parameterList_;
  parameterList_ = Teuchos::null;
  return(temp_param_list);
}


} // namespace Rythmos


#endif // Rythmos_INTERPOLATION_BUFFER_H
