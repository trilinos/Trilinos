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
#include "Rythmos_LinearInterpolator.hpp"
#include "Rythmos_DataStore.hpp"

#include "Thyra_VectorBase.hpp"

namespace Rythmos {

/** \brief class for defining linear interpolation buffer functionality. */
template<class Scalar> 
class InterpolationBuffer : virtual public InterpolationBufferBase<Scalar>
{
  public:

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;
    
    /** \brief. */
    InterpolationBuffer();
    InterpolationBuffer( const Teuchos::RefCountPtr<InterpolatorBase<Scalar> >& interpolator_, int storage_ );

    /// Initialize the buffer:
    void initialize( const Teuchos::RefCountPtr<InterpolatorBase<Scalar> >& interpolator_, int storage_ );

    /// Set the interpolator for this buffer
    void SetInterpolator(const Teuchos::RefCountPtr<InterpolatorBase<Scalar> >& interpolator_);

    /// Set the maximum storage of this buffer
    void SetStorage( int storage );
        
    /// Destructor
    ~InterpolationBuffer() {};

    /// Add point to buffer
    bool SetPoints(
      const std::vector<Scalar>& time_vec
      ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >& x_vec
      ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >& xdot_vec
      ,const std::vector<ScalarMag> & accuracy_vec 
      );

    bool SetPoints(
      const std::vector<Scalar>& time_vec
      ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >& x_vec
      ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >& xdot_vec);

    /// Get value from buffer
    bool GetPoints(
      const std::vector<Scalar>& time_vec
      ,std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >* x_vec
      ,std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >* xdot_vec
      ,std::vector<ScalarMag>* accuracy_vec) const;

    /// Fill data in from another interpolation buffer
    bool SetRange(
      const Scalar& time_lower
      ,const Scalar& time_upper
      ,const InterpolationBufferBase<Scalar>& IB);

    /// Get interpolation nodes
    bool GetNodes(std::vector<Scalar>* time_vec) const;

    /// Remove interpolation nodes
    bool RemoveNodes(std::vector<Scalar>& time_vec);

    /// Get order of interpolation
    int GetOrder() const;

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
    void setParameterList(Teuchos::RefCountPtr<Teuchos::ParameterList> const& paramList);

    /** \brief . */
    Teuchos::RefCountPtr<Teuchos::ParameterList> getParameterList();

    /** \brief . */
    Teuchos::RefCountPtr<Teuchos::ParameterList> unsetParameterList();

    enum IBPolicy {
      BUFFER_STATIC = 0
      ,BUFFER_KEEP_NEWEST = 1
    };
    
  private:

    Teuchos::RefCountPtr<InterpolatorBase<Scalar> > interpolator;
    int storage_limit;
    typename DataStore<Scalar>::DataStoreVector_t data_vec;

    Teuchos::RefCountPtr<Teuchos::ParameterList> parameterList;

    IBPolicy policy;

};

// ////////////////////////////
// Defintions
template<class Scalar>
InterpolationBuffer<Scalar>::InterpolationBuffer()
{
  initialize(Teuchos::null,0);
}

template<class Scalar>
InterpolationBuffer<Scalar>::InterpolationBuffer( 
    const Teuchos::RefCountPtr<InterpolatorBase<Scalar> >& interpolator_
    ,int storage_ 
    )
{
  initialize(interpolator_,storage_);
}

template<class Scalar>
void InterpolationBuffer<Scalar>::initialize( 
    const Teuchos::RefCountPtr<InterpolatorBase<Scalar> >& interpolator_
    ,int storage_
    )
{
  Teuchos::RefCountPtr<Teuchos::FancyOStream> out = this->getOStream();
  out->precision(15);
  out->setMaxLenLinePrefix(30);
  out->pushLinePrefix("Rythmos::InterpolationBuffer");
  out->setShowLinePrefix(true);
  out->setTabIndentStr("    ");
  *out << "Initializing InterpolationBuffer" << std::endl;
  Teuchos::OSTab ostab(out,1,"IB::initialize");
  if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
    *out << "Calling SetInterpolator..." << std::endl;
  SetInterpolator(interpolator_);
  if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
    *out << "Calling SetStorage..." << std::endl;
  SetStorage(storage_);
  policy = BUFFER_KEEP_NEWEST;
}

template<class Scalar>
void InterpolationBuffer<Scalar>::SetStorage( int storage_ )
{
  storage_limit = max(2,storage_); // Minimum of two points so interpolation is possible
  Teuchos::RefCountPtr<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"IB::SetStorage");
  if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
    *out << "storage_limit = " << storage_limit << std::endl;
}

template<class Scalar>
void InterpolationBuffer<Scalar>::SetInterpolator(
    const Teuchos::RefCountPtr<InterpolatorBase<Scalar> >& interpolator_
    )
{
  if (interpolator_ == Teuchos::null)
  {
    interpolator = Teuchos::rcp(new LinearInterpolator<Scalar>);
  }
  else
  {
    interpolator = interpolator_;
  }
  Teuchos::RefCountPtr<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"IB::SetInterpolator");
  if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
    *out << "interpolator = " << interpolator->description() << std::endl;
}

template<class Scalar>
bool InterpolationBuffer<Scalar>::SetPoints( 
    const std::vector<Scalar>& time_vec
    ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >& x_vec
    ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >& xdot_vec 
    ,const std::vector<ScalarMag> & accuracy_vec 
    )
{
  Teuchos::RefCountPtr<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"IB::SetPoints");
  if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
  {
    *out << "time_vec = " << std::endl;
    for (unsigned int i=0 ; i<time_vec.size() ; ++i)
      *out << "time_vec[" << i << "] = " << time_vec[i] << std::endl;
    *out << "x_vec = " << std::endl;
    for (unsigned int i=0 ; i<x_vec.size() ; ++i)
    {
      *out << "x_vec[" << i << "] = " << std::endl;
      x_vec[i]->describe(*out,Teuchos::VERB_EXTREME);
    }
    if (xdot_vec.size() == 0)
      *out << "xdot_vec = empty vector" << std::endl;
    else
      *out << "xdot_vec = " << std::endl;
    for (unsigned int i=0 ; i<xdot_vec.size() ; ++i)
    {
      if (xdot_vec[i] == Teuchos::null)
        *out << "xdot_vec[" << i << "] = Teuchos::null" << std::endl;
      else
      {
        *out << "xdot_vec[" << i << "] = " << std::endl;
        xdot_vec[i]->describe(*out,Teuchos::VERB_EXTREME);
      }
    }
    if (accuracy_vec.size() == 0)
      *out << "accuracy_vec = empty vector" << std::endl;
    else
      *out << "accuracy_vec = " << std::endl;
    for (unsigned int i=0 ; i<accuracy_vec.size() ; ++i)
      *out << "accuracy_vec[" << i << "] = " << accuracy_vec[i] << std::endl;
  }
  typename DataStore<Scalar>::DataStoreList_t input_data_list;
  VectorToDataStoreList<Scalar>(time_vec,x_vec,xdot_vec,accuracy_vec,&input_data_list);
  input_data_list.sort();
  if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
  {
    *out << "input_data_list after sorting = " << std::endl;
    typename DataStore<Scalar>::DataStoreList_t::iterator
      data_it = input_data_list.begin();
    int i=0;
    for (; data_it != input_data_list.end() ; data_it++)
    {
      *out << "item " << i << ":" << std::endl;
      data_it->describe(*out,Teuchos::VERB_EXTREME);
      i++;
    }
  }
  // Determine if time is already in list and if so, replace existing data with new data.
  typename DataStore<Scalar>::DataStoreList_t::iterator 
    input_it = input_data_list.begin();
  while (input_it != input_data_list.end())
  {
    if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
      *out << "Looking for time = " << (*input_it).time << " in data_vec to replace with new value" << std::endl;
    typename DataStore<Scalar>::DataStoreVector_t::iterator 
      node_it = std::find(data_vec.begin(),data_vec.end(),*input_it);
    if (node_it != data_vec.end())
    {
      int node_index = node_it - data_vec.begin(); // 10/17/06 tscoffe:  this
                                                   // is how you back out an
                                                   // element's index into a
                                                   // vector from its iterator.
      if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
      {
        *out << "Replacing data_vec[" << node_index << "] = " << std::endl;
        node_it->describe(*out,Teuchos::VERB_EXTREME);
        *out << "with:" << std::endl;
        input_it->describe(*out,Teuchos::VERB_EXTREME);
      }
      data_vec[node_index] = *input_it;
      input_it = input_data_list.erase(input_it);
      if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
      {
        *out << "input_data_list after removing an element = " << std::endl;
        typename DataStore<Scalar>::DataStoreList_t::iterator
          data_it = input_data_list.begin();
        int i=0;
        for (; data_it != input_data_list.end() ; data_it++)
        {
          *out << "item " << i << ":" << std::endl;
          data_it->describe(*out,Teuchos::VERB_EXTREME);
          i++;
        }
      }
    }
    else
    {
      input_it++;
    }
  }
  // Check that we're not going to exceed our storage limit:
  if ((data_vec.size()+input_data_list.size()) > static_cast<unsigned int>(storage_limit))
  { 
    if (policy == BUFFER_STATIC)
      return(false);
    else if (policy == BUFFER_KEEP_NEWEST)
    {
      if (input_data_list.front() > data_vec.back())
      {
        // Case:  all of new points are past end of existing points
        // Remove points from the beginning of data_vec, then add new points
        int num_extra_points = input_data_list.size();
        typename DataStore<Scalar>::DataStoreVector_t::iterator 
          data_it = data_vec.begin();
        for (int i=0 ; i < num_extra_points ; ++i)
          data_it++;
        if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
        {
          *out << "Removing " << num_extra_points 
            << " from beginning of data_vec to make room for new points." << std::endl;
        }
        data_vec.erase(data_vec.begin(),data_it);
      }
      else
      {
        // Case:  At least one new point is before end of existing points
        return(false);
      }
    }
    else
      return(false);
  }
  // Now add all the remaining points to data_vec
  data_vec.insert(data_vec.end(),input_data_list.begin(),input_data_list.end());
  // And sort data_vec:
  std::sort(data_vec.begin(),data_vec.end());
  if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
  {
    *out << "data_vec at end of SetPoints:" << std::endl;
    for (unsigned int i=0 ; i<data_vec.size() ; ++i)
    {
      *out << "data_vec[" << i << "] = " << std::endl;
      data_vec[i].describe(*out,Teuchos::VERB_EXTREME);
    }
  }
  return(true);
}

template<class Scalar>
bool InterpolationBuffer<Scalar>::SetPoints( 
    const std::vector<Scalar>& time_vec
    ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >& x_vec
    ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >& xdot_vec )
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  std::vector<ScalarMag> accuracy_vec;
  accuracy_vec.reserve(x_vec.size());
  for (unsigned int i=0;i<x_vec.size();++i)
    accuracy_vec[i] = ST::zero();
  return(this->SetPoints(time_vec,x_vec,xdot_vec,accuracy_vec));
}

template<class Scalar>
bool InterpolationBuffer<Scalar>::GetPoints(
    const std::vector<Scalar>& time_vec
    ,std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >* x_vec
    ,std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >* xdot_vec
    ,std::vector<ScalarMag>* accuracy_vec) const
{
  Teuchos::RefCountPtr<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"IB::GetPoints");
  if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
    *out << "Calling interpolate..." << std::endl;
  typename DataStore<Scalar>::DataStoreVector_t data_out;
  bool status = interpolator->interpolate(data_vec, time_vec, &data_out);
  if (!status) return(status);
  if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
    *out << "Interpolation successful" << std::endl;
  std::vector<Scalar> time_out_vec;
  DataStoreVectorToVector<Scalar>(data_out, &time_out_vec, x_vec, xdot_vec, accuracy_vec);
  // Double check that time_out_vec == time_vec
  if (time_vec.size() != time_out_vec.size()) return(false);
  if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
    *out << "Conversion of DataStoreVector to Vector successful" << std::endl;
  return(true);
}


template<class Scalar>
bool InterpolationBuffer<Scalar>::SetRange(
    const Scalar& time_lower
    ,const Scalar& time_upper
    ,const InterpolationBufferBase<Scalar>& IB )
{
  Teuchos::RefCountPtr<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"IB::SetRange");
  if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
  {
    *out << "time_lower = " << time_lower << std::endl;
    *out << "time_upper = " << time_upper << std::endl;
    *out << "IB = " << IB.description() << std::endl;
  }
  std::vector<Scalar> input_nodes;
  bool status = IB.GetNodes(&input_nodes);
  if (!status) return(status);
  std::sort(input_nodes.begin(),input_nodes.end());
  if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
  {
    *out << "input_nodes after sorting = " << std::endl;
    for (unsigned int i=0 ; i<input_nodes.size() ; ++i)
      *out << "input_nodes[" << i << "] = " << input_nodes[i] << std::endl;
  }
  // Remove nodes outside the range [time_lower,time_upper]
  typename std::vector<Scalar>::iterator input_it_lower = input_nodes.begin();
  for (; input_it_lower != input_nodes.end() ; input_it_lower++)
  {
    if (*input_it_lower >= time_lower)
    {
      break;
    }
  }
  if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
  {
    int n0 = 0;
    int n1 = input_it_lower - input_nodes.begin();
    *out << "Removing input_nodes before time_lower with indices: [" << n0 << "," << n1 << ")" << std::endl;
    for (int i=n0 ; i<n1; ++i)
    {
      *out << "  input_nodes[" << i << "] = " << input_nodes[i] << std::endl;
    }
  }
  // tscoffe 10/19/06 Note:  erase removes the range [it_begin,it_end)
  if (input_it_lower - input_nodes.begin() >= 0)
    input_nodes.erase(input_nodes.begin(),input_it_lower);
  typename std::vector<Scalar>::iterator input_it_upper = input_nodes.end();
  input_it_upper--;
  for (; input_it_upper != input_nodes.begin() ; input_it_upper--)
  {
    if (*input_it_upper <= time_upper)
    {
      input_it_upper++;
      break;
    }
  }
  // This is to handle the case of one element in the vector
  if (*input_it_upper <= time_upper)
    input_it_upper++;
  if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
  {
    int n0 = input_it_upper - input_nodes.begin();
    int n1 = input_nodes.size();
    *out << "Removing input_nodes after time_upper with indices [" << n0 << "," << n1 << ")" << std::endl;
    for (int i=n0 ; i<n1; ++i)
      *out << "  input_nodes[" << i << "] = " << input_nodes[i] << std::endl;
  }
  if (static_cast<unsigned int>(input_it_upper - input_nodes.begin()) < input_nodes.size())
    input_nodes.erase(input_it_upper,input_nodes.end());
  if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
  {
    *out << "input_nodes remaining:" << std::endl;
    for (unsigned int i=0 ; i<input_nodes.size() ; ++i)
      *out << "input_nodes[" << i << "] = " << input_nodes[i] << std::endl;
  }

  // Ask IB to interpolate more points if IB's order is higher than ours
  typedef Teuchos::ScalarTraits<Scalar> ST;
  Scalar h_safety = Scalar(2*ST::one());
  int IBOrder = IB.GetOrder();
  if (IBOrder >= interpolator->order())
  {
    std::list<Scalar> add_nodes;
    for (unsigned int i=0 ; i<input_nodes.size()-1 ; ++i)
    {
      Scalar h_0 = input_nodes[i+1] - input_nodes[i];
      Scalar h = pow(h_0,(IBOrder/interpolator->order())/h_safety);
      if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
      {
        *out << "i = " << i << std::endl;
        *out << "interpolator->order() = " << interpolator->order() << std::endl;
        *out << "IB.GetOrder() = " << IB.GetOrder() << std::endl;
        *out << "h = " << h << std::endl;
      }
      Scalar N = ceil(h_0/h);
      h = Scalar(h_0/N);
      if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
      {
        *out << "h_0 = " << h_0 << std::endl;
        *out << "N = " << N << std::endl;
        *out << "h = " << h << std::endl;
        *out << "Inserting an additional " << N-1 << " points to be interpolated:" << std::endl;
      }
      for (int j=1 ; j<N ; ++j)
      {
        if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
          *out << input_nodes[i]+j*h << std::endl;
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
  // Use SetPoints and check return value to make sure we observe storage_limit.

  std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > > input_x;
  std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > > input_xdot;
  std::vector<ScalarMag> input_accuracy;
  status = IB.GetPoints( input_nodes, &input_x, &input_xdot, &input_accuracy );
  if (!status) return(status);
  // We could check that the accuracy meets our criteria here.
  status = SetPoints( input_nodes, input_x, input_xdot, input_accuracy );
  return(status);
}

template<class Scalar>
bool InterpolationBuffer<Scalar>::GetNodes( std::vector<Scalar>* time_vec ) const
{
  int N = data_vec.size();
  time_vec->clear();
  time_vec->reserve(N);
  for (int i=0 ; i<N ; ++i)
    time_vec->push_back(data_vec[i].time);
  Teuchos::RefCountPtr<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"IB::GetNodes");
  if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
  {
    *out << this->description() << std::endl;
    for (unsigned int i=0 ; i<time_vec->size() ; ++i)
      *out << "time_vec[" << i << "] = " << (*time_vec)[i] << std::endl;
  }
  return(true);
}

template<class Scalar>
bool InterpolationBuffer<Scalar>::RemoveNodes( std::vector<Scalar>& time_vec ) 
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  ScalarMag z = ST::zero();
  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > vec_temp;
  int N = time_vec.size();
  for (int i=0; i<N ; ++i)
  {
    DataStore<Scalar> ds_temp(time_vec[i],vec_temp,vec_temp,z);
    typename DataStore<Scalar>::DataStoreVector_t::iterator 
      data_it = std::find(data_vec.begin(),data_vec.end(),ds_temp);
    if (data_it != data_vec.end())
      data_vec.erase(data_it);
  }
  return(true);
}

template<class Scalar>
int InterpolationBuffer<Scalar>::GetOrder() const
{
  return(interpolator->order());
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
  if ( (static_cast<int>(verbLevel) == static_cast<int>(Teuchos::VERB_DEFAULT) ) ||
       (static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW)     )
     )
  {
    out << description() << "::describe" << std::endl;
    out << "interpolator = " << interpolator->description() << std::endl;
    out << "storage_limit = " << storage_limit << std::endl;
  }
  else if (static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
  {
  }
  else if (static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_MEDIUM))
  {
  }
  else if (static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_HIGH))
  {
    out << "data_vec = " << std::endl;
    for (unsigned int i=0; i<data_vec.size() ; ++i)
    {
      out << "data_vec[" << i << "] = " << std::endl;
      data_vec[i].describe(out,this->getVerbLevel());
    }
  }
}

template <class Scalar>
void InterpolationBuffer<Scalar>::setParameterList(Teuchos::RefCountPtr<Teuchos::ParameterList> const& paramList)
{
  parameterList = paramList;
  int outputLevel = parameterList->get( "outputLevel", int(-1) );
  outputLevel = min(max(outputLevel,-1),4);
  this->setVerbLevel(static_cast<Teuchos::EVerbosityLevel>(outputLevel));
  int policyLevel = parameterList->get( "InterpolationBufferPolicy", int(1) );
  policyLevel = min(max(policyLevel,0),1);
  policy = static_cast<IBPolicy>(policyLevel);
}

template <class Scalar>
Teuchos::RefCountPtr<Teuchos::ParameterList> InterpolationBuffer<Scalar>::getParameterList()
{
  return(parameterList);
}

template <class Scalar>
Teuchos::RefCountPtr<Teuchos::ParameterList> InterpolationBuffer<Scalar>::unsetParameterList()
{
  Teuchos::RefCountPtr<Teuchos::ParameterList> temp_param_list = parameterList;
  parameterList = Teuchos::null;
  return(temp_param_list);
}

} // namespace Rythmos

#endif // Rythmos_INTERPOLATION_BUFFER_H
