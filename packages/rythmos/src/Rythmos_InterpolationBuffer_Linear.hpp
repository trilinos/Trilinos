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

#ifndef Rythmos_INTERPOLATION_BUFFER_LINEAR_H
#define Rythmos_INTERPOLATION_BUFFER_LINEAR_H

#include "Rythmos_InterpolationBuffer.hpp"
#include "Thyra_VectorBase.hpp"

namespace Rythmos {

/** \brief class for defining linear interpolation buffer functionality. */
template<class Scalar> 
class LinearInterpolationBuffer : virtual public InterpolationBuffer<Scalar>
{
  public:

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;
    
    /** \brief. */
    LinearInterpolationBuffer();
    LinearInterpolationBuffer( int storage );

    SetStorage( int storage );
        
    /// Destructor
    ~LinearInterpolationBuffer() {};

    /// Add point to buffer
    bool SetPoints(
      const std::vector<Scalar>& time_list
      ,const std::vector<Thyra::VectorBase<Scalar> >& x_list
      ,const std::vector<Thyra::VectorBase<Scalar> >& xdot_list);

    /// Get value from buffer
    bool GetPoints(
      const std::vector<Scalar>& time_list
      ,std::vector<Thyra::VectorBase<Scalar> >* x_list
      ,std::vector<Thyra::VectorBase<Scalar> >* xdot_list
      ,std::vector<ScalarMag>* accuracy_list) const;

    /// Fill data in from another interpolation buffer
    bool SetRange(
      const Scalar& time_lower
      ,const Scalar& time_upper
      ,const InterpolationBuffer<Scalar>& IB);

    /// Get interpolation nodes
    bool GetNodes(std::vector<Scalar>* time_list) const;

    /// Remove interpolation nodes
    bool RemoveNodes(std::vector<Scalar>* time_list) const;

    /// Get order of interpolation
    int GetOrder() const;

  private:

    int storage_limit;
    Thyra::VectorBase<Scalar> tmp_vec;
    std::list<Teuchos::RefCountPtr<DataStore<Scalar> > > node_list;

    template<class Scalar>
    class DataStore
    {
      public:
        ~DataStore();
        DataStore();
        DataStore(ScalarMag &time
          ,Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > &x
          ,Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > &xdot);
        ScalarMag time;
        Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > x;
        Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > xdot;
        bool operator< (const DataStore<Scalar>& d1, const DataStore<Scalar>& d2) const
        { return( d1.time < d2.time ); }
    }
    void VectorToDataStoreList(
      std::list<Teuchos::RefCountPtr<DataStore<Scalar> > > *list_ds
      ,const std::vector<ScalarMag> time_vec
      ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > > x_vec
      ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > > xdot_vec) const;
    void DataStoreListToVector(
      std::vector<ScalarMag> *time_vec
      ,std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > > *x_vec
      ,std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > > *xdot_vec
      ,const std::list<Teuchos::RefCountPtr<DataStore<Scalar> > > &list_ds) const;


};

// ////////////////////////////
// Defintions

template<class Scalar>
DataStore<Scalar>::DataStore(ScalarMag &time_
  ,Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > &x_
  ,Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > &xdot_)
{
  time = time_;
  x = x_;
  xdot = xdot_;
}

template<class Scalar>
void LinearInterpolationBuffer<Scalar>::VectorToDataStoreList(
      std::list<Teuchos::RefCountPtr<DataStore<Scalar> > > *list_ds
      ,const std::vector<ScalarMag> time_vec
      ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > > x_vec
      ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > > xdot_vec) const
{
  int N = time_vec.size();
  int Nx = x_vec.size();
  int Nxdot = xdot_vec.size();
  if ((N != Nx) || (N != Nxdot))
  {
    list_ds = NULL;
    return;
  }
  list_ds->clear();
  for (int i=0; i<N ; ++i)
  {
    list_ds->push_back(Teuchos::rcp(new DataStore(time_list[i],x_list[i],xdot_list[i])));
  }
}

template<class Scalar>
void LinearInterpolationBuffer<Scalar>::DataStoreListToVector(
      std::vector<ScalarMag> *time_vec
      ,std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > > *x_vec
      ,std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > > *xdot_vec
      ,const std::list<Teuchos::RefCountPtr<DataStore<Scalar> > > &list_ds) const
{
  int N = list_ds.size();
  time_list->reserve(N); time_list->clear();
  x_list->reserve(N);    x_list->clear();
  xdot_list->reserve(N); xdot_list->clear();
  std::list<DataStore<Scalar> >::iterator list_it = list_ds.begin();
  for (; list_it < list_ds.end() ; list_it++)
  {
    time_list->push_back(list_it->time);
    x_list->push_back(list_it->x);
    x_dot_list->push_back(list_it->xdot);
  }
}

template<class Scalar>
LinearInterpolationBuffer<Scalar>::LinearInterpolationBuffer()
{
  SetStorage(2);
}

template<class Scalar>
LinearInterpolationBuffer<Scalar>::LinearInterpolationBuffer( int storage_ )
{
  SetStorage(storage_);
}

template<class Scalar>
LinearInterpolationBuffer<Scalar>::SetStorage( int storage_ )
{
  storage_limit = min(2,storage_); // Minimum of two points so interpolation is possible
}

template<class Scalar>
bool LinearInterpolationBuffer<Scalar>::SetPoints( 
    const std::vector<Scalar>& time_list
    ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >& x_list
    ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >& xdot_list );
{
  // Sanity checking:
  if ((time_list.size() != x_list.size()) || (time_list.size() != xdot_list.size()))
    return(false);
  if ((t_values.size()+time_list.size()) > storage_limit)
    return(false);
  // Determine if time is already in list and if so, replace existing data with new data.
  num_nodes = t_values.size();
  for (int i=0 ; i < num_nodes ; ++i)
  {
    std::vector<Scalar>::iterator it_dup = std::find(t_values.begin(),t_values.end(),t_values[i]);
    if (it_dup != t_values.last())
    {
       // int idx = "index of it_dup into time_list"
       // Replace x_values[i] with x_list[idx]
       // Replace xdot_values[i] with xdot_list[idx]
       // Remove idx from list of time_values to insert into buffer
    }
  }
  // Determine where time should be in list and insert it along with x and xdot.
  return(false);
}

template<class Scalar>
bool LinearInterpolationBuffer<Scalar>::GetPoints(
    const std::vector<Scalar>& time_list
    ,std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >* x_list
    ,std::vector<Tuechos::RefCountPtr<Thyra::VectorBase<Scalar> > >* xdot_list
    ,std::vector<ScalarMag>* accuracy_list) const
{
  // Copy the const time_list to a local sorted time_list
  std::vector<Scalar> local_time_list = time_list;
  std::sort(local_time_list.begin(),local_time_list.end());
  // If time is outside range of t_values, then return failure
  if ( (local_time_list.begin() < t_values.begin()) || (local_time_list.end() > t_values.end()) )
  {
    return(false);
  }
  // Find t_values on either side of time
  int num_nodes = t_values.size();
  int num_time_list = time_list.size();
  std::vector<Scalar>::iterator it_t = local_time_list.begin();
  for (i=0; i < num_nodes-1 ; ++i)
  {
    while ((*it_t >= t_values[i]) && (*it_t <= t_values[i+1]))
    {
      // interpolate this point
      Scalar h = t_values[i+1]-t_values[i];
      // First we work on x.
      V_StVpStV(&*tmp_vec,Scalar(ST::one()/h),x_values[i],Scalar(-ST::one()/h),x_values[i+1]);
      Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > x = x_values[i].clone_v();
      V_StVpStV(&*x, ST::one(), x_values[i], (*it_t)-t_values[i], tmp_vec);
      x_list.pushback(x);
      // Then we work on xdot.
      V_StVpStV(&*tmp_vec,Scalar(ST::one()/h),xdot_values[i],Scalar(-ST::one()/h),xdot_values[i+1]);
      Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > xdot = xdot_values[i].clone_v();
      V_StVpStV(&*xdot, ST::one(), xdot_values[i], (*it_t)-t_values[i], tmp_vec);
      xdot_list.pushback(xdot);
      // And finally we estimate our order of accuracy
      accuracy_list.pushback(h); 
      // Now we increment iterator for local_time_list
      it_t++;
    }
  }
  return(true);
}

template<class Scalar>
bool LinearInterpolationBuffer<Scalar>::SetRange(
    const Scalar& time_lower
    ,const Scalar& time_upper
    ,const InterpolationBuffer<Scalar>& IB )
{
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
  return(false);
}

template<class Scalar>
bool LinearInterpolationBuffer<Scalar>::GetNodes( std::vector<Scalar>* time_list ) const
{
  time_list = t_values; // std::copy of data (this may be incorrect)
  return(true);
}
template<class Scalar>
bool LinearInterpolationBuffer<Scalar>::RemoveNodes( std::vector<Scalar>& time_list ) const
{
  // Search through t_values and when you find one that matches an item in
  // time_list, then remove the time and the x, x_dot values 
  return(false);
}

template<class Scalar>
int LinearInterpolationBuffer<Scalar>::GetOrder() const
{
  return(1);
}

} // namespace Rythmos

#endif //Rythmos_INTERPOLATION_BUFFER_LINEAR_H
