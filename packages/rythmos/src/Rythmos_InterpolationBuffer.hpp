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
    InterpolationBuffer( int storage );

    SetStorage( int storage );
        
    /// Destructor
    ~InterpolationBuffer() {};

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
      ,const InterpolationBufferBase<Scalar>& IB);

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
    int order;

};

// ////////////////////////////
// Defintions
template<class Scalar>
InterpolationBuffer<Scalar>::InterpolationBuffer()
{
  order = 1;
  SetStorage(2);
}

template<class Scalar>
InterpolationBuffer<Scalar>::InterpolationBuffer( int storage_ )
{
  order = 1;
  SetStorage(storage_);
}

template<class Scalar>
InterpolationBuffer<Scalar>::SetStorage( int storage_ )
{
  storage_limit = min(2,storage_); // Minimum of two points so interpolation is possible
}

template<class Scalar>
bool InterpolationBuffer<Scalar>::SetPoints( 
    const std::vector<Scalar>& time_vec
    ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >& x_vec
    ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >& xdot_vec );
{
  std::list<DataStore<Scalar> > input_list;
  VectorToDataStoreList(&input_list,time_vec,x_vec,xdot_vec);
  input_list.sort();
  // Determine if time is already in list and if so, replace existing data with new data.
  std::list<DataStore<Scalar> >::iterator node_it;
  std::list<DataStore<Scalar> >::iterator input_it = input_list.begin();
  for (; input_it != input_list.end() ; input_it++)
  {
    node_it = find(node_list.begin(),node_list.end(),*input_it);
    if (node_it != node_list.end())
    {
      node_list.splice(node_it,input_list,input_it);
      node_list.erase(node_it);
    }
  }
  // Check that we're not going to exceed our storage limit:
  if ((node_list.size()+input_list.size()) > storage_limit)
    return(false);
  // Now add all the remaining points to node_list
  node_list.merge(input_list);
  return(true);
}

template<class Scalar>
bool InterpolationBuffer<Scalar>::GetPoints(
    const std::vector<Scalar>& time_vec
    ,std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >* x_vec
    ,std::vector<Tuechos::RefCountPtr<Thyra::VectorBase<Scalar> > >* xdot_vec
    ,std::vector<ScalarMag>* accuracy_vec) const
{

}

template<class Scalar>
bool InterpolationBuffer<Scalar>::GetPoints(
    const std::vector<Scalar>& time_vec
    ,std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >* x_vec
    ,std::vector<Tuechos::RefCountPtr<Thyra::VectorBase<Scalar> > >* xdot_vec
    ,std::vector<ScalarMag>* accuracy_vec) const
{
  // Copy the const time_vec to a local sorted time_vec
  std::vector<Scalar> local_time_vec = time_vec;
  std::sort(local_time_vec.begin(),local_time_vec.end());
  // If there are fewer than 2 points in node_list, then return failure
  if (node_list.size() < 2)
    return(false);
  // If time is outside range of t_values, then return failure
  if ( (*local_time_vec.begin() < node_list.begin().t) || (*local_time_vec.end() > node_list.end().t) )
    return(false);
  // Find t values on either side of time
  std::vector<Scalar>::iterator input_it = local_time_vec.begin();
  std::vector<DataStore<Scalar> >::iterator node_it = node_list.begin();
  for ( ; node_it != node_list.end() ; node_it++ )
  {
    while ((*input_it >= node_it->t) && (*input_it <= (node_it+1)->t))
    {
      Scalar& t = *input_it;
      Scalar& ti = node_it->t;
      Scalar& tip1 = (node_it+1)->t;
      Thyra::VectorBase<Scalar>& xi = *(node_it->x);
      Thyra::VectorBase<Scalar>& xip1 = *((node_it+1)->x);
      Thyra::VectorBase<Scalar>& xdoti = *(node_it->xdot);
      Thyra::VectorBase<Scalar>& xdotip1 = *((node_it+1)->xdot);

      // interpolate this point
      Scalar h = tip1-ti;
      // First we work on x.
      V_StVpStV(&*tmp_vec,Scalar(ST::one()/h),xi,Scalar(-ST::one()/h),xip1);
      Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > x = (node_it->x).clone_v();
      V_StVpStV(&*x, ST::one(), xi, t-ti, tmp_vec);
      x_vec.pushback(x);
      // Then we work on xdot.
      V_StVpStV(&*tmp_vec,Scalar(ST::one()/h),xdoti,Scalar(-ST::one()/h),xdotip1);
      Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > xdot = (node_it->xdot).clone_v();
      V_StVpStV(&*xdot, ST::one(), xdoti, t-ti, tmp_vec);
      xdot_vec.pushback(xdot);
      // And finally we estimate our order of accuracy
      accuracy_vec.pushback(h); 
      // Now we increment iterator for local_time_vec
      input_it++;
    }
  }
  return(true);
}

template<class Scalar>
bool InterpolationBuffer<Scalar>::SetRange(
    const Scalar& time_lower
    ,const Scalar& time_upper
    ,const InterpolationBufferBase<Scalar>& IB )
{
  std::vector<ScalarMag> input_nodes;
  bool status = IB.GetNodes(&input_nodes);
  if (status == false) return(status);
  std::sort(input_nodes.begin(),input_nodes.end());
  // Remove nodes outside the range [time_lower,time_upper]
  std::vector<ScalarMag>::iterator input_it = input_nodes.begin();
  for (; input_it != input_nodes.end() ; input_it++)
  {
    if (*input_it >= time_lower)
    {
      input_it--;
      break;
    }
  }
  input_nodes.erase(input_nodes.begin(),input_it);
  input_it = input_nodes.end();
  for (; input_it != input_nodes.begin() ; input_it--)
  {
    if (*input_it <= time_upper)
    {
      input_it++;
      break;
    }
  }
  input_nodes.erase(input_it,input_nodes.end());

  // Ask IB to interpolate more points if IB's order is higher than ours
  ScalarMag h_safety = ScalarMag(2*ST::one());
  int IBOrder = IB.GetOrder();
  if (IBOrder >= order)
  {
    input_it = input_nodes.begin();
    for (; input_it != input_nodes.end() ; input_it++)
    {
      std::vector<ScalarMag>::iterator input_it_next = input_it++;
      if (input_it_next == input_nodes.end())
        break;
      ScalarMag h_0 = *input_it_next - *input_it;
      ScalarMag h = h_0^(IBOrder/order)/h_safety;
      int N = ceil(h_0/h);
      h = ScalarMag(h_0/N);
      for (int i=1 ; i<N ; ++i)
      {
        input_nodes.insert(input_it_next,*input_it+i*h);
      }
    }
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

  std::vector<Teuchos::RefCountPtr<Scalar> > input_x;
  std::vector<Teuchos::RefCountPtr<Scalar> > input_xdot;
  std::vector<ScalarMag> input_accuracy;
  status = IB.GetPoints( input_nodes, &input_x, &input_xdot &input_accuracy );
  if (status == false) return(status);
  // We could check that the accuracy meets our criteria here.
  status = SetPoints( input_nodes, input_x, input_xdot );
  return(status);
}

template<class Scalar>
bool InterpolationBuffer<Scalar>::GetNodes( std::vector<Scalar>* time_list ) const
{
  std::list<DataStore<Scalar> >::iterator list_it = node_list.begin();
  for (; list_it != node_list.end() ; list_it++)
  {
    time_list->push_back(list_it->t);
  }
  return(true);
}

template<class Scalar>
bool InterpolationBuffer<Scalar>::RemoveNodes( std::vector<Scalar>& time_vec ) const
{
  int N = time_vec.size();
  for (int i=0; i<N ; ++i)
  {
    node_list.remove(time_vec[i]);
  }
  return(true);
}

template<class Scalar>
int InterpolationBuffer<Scalar>::GetOrder() const
{
  return(order);
}

} // namespace Rythmos

#endif // Rythmos_INTERPOLATION_BUFFER_H
