//@HEADER
// ***********************************************************************
//
//                     Rythmos Package
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

#ifndef Rythmos_INTERPOLATOR_BASE_HELPERS_H
#define Rythmos_INTERPOLATOR_BASE_HELPERS_H

#include "Rythmos_InterpolatorBase.hpp"

namespace Rythmos {

/** \relates InterplatorBase . */
template<class Scalar>
void interpolate(
    InterpolatorBase<Scalar>& interp,
    const RCP<const typename DataStore<Scalar>::DataStoreVector_t> & nodes,
    const Array<Scalar> &t_values,
    typename DataStore<Scalar>::DataStoreVector_t *data_out
    )
{
  interp.setNodes(nodes);
  interp.interpolate(t_values,data_out);
}

/** \relates InterpolatorBase . */
template<class Scalar>
void assertBaseInterpolatePreconditions(
  const typename DataStore<Scalar>::DataStoreVector_t &data_in,
  const Array<Scalar> &t_values,
  typename DataStore<Scalar>::DataStoreVector_t *data_out
  )
{
  TEST_FOR_EXCEPTION(
      data_in.size()==0, std::logic_error,
      "Error, data_in.size() == 0!\n"
      );
  Array<Scalar> time_vec;
  dataStoreVectorToVector<Scalar>(data_in, &time_vec, 0, 0, 0);
  assertTimePointsAreSorted<Scalar>(time_vec);
  assertTimePointsAreSorted<Scalar>(t_values);
  if (data_in.size() == 1) {
    TEST_FOR_EXCEPTION(
      t_values.size()>1, std::logic_error,
      "Error, data_in.size() == 1, but t_values.size() > 1!\n"
      );
    TEST_FOR_EXCEPTION(
      t_values[0]!=data_in[0].time, std::logic_error,
      "Error, data_in.size) == 1, but t_values[0] = " << 
      t_values[0] << " != " << data_in[0].time << " = data_in[0].time!\n"
      );
  }
  TimeRange<Scalar> range(data_in.front().time,data_in.back().time);
  for (int i=0; i<Teuchos::as<int>(t_values.size()) ; ++i) {
    TEST_FOR_EXCEPTION(
      !range.isInRange(t_values[i]), std::out_of_range,
      "Error, t_values[" << i << "] = " << t_values[i] << 
      " is not in range of data_in = " << range << "!\n"
      );
  }
  TEST_FOR_EXCEPTION(
    data_out == 0, std::logic_error,
    "Error, data_out = NULL!\n"
    );
  for (int i=0; i<Teuchos::as<int>(data_in.size()) ; ++i) {
    TEST_FOR_EXCEPTION(
      data_in[i].x == Teuchos::null, std::logic_error,
      "Error, data_in[" << i << "].x == Teuchos::null.\n"
      );
  }
}

template<class Scalar>
void assertNodesUnChanged(
    const typename DataStore<Scalar>::DataStoreVector_t & nodes, 
    const typename DataStore<Scalar>::DataStoreVector_t & nodes_copy 
    ) 
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  int N = nodes.size();
  int Ncopy = nodes_copy.size();
  TEST_FOR_EXCEPTION( N != Ncopy, std::logic_error, 
      "Error!  The number of nodes passed in through setNodes has changed!"
      );
  if (N > 0) {
    RCP<Thyra::VectorBase<Scalar> > xdiff = nodes[0].x->clone_v();
    RCP<Thyra::VectorBase<Scalar> > xdotdiff = xdiff->clone_v();
    V_S(outArg(*xdiff),ST::one());
    V_S(outArg(*xdotdiff),ST::one());
    for (int i=0 ; i<N ; ++i) {
      V_StVpStV(outArg(*xdiff),ST::one(),*nodes[i].x,-ST::one(),*nodes_copy[i].x);
      if ((!Teuchos::is_null(nodes[i].xdot)) && (!Teuchos::is_null(nodes_copy[i].xdot))) {
        V_StVpStV(outArg(*xdotdiff),ST::one(),*nodes[i].xdot,-ST::one(),*nodes_copy[i].xdot);
      } else if (Teuchos::is_null(nodes[i].xdot) && Teuchos::is_null(nodes_copy[i].xdot)) {
        V_S(outArg(*xdotdiff),ST::zero());
      }
      Scalar xdiffnorm = norm_inf(*xdiff);
      Scalar xdotdiffnorm = norm_inf(*xdotdiff);
      TEST_FOR_EXCEPTION(
          ( ( nodes[i].time != nodes_copy[i].time ) ||
            ( xdiffnorm != ST::zero() ) ||
            ( xdotdiffnorm != ST::zero() ) ||
            ( nodes[i].accuracy != nodes_copy[i].accuracy ) ), 
          std::logic_error,
          "Error!  The data in the nodes passed through setNodes has changed!"
          );
    }
  }
}


} // namespace Rythmos 

#endif // Rythmos_INTERPOLATOR_BASE_HELPERS_H


