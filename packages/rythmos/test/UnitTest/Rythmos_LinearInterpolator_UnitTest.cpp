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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER

#include "Teuchos_UnitTestHarness.hpp"

#include "Rythmos_InterpolatorBaseHelpers.hpp"
#include "Rythmos_LinearInterpolator.hpp"
#include "Rythmos_UnitTestHelpers.hpp"

#include "Teuchos_Polynomial.hpp"
#include "Thyra_DetachedVectorView.hpp"

namespace Rythmos {

using Teuchos::Polynomial;
using Teuchos::RCP;

TEUCHOS_UNIT_TEST( Rythmos_LinearInterpolator, nonMemberConstructor ) {
  RCP<LinearInterpolator<double> > li = linearInterpolator<double>();
  TEST_EQUALITY_CONST( li->order(), 1 );
}

TEUCHOS_UNIT_TEST( Rythmos_LinearInterpolator, interpolate ) {
  RCP<InterpolatorBase<double> > li = linearInterpolator<double>();
  double maxOrder = li->order();
  for (int order = 0 ; order <= maxOrder+1 ; ++order) {
    TEST_EQUALITY( order, order );
    Polynomial<double> p(order,1.0); 
    RCP<DataStore<double>::DataStoreVector_t> data_in = rcp( new DataStore<double>::DataStoreVector_t );;
    double T0 = 0.0;
    double T1 = 1.0;
    int N = 5;
    for (int i=0 ; i < N ; ++i) {
      double t = ((T1-T0)/(N-1.0))*i+T0;
      double x = 0.0;
      p.evaluate(t,&x);
      RCP<Thyra::VectorBase<double> > xv, xvdot;
      double accuracy = 0.0;
      xv = createDefaultVector<double>(1,x);
      data_in->push_back(DataStore<double>(t,xv,xvdot,accuracy));
    }
    Array<double> t_values;
    DataStore<double>::DataStoreVector_t data_out;
    N = 2*N;
    for (int i=0 ; i < N ; ++i) {
      double t = ((T1-T0)/(N-1.0))*i+T0;
      t_values.push_back(t);
    }
    interpolate<double>(*li, data_in, t_values, &data_out);
    // Verify that the interpolated values are exactly the same as the polynomial values
    // unless the order of polynomial is greater than the order of the interpolator
    unsigned int N_out = data_out.size();
    for (unsigned int i=0 ; i < N_out ; ++i ) {
      double x = 0.0;
      double t = data_out[i].time;
      RCP<const Thyra::VectorBase<double> > xv = data_out[i].x;
      {
        Thyra::ConstDetachedVectorView<double> xv_view(*xv);
        x = xv_view[0];
      }
      double x_exact = 0.0;
      p.evaluate(t,&x_exact);
      double tol = 1.0e-15;
      if ((order <= maxOrder) || (i == 0) || (i == N_out-1)) {
        TEST_FLOATING_EQUALITY( x, x_exact, tol );
      } else {
        TEST_COMPARE( fabs((x-x_exact)/x_exact), >, tol );
      }
    }
  }
}

TEUCHOS_UNIT_TEST( Rythmos_LinearInterpolator, boundary_interpolate ) {
  RCP<DataStore<double>::DataStoreVector_t> data_in = rcp( new DataStore<double>::DataStoreVector_t );;
  int N = 2;
  for (int i=0 ; i<N ; ++i) {
    double time = 0.0 + 1.0*i;
    RCP<VectorBase<double> > x = createDefaultVector(2,1.0+1.0*i);
    RCP<VectorBase<double> > xdot = createDefaultVector(2,10.0+1.0*i);
    double accuracy = 0.0;
    data_in->push_back(DataStore<double>(time,x,xdot,accuracy));
  }
  RCP<LinearInterpolator<double> > li = linearInterpolator<double>();
  Array<double> t_values;
  t_values.push_back(0.0);
  t_values.push_back(1.0);
  DataStore<double>::DataStoreVector_t data_out;
  interpolate<double>(*li, data_in, t_values, &data_out);
  for (int i=0 ; i<N ; ++i) {
    TEST_EQUALITY( (*data_in)[i].time, data_out[i].time );
    TEST_EQUALITY( (*data_in)[i].x.get(), data_out[i].x.get() );
    TEST_EQUALITY( (*data_in)[i].xdot.get(), data_out[i].xdot.get() );
  }
}
TEUCHOS_UNIT_TEST( Rythmos_LinearInterpolator, interpolate_then_clone ) {
  // 04/20/09 tscoffe:  Based on the interpolator use in StepperHelpers:defaultGetPoints(...)
  double t_old = 1.0;
  RCP<VectorBase<double> > x_old = createDefaultVector(1,2.0);
  RCP<VectorBase<double> > xdot_old = createDefaultVector(1,3.0);
  double t = 2.0;
  RCP<VectorBase<double> > x = createDefaultVector(1,4.0);
  RCP<VectorBase<double> > xdot = createDefaultVector(1,5.0);
  DataStore<double>::DataStoreVector_t ds_nodes;
  DataStore<double>::DataStoreVector_t ds_out;
  {
    // t_old
    DataStore<double> ds;
    ds.time = t_old;
    ds.x = rcp(x_old.get(),false);
    ds.xdot = rcp(xdot_old.get(),false);
    ds_nodes.push_back(ds);
  }
  {
    // t
    DataStore<double> ds;
    ds.time = t;
    ds.x = rcp(x.get(),false);
    ds.xdot = rcp(xdot.get(),false);
    ds_nodes.push_back(ds);
  }
  Array<double> time_vec_in;
  time_vec_in.push_back(1.5);
  RCP<InterpolatorBase<double> > interpolator = linearInterpolator<double>();
  interpolate<double>(*interpolator,rcp(&ds_nodes,false),time_vec_in,&ds_out);
  Array<double> time_vec_out;
  Array<RCP<const VectorBase<double> > > x_vec_out;
  Array<RCP<const VectorBase<double> > > xdot_vec_out;
  Array<Teuchos::ScalarTraits<double>::magnitudeType> accuracy_vec_out;
  dataStoreVectorToVector(ds_out,&time_vec_out,&x_vec_out,&xdot_vec_out,&accuracy_vec_out);
  TEST_ASSERT( time_vec_out.length()==1 );
  RCP<VectorBase<double> > tmpVec = x_vec_out[0]->clone_v();
  RCP<VectorBase<double> > tmpVecDot = xdot_vec_out[0]->clone_v();
  {
    Thyra::ConstDetachedVectorView<double> x_view( tmpVec );
    Thyra::ConstDetachedVectorView<double> xdot_view( tmpVecDot );
    TEST_EQUALITY_CONST( x_view[0], 3.0 );
    TEST_EQUALITY_CONST( xdot_view[0], 4.0 );
    TEST_EQUALITY_CONST( accuracy_vec_out[0], 1.0 );
  }
}

} // namespace Rythmos



