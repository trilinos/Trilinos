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

#include "Teuchos_UnitTestHarness.hpp"

#include "Rythmos_Types.hpp"
#include "Rythmos_InterpolatorBaseHelpers.hpp"
#include "Rythmos_CubicSplineInterpolator.hpp"
#include "Rythmos_UnitTestHelpers.hpp"

#include "Teuchos_Polynomial.hpp"
#include "Thyra_DetachedVectorView.hpp"

namespace Rythmos {

using Teuchos::Polynomial;
using Teuchos::outArg;
using Teuchos::rcp;

TEUCHOS_UNIT_TEST( Rythmos_CubicSplineInterpolator, nonMemberConstructor ) {
  RCP<CubicSplineInterpolator<double> > csi = cubicSplineInterpolator<double>();
  TEST_EQUALITY_CONST( csi->order(), 1 );
}

TEUCHOS_UNIT_TEST( Rythmos_CubicSplineInterpolator, evaluateCubicSpline ) {
  CubicSplineCoeff<double> coeff;
  coeff.t.push_back(1.0);
  coeff.t.push_back(2.0);
  coeff.a.push_back(createDefaultVector(1,1.0));
  coeff.b.push_back(createDefaultVector(1,2.0));
  coeff.c.push_back(createDefaultVector(1,3.0));
  coeff.d.push_back(createDefaultVector(1,4.0));

  int index = 0;
  double t = 1.5;
  RCP<Thyra::VectorBase<double> > S   = createDefaultVector(1,0.0);
  RCP<Thyra::VectorBase<double> > Sp  = createDefaultVector(1,0.0);
  RCP<Thyra::VectorBase<double> > Spp = createDefaultVector(1,0.0);
  double tol = 1.0e-10;
  evaluateCubicSpline<double>(coeff, index, t, outArg(*S), outArg(*Sp), outArg(*Spp));
  {
    Thyra::ConstDetachedVectorView<double> S_view(*S);
    Thyra::ConstDetachedVectorView<double> Sp_view(*Sp);
    Thyra::ConstDetachedVectorView<double> Spp_view(*Spp);
    TEST_FLOATING_EQUALITY( S_view[0],   1.0 + 2.0*0.5 + 3.0*0.5*0.5 + 4.0*0.5*0.5*0.5, tol );
    TEST_FLOATING_EQUALITY( Sp_view[0],  2.0 + 3.0*2*0.5 + 4.0*3*0.5*0.5, tol );
    TEST_FLOATING_EQUALITY( Spp_view[0], 3.0*2 + 4.0*6*0.5, tol );
  }
  evaluateCubicSpline<double>(coeff, index, t, outArg(*S), outArg(*Sp));
  {
    Thyra::ConstDetachedVectorView<double> S_view(*S);
    Thyra::ConstDetachedVectorView<double> Sp_view(*Sp);
    TEST_FLOATING_EQUALITY( S_view[0],   1.0 + 2.0*0.5 + 3.0*0.5*0.5 + 4.0*0.5*0.5*0.5, tol );
    TEST_FLOATING_EQUALITY( Sp_view[0],  2.0 + 3.0*2*0.5 + 4.0*3*0.5*0.5, tol );
  }
  evaluateCubicSpline<double>(coeff, index, t, outArg(*S));
  {
    Thyra::ConstDetachedVectorView<double> S_view(*S);
    TEST_FLOATING_EQUALITY( S_view[0],   1.0 + 2.0*0.5 + 3.0*0.5*0.5 + 4.0*0.5*0.5*0.5, tol );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_CubicSplineInterpolator, bad_evaluateCubicSpline_1 ) {
    CubicSplineCoeff<double> coeff;
    coeff.t.push_back(0.0);
    coeff.a.push_back(createDefaultVector(1,1.0));
    coeff.b.push_back(createDefaultVector(1,2.0));
    coeff.c.push_back(createDefaultVector(1,3.0));
    coeff.d.push_back(createDefaultVector(1,4.0));
    int index = 0;
    double t = 1.5;
    RCP<Thyra::VectorBase<double> > S   = createDefaultVector(1,0.0);
    TEST_THROW(evaluateCubicSpline<double>(coeff, index, t, outArg(*S)), std::logic_error);
  }
TEUCHOS_UNIT_TEST( Rythmos_CubicSplineInterpolator, bad_evaluateCubicSpline_2 ) {
    CubicSplineCoeff<double> coeff;
    coeff.t.push_back(0.0);
    coeff.t.push_back(1.0);
    coeff.a.push_back(createDefaultVector(1,1.0));
    coeff.a.push_back(createDefaultVector(1,1.0));
    coeff.b.push_back(createDefaultVector(1,2.0));
    coeff.c.push_back(createDefaultVector(1,3.0));
    coeff.d.push_back(createDefaultVector(1,4.0));
    int index = 0;
    double t = 1.5;
    RCP<Thyra::VectorBase<double> > S   = createDefaultVector(1,0.0);
    TEST_THROW(evaluateCubicSpline<double>(coeff, index, t, outArg(*S)), std::logic_error);
  }
TEUCHOS_UNIT_TEST( Rythmos_CubicSplineInterpolator, bad_evaluateCubicSpline_3 ) {
    CubicSplineCoeff<double> coeff;
    coeff.t.push_back(0.0);
    coeff.t.push_back(1.0);
    coeff.a.push_back(createDefaultVector(1,1.0));
    coeff.b.push_back(createDefaultVector(1,2.0));
    coeff.b.push_back(createDefaultVector(1,2.0));
    coeff.c.push_back(createDefaultVector(1,3.0));
    coeff.d.push_back(createDefaultVector(1,4.0));
    int index = 0;
    double t = 1.5;
    RCP<Thyra::VectorBase<double> > S   = createDefaultVector(1,0.0);
    TEST_THROW(evaluateCubicSpline<double>(coeff, index, t, outArg(*S)), std::logic_error);
  }
TEUCHOS_UNIT_TEST( Rythmos_CubicSplineInterpolator, bad_evaluateCubicSpline_4 ) {
    CubicSplineCoeff<double> coeff;
    coeff.t.push_back(0.0);
    coeff.t.push_back(1.0);
    coeff.a.push_back(createDefaultVector(1,1.0));
    coeff.b.push_back(createDefaultVector(1,2.0));
    coeff.c.push_back(createDefaultVector(1,3.0));
    coeff.c.push_back(createDefaultVector(1,3.0));
    coeff.d.push_back(createDefaultVector(1,4.0));
    int index = 0;
    double t = 1.5;
    RCP<Thyra::VectorBase<double> > S   = createDefaultVector(1,0.0);
    TEST_THROW(evaluateCubicSpline<double>(coeff, index, t, outArg(*S)), std::logic_error);
  }
TEUCHOS_UNIT_TEST( Rythmos_CubicSplineInterpolator, bad_evaluateCubicSpline_5 ) {
    CubicSplineCoeff<double> coeff;
    coeff.t.push_back(0.0);
    coeff.t.push_back(1.0);
    coeff.a.push_back(createDefaultVector(1,1.0));
    coeff.b.push_back(createDefaultVector(1,2.0));
    coeff.c.push_back(createDefaultVector(1,3.0));
    coeff.d.push_back(createDefaultVector(1,4.0));
    coeff.d.push_back(createDefaultVector(1,4.0));
    int index = 0;
    double t = 1.5;
    RCP<Thyra::VectorBase<double> > S   = createDefaultVector(1,0.0);
    TEST_THROW(evaluateCubicSpline<double>(coeff, index, t, outArg(*S)), std::logic_error);
  }
TEUCHOS_UNIT_TEST( Rythmos_CubicSplineInterpolator, bad_evaluateCubicSpline_6 ) {
    CubicSplineCoeff<double> coeff;
    coeff.t.push_back(0.0);
    coeff.t.push_back(1.0);
    coeff.a.push_back(createDefaultVector(1,1.0));
    coeff.b.push_back(createDefaultVector(1,2.0));
    coeff.c.push_back(createDefaultVector(1,3.0));
    coeff.d.push_back(createDefaultVector(1,4.0));
    int index = 1;
    double t = 1.5;
    RCP<Thyra::VectorBase<double> > S   = createDefaultVector(1,0.0);
    TEST_THROW(evaluateCubicSpline<double>(coeff, index, t, outArg(*S)), std::out_of_range);
  }
//  {
//    ECHO("TEST 7")
//    CubicSplineCoeff<double> coeff;
//    coeff.t.push_back(0.0);
//    coeff.t.push_back(1.0);
//    coeff.a.push_back(createDefaultVector(1,1.0));
//    coeff.b.push_back(createDefaultVector(1,2.0));
//    coeff.c.push_back(createDefaultVector(1,3.0));
//    coeff.d.push_back(createDefaultVector(1,4.0));
//    int index = -1;
//    double t = 1.5;
//    RCP<Thyra::VectorBase<double> > S   = createDefaultVector(1,0.0);
//#ifdef RYTHMOS_DEBUG
//    TEST_THROW(evaluateCubicSpline<double>(coeff, index, t, outArg(*S)), Teuchos::RangeError);
//#else
//    TEST_THROW(evaluateCubicSpline<double>(coeff, index, t, outArg(*S)), std::out_of_range);
//#endif
//
//  }

TEUCHOS_UNIT_TEST( Rythmos_CubicSplineInterpolator, computeCubicSplineCoeff ) {
  RCP<DataStore<double>::DataStoreVector_t> data_in = rcp( new DataStore<double>::DataStoreVector_t );
  double t0 = 0.0;
  RCP<Thyra::VectorBase<double> > x0 = createDefaultVector<double>(1,0.0);
  RCP<Thyra::VectorBase<double> > xdot0;
  double accuracy0 = 0.0;
  data_in->push_back(DataStore<double>(t0,x0,xdot0,accuracy0));
  double t1 = 1.0;
  RCP<Thyra::VectorBase<double> > x1 = createDefaultVector<double>(1,1.0);
  RCP<Thyra::VectorBase<double> > xdot1;
  double accuracy1 = 0.0;
  data_in->push_back(DataStore<double>(t1,x1,xdot1,accuracy1));
  double t2 = 2.0;
  RCP<Thyra::VectorBase<double> > x2 = createDefaultVector<double>(1,0.0);
  RCP<Thyra::VectorBase<double> > xdot2;
  double accuracy2 = 0.0;
  data_in->push_back(DataStore<double>(t2,x2,xdot2,accuracy2));
  RCP<CubicSplineCoeff<double> > coeff = rcp(new CubicSplineCoeff<double>);
  computeCubicSplineCoeff(*data_in,outArg(*coeff));
  TEST_EQUALITY_CONST( coeff->t.size(), 3 );
  TEST_EQUALITY_CONST( coeff->a.size(), 2 );
  TEST_EQUALITY_CONST( coeff->b.size(), 2 );
  TEST_EQUALITY_CONST( coeff->c.size(), 2 );
  TEST_EQUALITY_CONST( coeff->d.size(), 2 );
  // 11/17/08 tscoffe:  Computed by hand on whiteboard.  The computed
  // polynomials satisfy all the requirements for a cubic spline.
  // i.e. they match the function values, the first and second derivatives of
  // the splines match at the nodes, and the second derivative of boundary
  // splines are zero on the boundaries.
  double tol = 1.0e-12;
  TEST_EQUALITY_CONST( coeff->t[0],  0.0 );
  TEST_EQUALITY_CONST( coeff->t[1],  1.0 );
  TEST_EQUALITY_CONST( coeff->t[2],  2.0 );
  {
    Thyra::ConstDetachedVectorView<double> a0_view(*(coeff->a[0]));
    Thyra::ConstDetachedVectorView<double> b0_view(*(coeff->b[0]));
    Thyra::ConstDetachedVectorView<double> c0_view(*(coeff->c[0]));
    Thyra::ConstDetachedVectorView<double> d0_view(*(coeff->d[0]));
    Thyra::ConstDetachedVectorView<double> a1_view(*(coeff->a[1]));
    Thyra::ConstDetachedVectorView<double> b1_view(*(coeff->b[1]));
    Thyra::ConstDetachedVectorView<double> c1_view(*(coeff->c[1]));
    Thyra::ConstDetachedVectorView<double> d1_view(*(coeff->d[1]));
    TEST_FLOATING_EQUALITY( a0_view[0]+1,  1.0, tol );
    TEST_FLOATING_EQUALITY( b0_view[0],    1.5, tol );
    TEST_FLOATING_EQUALITY( c0_view[0]+1,  1.0, tol );
    TEST_FLOATING_EQUALITY( d0_view[0],   -0.5, tol );
    TEST_FLOATING_EQUALITY( a1_view[0],    1.0, tol );
    TEST_FLOATING_EQUALITY( b1_view[0]+1,  1.0, tol );
    TEST_FLOATING_EQUALITY( c1_view[0],   -1.5, tol );
    TEST_FLOATING_EQUALITY( d1_view[0],    0.5, tol );
  }
  // Verify cubic spline properties:
  RCP<Thyra::VectorBase<double> > S = createDefaultVector(1,0.0);
  RCP<Thyra::VectorBase<double> > Sp = createDefaultVector(1,0.0);
  RCP<Thyra::VectorBase<double> > Spp = createDefaultVector(1,0.0);
  int j = 0;
  double time = 0.0;
  evaluateCubicSpline<double>(*coeff, j, time, outArg(*S), outArg(*Sp), outArg(*Spp));
  {
    Thyra::ConstDetachedVectorView<double> S_view(*S);
    Thyra::ConstDetachedVectorView<double> Sp_view(*Sp);
    Thyra::ConstDetachedVectorView<double> Spp_view(*Spp);
    TEST_FLOATING_EQUALITY( S_view[0]+1,   1.0, tol ); // function value at t0
    TEST_FLOATING_EQUALITY( Sp_view[0],    1.5, tol ); // irrelevant number
    TEST_FLOATING_EQUALITY( Spp_view[0]+1, 1.0, tol ); // second derivative is zero at boundary
  }
  j = 0;
  time = 1.0;
  evaluateCubicSpline<double>(*coeff, j, time, outArg(*S), outArg(*Sp), outArg(*Spp));
  {
    Thyra::ConstDetachedVectorView<double> S_view(*S);
    Thyra::ConstDetachedVectorView<double> Sp_view(*Sp);
    Thyra::ConstDetachedVectorView<double> Spp_view(*Spp);
    TEST_FLOATING_EQUALITY( S_view[0],    1.0, tol ); // function value at t1
    TEST_FLOATING_EQUALITY( Sp_view[0]+1, 1.0, tol ); // first derivatives match at nodes
    TEST_FLOATING_EQUALITY( Spp_view[0], -3.0, tol ); // second derivatives match at nodes
  }
  j = 1;
  time = 1.0;
  evaluateCubicSpline<double>(*coeff, j, time, outArg(*S), outArg(*Sp), outArg(*Spp));
  {
    Thyra::ConstDetachedVectorView<double> S_view(*S);
    Thyra::ConstDetachedVectorView<double> Sp_view(*Sp);
    Thyra::ConstDetachedVectorView<double> Spp_view(*Spp);
    TEST_FLOATING_EQUALITY( S_view[0],    1.0, tol ); // function value at t1
    TEST_FLOATING_EQUALITY( Sp_view[0]+1, 1.0, tol ); // first derivatives match at nodes
    TEST_FLOATING_EQUALITY( Spp_view[0], -3.0, tol ); // second derivatives match at nodes
  }
  j = 1;
  time = 2.0;
  evaluateCubicSpline<double>(*coeff, j, time, outArg(*S), outArg(*Sp), outArg(*Spp));
  {
    Thyra::ConstDetachedVectorView<double> S_view(*S);
    Thyra::ConstDetachedVectorView<double> Sp_view(*Sp);
    Thyra::ConstDetachedVectorView<double> Spp_view(*Spp);
    TEST_FLOATING_EQUALITY( S_view[0]+1,   1.0, tol ); // function value at t2
    TEST_FLOATING_EQUALITY( Sp_view[0],   -1.5, tol ); // irrelevant number
    TEST_FLOATING_EQUALITY( Spp_view[0]+1, 1.0, tol ); // second derivative is zero at boundary
  }
}

TEUCHOS_UNIT_TEST( Rythmos_CubicSplineInterpolator, bad_computeCubicSplineCoeff) {
  RCP<DataStore<double>::DataStoreVector_t> data_in = rcp( new DataStore<double>::DataStoreVector_t );
  double t0 = 0.0;
  RCP<Thyra::VectorBase<double> > x0 = createDefaultVector<double>(1,0.0);
  RCP<Thyra::VectorBase<double> > xdot0;
  double accuracy0 = 0.0;
  data_in->push_back(DataStore<double>(t0,x0,xdot0,accuracy0));
  RCP<CubicSplineCoeff<double> > coeff = rcp(new CubicSplineCoeff<double>);
  // not enough points
  TEST_THROW( computeCubicSplineCoeff(*data_in,outArg(*coeff)), std::logic_error);
  double t1 = 1.0;
  RCP<Thyra::VectorBase<double> > x1 = createDefaultVector<double>(1,1.0);
  RCP<Thyra::VectorBase<double> > xdot1;
  double accuracy1 = 0.0;
  data_in->push_back(DataStore<double>(t1,x1,xdot1,accuracy1));
  coeff = rcp(new CubicSplineCoeff<double>);
  // enough points for linear special case
  TEST_NOTHROW( computeCubicSplineCoeff(*data_in,outArg(*coeff)));
  double t2 = 2.0;
  RCP<Thyra::VectorBase<double> > x2 = createDefaultVector<double>(1,0.0);
  RCP<Thyra::VectorBase<double> > xdot2;
  double accuracy2 = 0.0;
  data_in->push_back(DataStore<double>(t2,x2,xdot2,accuracy2));
  coeff = rcp(new CubicSplineCoeff<double>);
  // enough points for cubic spline
  TEST_NOTHROW( computeCubicSplineCoeff(*data_in,outArg(*coeff)) );

  double t3 = 2.0;
  RCP<Thyra::VectorBase<double> > x3 = createDefaultVector<double>(1,0.0);
  RCP<Thyra::VectorBase<double> > xdot3;
  double accuracy3 = 0.0;
  data_in->push_back(DataStore<double>(t3,x3,xdot3,accuracy3));
  coeff = rcp(new CubicSplineCoeff<double>);
  // non unique time values
#ifdef RYTHMOS_DEBUG
  TEST_THROW( computeCubicSplineCoeff(*data_in,outArg(*coeff)), std::logic_error);
#else 
  TEST_NOTHROW( computeCubicSplineCoeff(*data_in,outArg(*coeff)) );
  TEST_EQUALITY_CONST( coeff->t.size(), 4 );
  for (int i=0; i<3 ; ++i) {
    Thyra::ConstDetachedVectorView<double> b_view(*(coeff->b[i]));
    Thyra::ConstDetachedVectorView<double> c_view(*(coeff->c[i]));
    Thyra::ConstDetachedVectorView<double> d_view(*(coeff->d[i]));
    typedef Teuchos::ScalarTraits<double> ST;
    TEST_EQUALITY_CONST( ST::isnaninf(b_view[0]), true );
    TEST_EQUALITY_CONST( ST::isnaninf(c_view[0]), true );
    TEST_EQUALITY_CONST( ST::isnaninf(d_view[0]), true );
  }
#endif // RYTHMOS_DEBUG

  data_in->pop_back();

  coeff = rcp(new CubicSplineCoeff<double>);
  // back to normal
  TEST_NOTHROW( computeCubicSplineCoeff(*data_in,outArg(*coeff)) );

  double t4 = 1.5;
  RCP<Thyra::VectorBase<double> > x4 = createDefaultVector<double>(1,0.0);
  RCP<Thyra::VectorBase<double> > xdot4;
  double accuracy4 = 0.0;
  data_in->push_back(DataStore<double>(t4,x4,xdot4,accuracy4));
  coeff = rcp(new CubicSplineCoeff<double>);
  // non sorted time values
#ifdef RYTHMOS_DEBUG
  TEST_THROW( computeCubicSplineCoeff(*data_in,outArg(*coeff)), std::logic_error);
#else
  TEST_NOTHROW( computeCubicSplineCoeff(*data_in,outArg(*coeff)) );
  TEST_EQUALITY_CONST( coeff->t.size(), 4 );
  {
    Thyra::ConstDetachedVectorView<double> b_view(*(coeff->b[2]));
    Thyra::ConstDetachedVectorView<double> c_view(*(coeff->c[2]));
    Thyra::ConstDetachedVectorView<double> d_view(*(coeff->d[2]));
    TEST_EQUALITY_CONST( b_view[0],  2.0 );
    TEST_EQUALITY_CONST( c_view[0],  6.0 );
    TEST_EQUALITY_CONST( d_view[0],  4.0 );
  }
#endif // RYTHMOS_DEBUG

}

TEUCHOS_UNIT_TEST( Rythmos_CubicSplineInterpolator, interpolate ) {
  RCP<InterpolatorBase<double> > csi = cubicSplineInterpolator<double>();
  double maxOrder = csi->order();
  for (int order = 0 ; order <= maxOrder+1 ; ++order) {
    TEST_EQUALITY( order, order );
    Polynomial<double> p(order,1.0); 
    RCP<DataStore<double>::DataStoreVector_t> data_in = rcp( new DataStore<double>::DataStoreVector_t );
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
    interpolate<double>(*csi, data_in, t_values, &data_out);
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

TEUCHOS_UNIT_TEST( Rythmos_CubicSplineInterpolator, bad_interpolate ) {
  RCP<InterpolatorBase<double> > csi = cubicSplineInterpolator<double>();
  Array<double> t_values;
  DataStore<double>::DataStoreVector_t data_out;
  // must call setNodes before interpolate
  TEST_THROW(csi->interpolate(t_values,&data_out),std::logic_error);

  RCP<DataStore<double>::DataStoreVector_t> data_in = rcp( new DataStore<double>::DataStoreVector_t );

  double t0 = 0.0;
  RCP<Thyra::VectorBase<double> > x0 = createDefaultVector<double>(1,0.0);
  RCP<Thyra::VectorBase<double> > xdot0;
  double accuracy0 = 0.0;
  data_in->push_back(DataStore<double>(t0,x0,xdot0,accuracy0));
  csi->setNodes(data_in);

  t_values.push_back(0.5);
  // Not enough nodes to evaluate t=0.5
#ifdef RYTHMOS_DEBUG
  TEST_THROW(csi->interpolate(t_values,&data_out),std::logic_error);
#else // RYTHMOS_DEBUG
  // In this case, we're not checking the preconditions, so we assume the point you asked for is the only point we have
  TEST_NOTHROW(csi->interpolate(t_values,&data_out));
  {
    Thyra::ConstDetachedVectorView<double> x_view(*data_out[0].x);
    TEST_EQUALITY_CONST( x_view[0], 0.0 );
  }
#endif // RYTHMOS_DEBUG

  double t1 = 1.0;
  RCP<Thyra::VectorBase<double> > x1 = createDefaultVector<double>(1,1.0);
  RCP<Thyra::VectorBase<double> > xdot1;
  double accuracy1 = 0.0;
  data_in->push_back(DataStore<double>(t1,x1,xdot1,accuracy1));
  csi->setNodes(data_in);

  // Linear interpolation used with only two points
  TEST_NOTHROW(csi->interpolate(t_values,&data_out));
  double tol = 1.0e-15;
  {
    Thyra::ConstDetachedVectorView<double> x_view(*data_out[0].x);
    TEST_FLOATING_EQUALITY( x_view[0], 0.5, tol );
  }
  t_values[0] = 0.0;
  TEST_NOTHROW(csi->interpolate(t_values,&data_out));
  {
    Thyra::ConstDetachedVectorView<double> x_view(*data_out[0].x);
    TEST_FLOATING_EQUALITY( x_view[0]+1, 1.0, tol );
  }
  t_values[0] = 1.0;
  TEST_NOTHROW(csi->interpolate(t_values,&data_out));
  {
    Thyra::ConstDetachedVectorView<double> x_view(*data_out[0].x);
    TEST_FLOATING_EQUALITY( x_view[0], 1.0, tol );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_CubicSplineInterpolator, changedData ) {
  // Create three points of data and call setNodes & interpolate on csi
  // Delete the 3rd point and add a different value and call interpolate again
  // Delete the 3rd point and add a different time point and call interpolate again
  // Add a fourth point and call interpolate again
  // Delete the 1st point and call interpolate again
  RCP<InterpolatorBase<double> > csi = cubicSplineInterpolator<double>();
  RCP<DataStore<double>::DataStoreVector_t> data_in = rcp( new DataStore<double>::DataStoreVector_t );

  double t0 = 0.0;
  RCP<Thyra::VectorBase<double> > x0 = createDefaultVector<double>(1,0.0);
  RCP<Thyra::VectorBase<double> > xdot0;
  double accuracy0 = 0.0;
  data_in->push_back(DataStore<double>(t0,x0,xdot0,accuracy0));

  double t1 = 1.0;
  RCP<Thyra::VectorBase<double> > x1 = createDefaultVector<double>(1,1.0);
  RCP<Thyra::VectorBase<double> > xdot1;
  double accuracy1 = 0.0;
  data_in->push_back(DataStore<double>(t1,x1,xdot1,accuracy1));

  double t2 = 2.0;
  RCP<Thyra::VectorBase<double> > x2 = createDefaultVector<double>(1,0.0);
  RCP<Thyra::VectorBase<double> > xdot2;
  double accuracy2 = 0.0;
  data_in->push_back(DataStore<double>(t2,x2,xdot2,accuracy2));
  csi->setNodes(data_in);
  //ECHO( "data_in = [(0,0), (1,1), (2,0)]" );

  Array<double> t_values;
  DataStore<double>::DataStoreVector_t data_out;
  t_values.push_back(0.5);
  //ECHO( "Calling interpolate at t=0.5" );
  // This is a valid interpolate call
  TEST_NOTHROW(csi->interpolate(t_values,&data_out));
  double tol = 1.0e-15;
  double t05value = 0.6875;
  {
    Thyra::ConstDetachedVectorView<double> x_view(*data_out[0].x);
    TEST_FLOATING_EQUALITY( x_view[0], t05value, tol );
  }

  // Now we have computed the cubic spline coefficients.
  // Delete the last node:
  data_in->pop_back();
  //ECHO( "data_in = [(0,0), (1,1)]" );
  //ECHO( "Calling interpolate at t=0.5 after deleting last node, should throw" );
#ifdef RYTHMOS_DEBUG
  TEST_THROW(csi->interpolate(t_values,&data_out),std::logic_error);
#else // RYTHMOS_DEBUG
  // spline coefficients did not change, so you get the same output
  TEST_NOTHROW(csi->interpolate(t_values,&data_out));
  {
    Thyra::ConstDetachedVectorView<double> x_view(*data_out[0].x);
    TEST_FLOATING_EQUALITY( x_view[0], t05value, tol );
  }
#endif // RYTHMOS_DEBUG

  // Change the value of the third point
  V_S(outArg(*x2),1.0);
  data_in->push_back(DataStore<double>(t2,x2,xdot2,accuracy2));
  //ECHO( "data_in = [(0,0), (1,1), (2,1)]" );
  //ECHO( "Calling interpolate at t=1.5 after adding different last node, should throw" );
#ifdef RYTHMOS_DEBUG
  TEST_THROW(csi->interpolate(t_values,&data_out),std::logic_error);
#else // RYTHMOS_DEBUG
  // spline coefficients did not change, so you get the same output
  TEST_NOTHROW(csi->interpolate(t_values,&data_out));
  {
    Thyra::ConstDetachedVectorView<double> x_view(*data_out[0].x);
    TEST_FLOATING_EQUALITY( x_view[0], t05value, tol );
  }
#endif // RYTHMOS_DEBUG


  // Change the time-point of the third point
  data_in->pop_back();
  //ECHO( "data_in = [(0,0), (1,1)]" );
  t2 = 2.5;
  data_in->push_back(DataStore<double>(t2,x2,xdot2,accuracy2));
  //ECHO( "data_in = [(0,0), (1,1), (2.5,1)]" );
  //ECHO( "Calling interpolate at t=1.5 after changing last node to t=2.5, should throw" );
#ifdef RYTHMOS_DEBUG
  TEST_THROW(csi->interpolate(t_values,&data_out),std::logic_error);
#else // RYTHMOS_DEBUG
  // spline coefficients did not change, so you get the same output
  TEST_NOTHROW(csi->interpolate(t_values,&data_out));
  {
    Thyra::ConstDetachedVectorView<double> x_view(*data_out[0].x);
    TEST_FLOATING_EQUALITY( x_view[0], t05value, tol );
  }
#endif // RYTHMOS_DEBUG

  // Add a fourth point:
  double t3 = 3.0;
  RCP<Thyra::VectorBase<double> > x3 = createDefaultVector<double>(1,3.0);
  RCP<Thyra::VectorBase<double> > xdot3;
  double accuracy3 = 0.0;
  data_in->push_back(DataStore<double>(t3,x3,xdot3,accuracy3));
  //ECHO( "data_in = [(0,0), (1,1), (2.5,1), (3,3)]" );

  t_values[0] = 2.7;
  //ECHO( "Calling interpolate at t=2.7 after adding a fourth point, should throw" );
#ifdef RYTHMOS_DEBUG
  TEST_THROW(csi->interpolate(t_values,&data_out), std::logic_error);
#else // RYTHMOS_DEBUG
  // In this case, t_values[0]=2.7 is out of range for the original [0,2] range
  // that the spline coefficients were computed for.  So we should get an
  // out_of_range exception from the cubic spline interpolator.
  TEST_THROW(csi->interpolate(t_values,&data_out), std::out_of_range);
#endif // RYTHMOS_DEBUG

  // Delete the 1st point:
  data_in->erase(data_in->begin());
  //ECHO( "data_in = [(1,1), (2.5,1), (3,3)]" );
  t_values[0] = 0.5;
  //ECHO( "Calling interpolate at t=0.5 after deleting the first point, should throw" );
#ifdef RYTHMOS_DEBUG
  TEST_THROW(csi->interpolate(t_values,&data_out), std::logic_error);
#else // RYTHMOS_DEBUG
  // In this case, the InterpolationBuffer searches for an interval containing
  // this point and can't find any, so it returns no points.
  TEST_NOTHROW(csi->interpolate(t_values,&data_out));
  TEST_EQUALITY_CONST( data_out.size(), 0 );
#endif // RYTHMOS_DEBUG

}

} // namespace Rythmos



