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

#include "Rythmos_UnitTestHelpers.hpp"

#include "Rythmos_InterpolationBuffer.hpp"
#include "Rythmos_HermiteInterpolator.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Rythmos_LinearInterpolator.hpp"

namespace Rythmos {

// member functions to test:
// get_x_space
// constructor TESTED
// initialize TESTED through constructor
// setInterpolator TESTED
// unSetInterpolator TESTED
// setStorage TESTED
// getStorage TESTED
// addPoints TESTED
//   interaction with IBPolicy TESTED
// getPoints  TESTED
// getTimeRange TESTED
// getNodes
// getOrder TESTED
// removeNodes
// description TESTED
// describe
// setParameterList (StorageLimit [TESTED], InterpolationBufferPolicy [TESTED])
// getNonconstParameterList
// unsetParameterList
// getValidParameters TESTED
// IBPolicy (set through ParameterList) TESTED
// getIBPolicy TESTED
// high verbosity output


using Teuchos::as;


TEUCHOS_UNIT_TEST( Rythmos_InterpolationBuffer, newBuffer ) {
  InterpolationBuffer<double> ib;
  TEST_EQUALITY_CONST( ib.getStorage(), 2 );
  TEST_EQUALITY_CONST( ib.getTimeRange().isValid(), false ); 
  TEST_EQUALITY_CONST( ib.getOrder(), 1 ); // linear interpolator by default
  TEST_EQUALITY( ib.getParameterList(), Teuchos::null );
}

TEUCHOS_UNIT_TEST( Rythmos_InterpolationBuffer, nonMemberConstructor ) {
  // Basic initial test
  RCP<InterpolationBuffer<double> > ib = interpolationBuffer<double>();
  TEST_EQUALITY_CONST( ib->getOrder(), 1 );
  TEST_EQUALITY_CONST( ib->getStorage(), 2 );
  TEST_EQUALITY_CONST( ib->getTimeRange().isValid(), false );
  TEST_EQUALITY_CONST( ib->getParameterList(), Teuchos::null );

  // modified initial test
  ib = interpolationBuffer<double>(Teuchos::null,6);
  TEST_EQUALITY_CONST( ib->getStorage(), 6 );

  // Set up specfic interpolator with different order
  RCP<HermiteInterpolator<double> > hi = hermiteInterpolator<double>();
  ib = interpolationBuffer<double>(hi,13);
  TEST_EQUALITY_CONST( ib->getOrder(), 3 );
}

TEUCHOS_UNIT_TEST( Rythmos_InterpolationBuffer, setInterpolator ) {
  RCP<InterpolationBuffer<double> > ib = interpolationBuffer<double>(Teuchos::null,13);
  RCP<HermiteInterpolator<double> > hi = hermiteInterpolator<double>();
  TEST_EQUALITY_CONST( ib->getOrder(), 1 );
  ib->setInterpolator(hi);
  TEST_EQUALITY_CONST( ib->getOrder(), 3 );

  // Test unSetInterpolator
  RCP<InterpolatorBase<double> > interp = ib->unSetInterpolator();
  RCP<HermiteInterpolator<double> > hi2 = Teuchos::rcp_dynamic_cast<HermiteInterpolator<double> >(interp);
  TEST_EQUALITY_CONST( Teuchos::is_null(hi2), false );
  TEST_EQUALITY( hi2, hi );
  TEST_EQUALITY_CONST( ib->getOrder(), 1 );

  // Test setInterpolator(Teuchos::null)
  ib = interpolationBuffer<double>(hi,8);
  TEST_EQUALITY_CONST( ib->getOrder(), 3 );
  ib->setInterpolator(Teuchos::null);
  TEST_EQUALITY_CONST( ib->getOrder(), 1 );
}

TEUCHOS_UNIT_TEST( Rythmos_InterpolationBuffer, setStorage) {
  InterpolationBuffer<double> ib;
  ib.setStorage( 5 );
  TEST_EQUALITY_CONST( ib.getStorage(), 5 );
  ib.setStorage( 1 );
  TEST_EQUALITY_CONST( ib.getStorage(), 2 );

  RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
  pl->set("StorageLimit",13);
  ib.setParameterList(pl);
  TEST_EQUALITY_CONST(ib.getStorage(), 13 );
  pl->set("StorageLimit",1);
  ib.setParameterList(pl);
  TEST_EQUALITY_CONST(ib.getStorage(), 2 );
  pl->set("StorageLimit",0);
  ib.setParameterList(pl);
  TEST_EQUALITY_CONST(ib.getStorage(), 2 );
}


// Verify we can append entries to an IB correctly
//   * try passing a range of points that includes the last point in the IB
TEUCHOS_UNIT_TEST( Rythmos_InterpolationBuffer, invalidAddVectors ) {
  Array<double> time_vec;
  time_vec.push_back(2.0);
  time_vec.push_back(3.0);

  Array<RCP<const VectorBase<double> > > v_vec;
  RCP<VectorBase<double> > v = createDefaultVector(2,2.0);
  v_vec.push_back(v);
  v = createDefaultVector(2,3.0);
  v_vec.push_back(v);

  Array<RCP<const VectorBase<double> > > vdot_vec;
  RCP<VectorBase<double> > vdot = createDefaultVector(2,4.0);
  vdot_vec.push_back(vdot);
  vdot = createDefaultVector(2,5.0);
  vdot_vec.push_back(vdot);

  InterpolationBuffer<double> ib;
  ib.addPoints(time_vec,v_vec,vdot_vec);

  Array<double> newTime_vec;
  newTime_vec.push_back(3.0);
  newTime_vec.push_back(4.0);

  Array<RCP<const VectorBase<double> > > newV_vec;
  RCP<VectorBase<double> > newV = createDefaultVector(2,6.0);
  newV_vec.push_back(newV);
  newV = createDefaultVector(2,7.0);
  newV_vec.push_back(newV);

  Array<RCP<const VectorBase<double> > > newVdot_vec;
  RCP<VectorBase<double> > newVdot = createDefaultVector(2,8.0);
  newVdot_vec.push_back(newVdot);
  newVdot = createDefaultVector(2,9.0);
  newVdot_vec.push_back(newVdot);

  TEST_THROW( ib.addPoints(newTime_vec, newV_vec, newVdot_vec), std::logic_error);
}

// Verify that the IB is copying the vectors rather than just storing the pointers
TEUCHOS_UNIT_TEST( Rythmos_InterpolationBuffer, copyVectors ) {
  double t = 2.5;
  Array<double> time_vec;
  time_vec.push_back(t);

  RCP<VectorBase<double> > v = createDefaultVector(2,2.0);
  Array<RCP<const VectorBase<double> > > v_vec;
  v_vec.push_back(v);

  RCP<VectorBase<double> > v_dot = createDefaultVector(2,3.0);
  Array<RCP<const VectorBase<double> > > v_dot_vec;
  v_dot_vec.push_back(v_dot);

  InterpolationBuffer<double> ib;
  ib.addPoints(time_vec, v_vec, v_dot_vec);

  Thyra::V_S(v.ptr(), 4.0);
  Thyra::V_S(v_dot.ptr(), 5.0);

  Array<RCP<const VectorBase<double> > > v_vec_out;
  Array<RCP<const VectorBase<double> > > v_dot_vec_out;
  Array<double> accuracy_vec_out;
  ib.getPoints(time_vec, &v_vec_out, &v_dot_vec_out, &accuracy_vec_out);

  TEST_EQUALITY_CONST( get_ele(*(v_vec_out[0]),0), 2.0 );
  TEST_EQUALITY_CONST( get_ele(*(v_dot_vec_out[0]),0), 3.0 );
}

TEUCHOS_UNIT_TEST( Rythmos_InterpolationBuffer, setIBPolicy ) {
  RCP<InterpolationBuffer<double> > ib = interpolationBuffer<double>();
  TEST_EQUALITY_CONST( ib->getIBPolicy(), BUFFER_POLICY_KEEP_NEWEST );
  RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
  pl->set("InterpolationBufferPolicy", "Static Policy");
  ib->setParameterList(pl);
  TEST_EQUALITY_CONST( ib->getIBPolicy(), BUFFER_POLICY_STATIC );
  pl->set("InterpolationBufferPolicy", "Invalid Policy");
  ib->setParameterList(pl);
  TEST_EQUALITY_CONST( ib->getIBPolicy(), BUFFER_POLICY_STATIC );
}

TEUCHOS_UNIT_TEST( Rythmos_InterpolationBuffer, setParameterList ) {
  RCP<InterpolationBuffer<double> > ib = interpolationBuffer<double>();
  TEST_EQUALITY_CONST( ib->getIBPolicy(), BUFFER_POLICY_KEEP_NEWEST );
  TEST_EQUALITY_CONST( ib->getStorage(), 2 );
  {
    // Verify setting empty parameter list on default IB does not change values
    RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    ib->setParameterList(pl);
    TEST_EQUALITY_CONST( ib->getIBPolicy(), BUFFER_POLICY_KEEP_NEWEST );
    TEST_EQUALITY_CONST( ib->getStorage(), 2 );
  }
  {
    // Verify setting empty parameter list on non-default IB does not change values
    RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    pl->set("InterpolationBufferPolicy", "Static Policy");
    pl->set("StorageLimit", 15);
    ib->setParameterList(pl);
    TEST_EQUALITY_CONST( ib->getIBPolicy(), BUFFER_POLICY_STATIC );
    TEST_EQUALITY_CONST( ib->getStorage(), 15 );
    pl = Teuchos::parameterList();
    ib->setParameterList(pl);
    TEST_EQUALITY_CONST( ib->getIBPolicy(), BUFFER_POLICY_KEEP_NEWEST);
    TEST_EQUALITY_CONST( ib->getStorage(), 2 );
  }
  {
    // Verify that invalid policy is not accepted
    RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    pl->set("InterpolationBufferPolicy", "Invalid Policy");
    ib->setParameterList(pl);
    TEST_EQUALITY_CONST( ib->getIBPolicy(), BUFFER_POLICY_KEEP_NEWEST );
  }
  {
    // Verify that invalid storage is fixed
    RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    pl->set("StorageLimit", 1);
    ib->setParameterList(pl);
    TEST_EQUALITY_CONST( ib->getStorage(), 2 );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_InterpolationBuffer, getValidParameters ) {
  RCP<InterpolationBuffer<double> > ib = interpolationBuffer<double>();
  RCP<const ParameterList> pl = ib->getValidParameters();
  TEST_EQUALITY_CONST( is_null(pl), false );
  TEST_EQUALITY_CONST( pl->isSublist("VerboseObject"), true );
  TEST_EQUALITY_CONST( pl->isParameter("InterpolationBufferPolicy"), true );
  TEST_EQUALITY_CONST( pl->isParameter("StorageLimit"), true );
}

TEUCHOS_UNIT_TEST( Rythmos_InterpolationBuffer, addPoints_basic ) {
  RCP<InterpolationBuffer<double> > ib = interpolationBuffer<double>();
  RCP<VectorBase<double> > v = createDefaultVector(2,1.0);
  RCP<VectorBase<double> > vdot = createDefaultVector(2,2.0);
  Array<double> time_vec;
  time_vec.push_back(0.0);
  Array<RCP<const VectorBase<double> > > x_vec;
  x_vec.push_back(v);
  Array<RCP<const VectorBase<double> > > xdot_vec;
  xdot_vec.push_back(vdot);
  ib->addPoints(
      time_vec,
      x_vec,
      xdot_vec
      );
  Array<double> nodes;
  ib->getNodes( &nodes );
  TEST_EQUALITY_CONST( nodes.size(), 1 );
  TEST_EQUALITY_CONST( nodes[0], 0.0 );
  Array<RCP<const VectorBase<double> > > x_vec_out;
  Array<RCP<const VectorBase<double> > > xdot_vec_out;
  Array<double> accuracy_out;
  ib->getPoints(time_vec, &x_vec_out, &xdot_vec_out, &accuracy_out);
  TEST_EQUALITY_CONST( x_vec_out.size(), 1 );
  {
    Thyra::ConstDetachedVectorView<double> x_vec_out_view( *x_vec_out[0] );
    TEST_EQUALITY_CONST( x_vec_out_view[0], 1.0 );
    TEST_EQUALITY_CONST( x_vec_out_view[1], 1.0 );
    Thyra::ConstDetachedVectorView<double> xdot_vec_out_view( *xdot_vec_out[0] );
    TEST_EQUALITY_CONST( xdot_vec_out_view[0], 2.0 );
    TEST_EQUALITY_CONST( xdot_vec_out_view[1], 2.0 );
  }
  TimeRange<double> tr = ib->getTimeRange();
  TEST_EQUALITY_CONST( tr.lower(), 0.0 );
  TEST_EQUALITY_CONST( tr.upper(), 0.0 );
}

TEUCHOS_UNIT_TEST( Rythmos_InterpolationBuffer, addPoints_more ) {
  RCP<InterpolationBuffer<double> > ib = interpolationBuffer<double>();
  int N = 5;
  ib->setStorage(N);
  for (int i=0 ; i<N ; ++i) {
    Array<double> time_vec;
    time_vec.push_back(0.5+0.25*i);
    RCP<VectorBase<double> > v = createDefaultVector(2,1.0+i*1.0);
    RCP<VectorBase<double> > vdot = createDefaultVector(2,2.0+i*1.0);
    Array<RCP<const VectorBase<double> > > x_vec;
    x_vec.push_back(v);
    Array<RCP<const VectorBase<double> > > xdot_vec;
    xdot_vec.push_back(vdot);
    ib->addPoints(
        time_vec,
        x_vec,
        xdot_vec
        );
    TimeRange<double> tr = ib->getTimeRange();
    TEST_EQUALITY_CONST( tr.lower(), 0.5 );
    TEST_EQUALITY_CONST( tr.upper(), 0.5+0.25*i );
  }
  Array<double> time_vec;
  time_vec.push_back(0.50);
  time_vec.push_back(0.75);
  time_vec.push_back(1.00);
  time_vec.push_back(1.25);
  time_vec.push_back(1.50);
  Array<RCP<const VectorBase<double> > > x_vec_out;
  Array<RCP<const VectorBase<double> > > xdot_vec_out;
  Array<double> accuracy_out;
  ib->getPoints(time_vec, &x_vec_out, &xdot_vec_out, &accuracy_out);
  TEST_EQUALITY_CONST( as<int>(x_vec_out.size()), N );
  for (int i=0 ; i<N ; ++i) {
    Thyra::ConstDetachedVectorView<double> x_vec_out_view( *x_vec_out[i] );
    TEST_EQUALITY_CONST( x_vec_out_view[0], 1.0+i*1.0 );
    TEST_EQUALITY_CONST( x_vec_out_view[1], 1.0+i*1.0 );
    Thyra::ConstDetachedVectorView<double> xdot_vec_out_view( *xdot_vec_out[i] );
    TEST_EQUALITY_CONST( xdot_vec_out_view[0], 2.0+i*1.0 );
    TEST_EQUALITY_CONST( xdot_vec_out_view[1], 2.0+i*1.0 );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_InterpolationBuffer, addPoints_different ) {
  RCP<InterpolationBuffer<double> > ib = interpolationBuffer<double>();
  int N = 5;
  ib->setStorage(N);
  Array<double> time_vec;
  Array<RCP<const VectorBase<double> > > x_vec;
  Array<RCP<const VectorBase<double> > > xdot_vec;
  for (int i=0 ; i<N ; ++i) {
    time_vec.push_back(0.5+0.25*i);
    RCP<VectorBase<double> > v = createDefaultVector(2,1.0+i*1.0);
    RCP<VectorBase<double> > vdot = createDefaultVector(2,2.0+i*1.0);
    x_vec.push_back(v);
    xdot_vec.push_back(vdot);
  }
  ib->addPoints(
      time_vec,
      x_vec,
      xdot_vec
      );
  TimeRange<double> tr = ib->getTimeRange();
  TEST_EQUALITY_CONST( tr.lower(), 0.5 );
  TEST_EQUALITY_CONST( tr.upper(), 1.5 );

  Array<RCP<const VectorBase<double> > > x_vec_out;
  Array<RCP<const VectorBase<double> > > xdot_vec_out;
  Array<double> accuracy_out;
  ib->getPoints(time_vec, &x_vec_out, &xdot_vec_out, &accuracy_out);
  TEST_EQUALITY_CONST( as<int>(x_vec_out.size()), N );
  for (int i=0 ; i<N ; ++i) {
    Thyra::ConstDetachedVectorView<double> x_vec_out_view( *x_vec_out[i] );
    TEST_EQUALITY_CONST( x_vec_out_view[0], 1.0+i*1.0 );
    TEST_EQUALITY_CONST( x_vec_out_view[1], 1.0+i*1.0 );
    Thyra::ConstDetachedVectorView<double> xdot_vec_out_view( *xdot_vec_out[i] );
    TEST_EQUALITY_CONST( xdot_vec_out_view[0], 2.0+i*1.0 );
    TEST_EQUALITY_CONST( xdot_vec_out_view[1], 2.0+i*1.0 );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_InterpolationBuffer, addPoints_bad ) {
  RCP<InterpolationBuffer<double> > ib = interpolationBuffer<double>();
  int N = 5;
  ib->setStorage(N+1);
  Array<double> time_vec;
  Array<RCP<const VectorBase<double> > > x_vec;
  Array<RCP<const VectorBase<double> > > xdot_vec;
  for (int i=0 ; i<N ; ++i) {
    RCP<VectorBase<double> > v = createDefaultVector(2,1.0+i*1.0);
    RCP<VectorBase<double> > vdot = createDefaultVector(2,2.0+i*1.0);
    time_vec.push_back(0.5+0.25*i);
    x_vec.push_back(v);
    xdot_vec.push_back(vdot);
  }
  ib->addPoints(
      time_vec,
      x_vec,
      xdot_vec
      );
  time_vec.clear();
  x_vec.clear();
  xdot_vec.clear();
  time_vec.push_back(0.51);
  x_vec.push_back(createDefaultVector(2,15.0));
  xdot_vec.push_back(createDefaultVector(2,30.0));
#ifdef RYTHMOS_DEBUG
  // inside time range at beginning
  TEST_THROW( ib->addPoints(time_vec,x_vec,xdot_vec), std::out_of_range );
#else // RYTHMOS_DEBUG
  // TODO:  What happens in this case?
#endif // RYTHMOS_DEBUG
  time_vec[0] = 1.49;
  // inside time range at end
#ifdef RYTHMOS_DEBUG
  TEST_THROW( ib->addPoints(time_vec,x_vec,xdot_vec), std::out_of_range );
#else // RYTHMOS_DEBUG
  // TODO:  What happens in this case?
#endif // RYTHMOS_DEBUG
  time_vec.clear();
  // time_vec.size() == 0
#ifdef RYTHMOS_DEBUG
  TEST_THROW( ib->addPoints(time_vec,x_vec,xdot_vec), std::logic_error );
#else // RYTHMOS_DEBUG
  // TODO:  What happens in this case?
#endif // RYTHMOS_DEBUG
  time_vec.push_back(0.0);
  x_vec.clear();
  // x_vec.size() == 0
#ifdef RYTHMOS_DEBUG
  TEST_THROW( ib->addPoints(time_vec,x_vec,xdot_vec), std::logic_error );
#else // RYTHMOS_DEBUG
  // TODO:  What happens in this case?
#endif // RYTHMOS_DEBUG
  x_vec.push_back(createDefaultVector(2,16.0));
  xdot_vec.clear();
  // xdot_vec.size() == 0
#ifdef RYTHMOS_DEBUG
  TEST_THROW( ib->addPoints(time_vec,x_vec,xdot_vec), std::logic_error );
#else // RYTHMOS_DEBUG
  // TODO:  What happens in this case?
#endif // RYTHMOS_DEBUG
  xdot_vec.push_back(createDefaultVector(2,31.0));
  time_vec[0] = 0.5;
  // Replace first value
#ifdef RYTHMOS_DEBUG
  TEST_THROW( ib->addPoints(time_vec,x_vec,xdot_vec), std::out_of_range );
#else // RYTHMOS_DEBUG
  // TODO:  What happens in this case?
#endif // RYTHMOS_DEBUG
  time_vec[0] = 1.5;
  // Replace last value
#ifdef RYTHMOS_DEBUG
  TEST_THROW( ib->addPoints(time_vec,x_vec,xdot_vec), std::out_of_range );
#else // RYTHMOS_DEBUG
  // TODO:  What happens in this case?
#endif // RYTHMOS_DEBUG
  time_vec[0] = 2.0;
  x_vec[0] = Teuchos::null;
  // x_vec pointers are null
#ifdef RYTHMOS_DEBUG
  TEST_THROW( ib->addPoints(time_vec,x_vec,xdot_vec), std::logic_error );
#else // RYTHMOS_DEBUG
  // TODO:  What happens in this case?
#endif // RYTHMOS_DEBUG
  x_vec[0] = createDefaultVector(2,17.0);
  xdot_vec[0] = Teuchos::null;
  // xdot_vec pointers are null (this is okay)
  ib->addPoints(time_vec,x_vec,xdot_vec);
  time_vec[0] = 2.1;
  xdot_vec[0] = createDefaultVector(2,32.0);
  RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
  pl->set("InterpolationBufferPolicy", "Static Policy");
  pl->set("StorageLimit",N+1);
  ib->setParameterList(pl);
  // Exceed storage due to IBPolicy == static 
  TEST_THROW( ib->addPoints(time_vec,x_vec,xdot_vec), std::logic_error );
  ib->setStorage(ib->getStorage()+1);
  // Verify we can add the point if the storage limit is one higher
  TEST_NOTHROW( ib->addPoints(time_vec,x_vec,xdot_vec) );

  ib->setStorage(ib->getStorage()+2);
  time_vec[0] = 2.5;
  time_vec.push_back(2.25);
  x_vec.push_back(createDefaultVector(2,18.0));
  xdot_vec.push_back(createDefaultVector(2,33.0));
  // addPoints with non-sorted time_vec
#ifdef RYTHMOS_DEBUG
  TEST_THROW( ib->addPoints(time_vec,x_vec,xdot_vec), std::logic_error );
#else // RYTHMOS_DEBUG
  // TODO:  What happens in this case?
#endif // RYTHMOS_DEBUG
  time_vec[1] = 2.5;
  // addPoints with non-unique time_vec
#ifdef RYTHMOS_DEBUG
  TEST_THROW( ib->addPoints(time_vec,x_vec,xdot_vec), std::logic_error );
#else // RYTHMOS_DEBUG
  // TODO:  What happens in this case?
#endif // RYTHMOS_DEBUG
  // Verify we can add this points with valid time_vec data.
  time_vec[1] = 2.75;
  TEST_NOTHROW( ib->addPoints(time_vec,x_vec,xdot_vec));
  // Double check we know what the timerange should be now.
  TimeRange<double> tr = ib->getTimeRange();
  TEST_EQUALITY_CONST( tr.lower(), 0.5 );
  TEST_EQUALITY_CONST( tr.upper(), 2.75 );
}

TEUCHOS_UNIT_TEST( Rythmos_InterpolationBuffer, addPoints_IBPolicy ) {
  int N = 5;
  Array<double> time_vec;
  Array<RCP<const VectorBase<double> > > x_vec;
  Array<RCP<const VectorBase<double> > > xdot_vec;
  for (int i=0 ; i<N ; ++i) {
    RCP<VectorBase<double> > v = createDefaultVector(2,1.0+i*1.0);
    RCP<VectorBase<double> > vdot = createDefaultVector(2,2.0+i*1.0);
    time_vec.push_back(10.0+i*1.0);
    x_vec.push_back(v);
    xdot_vec.push_back(vdot);
  }
  {
    // Default Policy is keep newest
    RCP<InterpolationBuffer<double> > ib = interpolationBuffer<double>();
    ib->setStorage(N);
    ib->addPoints(
        time_vec,
        x_vec,
        xdot_vec
        );
    TimeRange<double> tr = ib->getTimeRange();
    TEST_EQUALITY_CONST( tr.lower(), 10.0 );
    TEST_EQUALITY_CONST( tr.upper(), 14.0 );
    TEST_EQUALITY_CONST( ib->getIBPolicy(), BUFFER_POLICY_KEEP_NEWEST );
    // Case 1.  new points are at end
    {
      //   Case 1a.  one space left, add two, removes one & adds two
      ib->setStorage(ib->getStorage()+1);
      Array<double> time_vec_new;
      Array<RCP<const VectorBase<double> > > x_vec_new;
      Array<RCP<const VectorBase<double> > > xdot_vec_new;
      time_vec_new.push_back(15.0);
      x_vec_new.push_back(createDefaultVector(2,6.0));
      xdot_vec_new.push_back(createDefaultVector(2,7.0));
      time_vec_new.push_back(16.0);
      x_vec_new.push_back(createDefaultVector(2,7.0));
      xdot_vec_new.push_back(createDefaultVector(2,8.0));
      ib->addPoints(
          time_vec_new,
          x_vec_new,
          xdot_vec_new
          );
      tr = ib->getTimeRange();
      TEST_EQUALITY_CONST( tr.lower(), 11.0 );
      TEST_EQUALITY_CONST( tr.upper(), 16.0 );
    }
    {
      //   Case 1b.  no space left, add two, removes two & adds two
      Array<double> time_vec_new;
      Array<RCP<const VectorBase<double> > > x_vec_new;
      Array<RCP<const VectorBase<double> > > xdot_vec_new;
      time_vec_new.push_back(17.0);
      x_vec_new.push_back(createDefaultVector(2,7.0));
      xdot_vec_new.push_back(createDefaultVector(2,8.0));
      time_vec_new.push_back(18.0);
      x_vec_new.push_back(createDefaultVector(2,8.0));
      xdot_vec_new.push_back(createDefaultVector(2,9.0));
      ib->addPoints(
          time_vec_new,
          x_vec_new,
          xdot_vec_new
          );
      tr = ib->getTimeRange();
      TEST_EQUALITY_CONST( tr.lower(), 13.0 );
      TEST_EQUALITY_CONST( tr.upper(), 18.0 );
    }
    // Case 2.  new points are at beginning
    {
      //   Case 1a.  one space left, add two, removes one & adds two
      ib->setStorage(ib->getStorage()+1);
      Array<double> time_vec_new;
      Array<RCP<const VectorBase<double> > > x_vec_new;
      Array<RCP<const VectorBase<double> > > xdot_vec_new;
      time_vec_new.push_back(8.0);
      x_vec_new.push_back(createDefaultVector(2,8.0));
      xdot_vec_new.push_back(createDefaultVector(2,9.0));
      time_vec_new.push_back(9.0);
      x_vec_new.push_back(createDefaultVector(2,9.0));
      xdot_vec_new.push_back(createDefaultVector(2,10.0));
      ib->addPoints(
          time_vec_new,
          x_vec_new,
          xdot_vec_new
          );
      tr = ib->getTimeRange();
      TEST_EQUALITY_CONST( tr.lower(), 8.0 );
      TEST_EQUALITY_CONST( tr.upper(), 17.0 );
    }
    {
      //   Case 1b.  no space left, add two, removes two & adds two
      Array<double> time_vec_new;
      Array<RCP<const VectorBase<double> > > x_vec_new;
      Array<RCP<const VectorBase<double> > > xdot_vec_new;
      time_vec_new.push_back(6.0);
      x_vec_new.push_back(createDefaultVector(2,8.0));
      xdot_vec_new.push_back(createDefaultVector(2,9.0));
      time_vec_new.push_back(7.0);
      x_vec_new.push_back(createDefaultVector(2,9.0));
      xdot_vec_new.push_back(createDefaultVector(2,10.0));
      ib->addPoints(
          time_vec_new,
          x_vec_new,
          xdot_vec_new
          );
      tr = ib->getTimeRange();
      TEST_EQUALITY_CONST( tr.lower(), 6.0 );
      TEST_EQUALITY_CONST( tr.upper(), 15.0 );
    }
    {
      // Case 3.  new points are both at beginning and end, throw
      Array<double> time_vec_new;
      Array<RCP<const VectorBase<double> > > x_vec_new;
      Array<RCP<const VectorBase<double> > > xdot_vec_new;
      time_vec_new.push_back(5.0);
      x_vec_new.push_back(createDefaultVector(2,8.0));
      xdot_vec_new.push_back(createDefaultVector(2,9.0));
      time_vec_new.push_back(16.0);
      x_vec_new.push_back(createDefaultVector(2,9.0));
      xdot_vec_new.push_back(createDefaultVector(2,10.0));
      TEST_THROW( ib->addPoints(time_vec_new,x_vec_new,xdot_vec_new), std::logic_error );
    }
  }

  {
    // Policy = static
    RCP<InterpolationBuffer<double> > ib = interpolationBuffer<double>();
    RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    pl->set("InterpolationBufferPolicy", "Static Policy");
    ib->setParameterList(pl);
    ib->setStorage(N);
    ib->addPoints(
        time_vec,
        x_vec,
        xdot_vec
        );
    TimeRange<double> tr = ib->getTimeRange();
    TEST_EQUALITY_CONST( tr.lower(), 10.0 );
    TEST_EQUALITY_CONST( tr.upper(), 14.0 );
    TEST_EQUALITY_CONST( ib->getIBPolicy(), BUFFER_POLICY_STATIC );
    {
      // Case 1.  one space left, add one, okay.
      ib->setStorage(ib->getStorage()+1);
      Array<double> time_vec_new;
      Array<RCP<const VectorBase<double> > > x_vec_new;
      Array<RCP<const VectorBase<double> > > xdot_vec_new;
      time_vec_new.push_back(15.0);
      x_vec_new.push_back(createDefaultVector(2,10.0));
      xdot_vec_new.push_back(createDefaultVector(2,20.0));
      TEST_NOTHROW( ib->addPoints(time_vec_new,x_vec_new,xdot_vec_new) );
    }
    {
      // Case 2.  no space left, add one, throws.
      Array<double> time_vec_new;
      Array<RCP<const VectorBase<double> > > x_vec_new;
      Array<RCP<const VectorBase<double> > > xdot_vec_new;
      time_vec_new.push_back(16.0);
      x_vec_new.push_back(createDefaultVector(2,11.0));
      xdot_vec_new.push_back(createDefaultVector(2,21.0));
      TEST_THROW( ib->addPoints(time_vec_new,x_vec_new,xdot_vec_new), std::logic_error );
    }
    {
      // Case 3.  no space left, add two, throws.
      Array<double> time_vec_new;
      Array<RCP<const VectorBase<double> > > x_vec_new;
      Array<RCP<const VectorBase<double> > > xdot_vec_new;
      time_vec_new.push_back(16.0);
      x_vec_new.push_back(createDefaultVector(2,11.0));
      xdot_vec_new.push_back(createDefaultVector(2,21.0));
      time_vec_new.push_back(17.0);
      x_vec_new.push_back(createDefaultVector(2,12.0));
      xdot_vec_new.push_back(createDefaultVector(2,22.0));
      TEST_THROW( ib->addPoints(time_vec_new,x_vec_new,xdot_vec_new), std::logic_error );
    }
    {
      // Case 4.  one space left, adding two, throws.
      ib->setStorage(ib->getStorage()+1);
      Array<double> time_vec_new;
      Array<RCP<const VectorBase<double> > > x_vec_new;
      Array<RCP<const VectorBase<double> > > xdot_vec_new;
      time_vec_new.push_back(16.0);
      x_vec_new.push_back(createDefaultVector(2,11.0));
      xdot_vec_new.push_back(createDefaultVector(2,21.0));
      time_vec_new.push_back(17.0);
      x_vec_new.push_back(createDefaultVector(2,12.0));
      xdot_vec_new.push_back(createDefaultVector(2,22.0));
      TEST_THROW( ib->addPoints(time_vec_new,x_vec_new,xdot_vec_new), std::logic_error );
    }
  }
}

TEUCHOS_UNIT_TEST( Rythmos_InterpolationBuffer, getPoints ) {
  {
    int N = 5;
    Array<double> time_vec;
    Array<RCP<const VectorBase<double> > > x_vec;
    Array<RCP<const VectorBase<double> > > xdot_vec;
    for (int i=0 ; i<N ; ++i) {
      RCP<VectorBase<double> > v = createDefaultVector(2,1.0+i*1.0);
      RCP<VectorBase<double> > vdot = createDefaultVector(2,2.0+i*1.0);
      time_vec.push_back(10.0+i*1.0);
      x_vec.push_back(v);
      xdot_vec.push_back(vdot);
    }
    RCP<InterpolationBuffer<double> > ib = interpolationBuffer<double>();
    ib->setStorage(N);
    ib->addPoints(
        time_vec,
        x_vec,
        xdot_vec
        );
    TimeRange<double> tr = ib->getTimeRange();
    TEST_EQUALITY_CONST( tr.lower(), 10.0 );
    TEST_EQUALITY_CONST( tr.upper(), 14.0 );
  }

  // Case A.  1 point in buffer
  {
    int N = 1;
    Array<double> time_vec;
    Array<RCP<const VectorBase<double> > > x_vec;
    Array<RCP<const VectorBase<double> > > xdot_vec;
    for (int i=0 ; i<N ; ++i) {
      RCP<VectorBase<double> > v = createDefaultVector(2,1.0+i*1.0);
      RCP<VectorBase<double> > vdot = createDefaultVector(2,2.0+i*1.0);
      time_vec.push_back(10.0+i*1.0);
      x_vec.push_back(v);
      xdot_vec.push_back(vdot);
    }
    RCP<InterpolationBuffer<double> > ib = interpolationBuffer<double>();
    ib->setStorage(N);
    ib->addPoints(
        time_vec,
        x_vec,
        xdot_vec
        );
    TimeRange<double> tr = ib->getTimeRange();
    TEST_EQUALITY_CONST( tr.lower(), 10.0 );
    TEST_EQUALITY_CONST( tr.upper(), 10.0 );
    {
      //   Case A.1.  asking for anything but one point throws
      Array<double> time_vec_out;
      Array<RCP<const VectorBase<double> > > x_vec_out;
      Array<RCP<const VectorBase<double> > > xdot_vec_out;
      Array<double> accuracy_vec_out;
      time_vec_out.push_back(9.0);
#ifdef RYTHMOS_DEBUG
      TEST_THROW(
          ib->getPoints(time_vec_out, &x_vec_out, &xdot_vec_out, &accuracy_vec_out), 
          std::logic_error
          );
#else // RYTHMOS_DEBUG
      TEST_NOTHROW(
          ib->getPoints(time_vec_out, &x_vec_out, &xdot_vec_out, &accuracy_vec_out)
          );
      // TODO:  What do we get out in this case?
#endif // RYTHMOS_DEBUG
      time_vec_out[0] = 11.0;
#ifdef RYTHMOS_DEBUG
      TEST_THROW(
          ib->getPoints(time_vec_out, &x_vec_out, &xdot_vec_out, &accuracy_vec_out), 
          std::logic_error
          );
#else // RYTHMOS_DEBUG
      TEST_NOTHROW(
          ib->getPoints(time_vec_out, &x_vec_out, &xdot_vec_out, &accuracy_vec_out)
          );
      // TODO:  What do we get out in this case?
#endif // RYTHMOS_DEBUG
      //   Case A.2.  asking for exact time, okay
      time_vec_out[0] = 10.0;
      TEST_NOTHROW(
          ib->getPoints(time_vec_out, &x_vec_out, &xdot_vec_out, &accuracy_vec_out)
          );
      {
        Thyra::ConstDetachedVectorView<double> x_view( *x_vec_out[0] );
        TEST_EQUALITY_CONST( x_view[0], 1.0 );
        TEST_EQUALITY_CONST( x_view[1], 1.0 );
        Thyra::ConstDetachedVectorView<double> xdot_view( *xdot_vec_out[0] );
        TEST_EQUALITY_CONST( xdot_view[0], 2.0 );
        TEST_EQUALITY_CONST( xdot_view[1], 2.0 );
      }
      //   Case A.3  asking for multiple points, throws
      time_vec_out.push_back(10.1);
      TEST_THROW(
          ib->getPoints(time_vec_out, &x_vec_out, &xdot_vec_out, &accuracy_vec_out), 
          std::logic_error
          );
    }
  }
  {
    // Case B.  >=2 points in buffer
    int N = 3;
    Array<double> time_vec;
    Array<RCP<const VectorBase<double> > > x_vec;
    Array<RCP<const VectorBase<double> > > xdot_vec;
    for (int i=0 ; i<N ; ++i) {
      RCP<VectorBase<double> > v = createDefaultVector(2,1.0+i*1.0);
      RCP<VectorBase<double> > vdot = createDefaultVector(2,2.0+i*1.0);
      time_vec.push_back(10.0+i*1.0);
      x_vec.push_back(v);
      xdot_vec.push_back(vdot);
    }
    RCP<InterpolationBuffer<double> > ib = interpolationBuffer<double>();
    ib->setStorage(N);
    ib->addPoints(
        time_vec,
        x_vec,
        xdot_vec
        );
    TimeRange<double> tr = ib->getTimeRange();
    TEST_EQUALITY_CONST( tr.lower(), 10.0 );
    TEST_EQUALITY_CONST( tr.upper(), 12.0 );
    {
      //   Case B.1.  one point before TimeRange, throws
      Array<double> time_vec_out;
      Array<RCP<const VectorBase<double> > > x_vec_out;
      Array<RCP<const VectorBase<double> > > xdot_vec_out;
      Array<double> accuracy_vec_out;
      time_vec_out.push_back(9.0);
      TEST_THROW(
          ib->getPoints(time_vec_out, &x_vec_out, &xdot_vec_out, &accuracy_vec_out), 
          std::logic_error
          );
      //   Case B.2.  one point after TimeRange, throws
      time_vec_out[0] = 13.0;
      TEST_THROW(
          ib->getPoints(time_vec_out, &x_vec_out, &xdot_vec_out, &accuracy_vec_out), 
          std::logic_error
          );
      //   Case B.3.  multiple points, with some after TimeRange, throws
      time_vec_out[0] = 11.5;
      time_vec_out.push_back(12.5);
      TEST_THROW(
          ib->getPoints(time_vec_out, &x_vec_out, &xdot_vec_out, &accuracy_vec_out), 
          std::logic_error
          );
      //   Case B.4.  multiple points, with some before TimeRange, throws
      time_vec_out[0] = 9.5;
      time_vec_out[1] = 10.5;
      TEST_THROW(
          ib->getPoints(time_vec_out, &x_vec_out, &xdot_vec_out, &accuracy_vec_out), 
          std::logic_error
          );
      //   Case B.5.  multiple points, with some before and after TimeRange, throws
      time_vec_out[0] = 9.5;
      time_vec_out[1] = 10.5;
      time_vec_out.push_back(12.5);
      TEST_THROW(
          ib->getPoints(time_vec_out, &x_vec_out, &xdot_vec_out, &accuracy_vec_out), 
          std::logic_error
          );
      //   Case B.6.  nodes, return exact value
      time_vec_out.clear();
      time_vec_out.push_back(11.0);
      TEST_NOTHROW(
          ib->getPoints(time_vec_out, &x_vec_out, &xdot_vec_out, &accuracy_vec_out)
          );
      {
        Thyra::ConstDetachedVectorView<double> x_view( *x_vec_out[0] );
        TEST_EQUALITY_CONST( x_view[0], 2.0 );
        TEST_EQUALITY_CONST( x_view[1], 2.0 );
        Thyra::ConstDetachedVectorView<double> xdot_view( *xdot_vec_out[0] );
        TEST_EQUALITY_CONST( xdot_view[0], 3.0 );
        TEST_EQUALITY_CONST( xdot_view[1], 3.0 );
      }
      //   Case B.7.  in between nodes, interpolated values
      //     In this case, the interpolator is linear (default)
      time_vec_out.clear();
      time_vec_out.push_back(10.5);
      time_vec_out.push_back(11.5);
      TEST_NOTHROW(
          ib->getPoints(time_vec_out, &x_vec_out, &xdot_vec_out, &accuracy_vec_out)
          );
      {
        Thyra::ConstDetachedVectorView<double> x_view0( *x_vec_out[0] );
        TEST_EQUALITY_CONST( x_view0[0], 1.5 );
        TEST_EQUALITY_CONST( x_view0[1], 1.5 );
        Thyra::ConstDetachedVectorView<double> xdot_view0( *xdot_vec_out[0] );
        TEST_EQUALITY_CONST( xdot_view0[0], 2.5 );
        TEST_EQUALITY_CONST( xdot_view0[1], 2.5 );

        Thyra::ConstDetachedVectorView<double> x_view1( *x_vec_out[1] );
        TEST_EQUALITY_CONST( x_view1[0], 2.5 );
        TEST_EQUALITY_CONST( x_view1[1], 2.5 );
        Thyra::ConstDetachedVectorView<double> xdot_view1( *xdot_vec_out[1] );
        TEST_EQUALITY_CONST( xdot_view1[0], 3.5 );
        TEST_EQUALITY_CONST( xdot_view1[1], 3.5 );
      }
    }
  }
}

TEUCHOS_UNIT_TEST( Rythmos_InterpolationBuffer, description ) {
  RCP<InterpolationBuffer<double> > ib = interpolationBuffer<double>();
  std::string desc = ib->description();
  TEST_EQUALITY_CONST( desc, "Rythmos::InterpolationBuffer" );
}

//TEUCHOS_UNIT_TEST( Rythmos_InterpolationBuffer, describe ) {
//  RCP<InterpolationBuffer<double> > ib = interpolationBuffer<double>();
//}

TEUCHOS_UNIT_TEST( Rythmos_InterpolationBuffer, add_get_points_1 ) {
  RCP<InterpolationBuffer<double> > ib = interpolationBuffer<double>();
  double t = 1.0;
  RCP<VectorBase<double> > x = createDefaultVector(2,1.0);
  RCP<VectorBase<double> > xdot = createDefaultVector(2,1.0);
  {
    Array<double> time_vec;
    Array<RCP<const VectorBase<double> > > x_vec;
    Array<RCP<const VectorBase<double> > > xdot_vec;
    time_vec.push_back(t);
    x_vec.push_back(x);
    xdot_vec.push_back(xdot);
    ib->addPoints(time_vec,x_vec,xdot_vec);
  }

  {
    Array<double> time_vec;
    Array<RCP<const VectorBase<double> > > x_vec;
    Array<RCP<const VectorBase<double> > > xdot_vec;
    Array<double> accuracy_vec;
    time_vec.push_back(t);
    ib->getPoints(time_vec,&x_vec,&xdot_vec,&accuracy_vec);
    RCP<const VectorBase<double> > x0 = x_vec[0];
    RCP<VectorBase<double> > y = x0->clone_v();
    TEST_ASSERT( !is_null(y) );
  }
}
TEUCHOS_UNIT_TEST( Rythmos_InterpolationBuffer, add_get_points_2 ) {
  RCP<InterpolationBuffer<double> > ib = interpolationBuffer<double>();
  double t1 = 0.0;
  RCP<VectorBase<double> > x1 = createDefaultVector(2,1.0);
  RCP<VectorBase<double> > xdot1 = createDefaultVector(2,2.0);
  double t2 = 1.0;
  RCP<VectorBase<double> > x2 = createDefaultVector(2,3.0);
  RCP<VectorBase<double> > xdot2 = createDefaultVector(2,4.0);
  {
    Array<double> time_vec;
    Array<RCP<const VectorBase<double> > > x_vec;
    Array<RCP<const VectorBase<double> > > xdot_vec;
    time_vec.push_back(t1);
    time_vec.push_back(t2);
    x_vec.push_back(x1);
    x_vec.push_back(x2);
    xdot_vec.push_back(xdot1);
    xdot_vec.push_back(xdot2);
    ib->addPoints(time_vec,x_vec,xdot_vec);
  }

  {
    Array<double> time_vec;
    Array<RCP<const VectorBase<double> > > x_vec;
    Array<RCP<const VectorBase<double> > > xdot_vec;
    Array<double> accuracy_vec;
    time_vec.push_back(t1);
    time_vec.push_back(t2);
    ib->getPoints(time_vec,&x_vec,&xdot_vec,&accuracy_vec);
    RCP<VectorBase<double> > y = x_vec[0]->clone_v();
    TEST_ASSERT( !is_null(y) );
  }
}

#ifdef RYTHMOS_BROKEN_TEST
// BUG 4388
TEUCHOS_UNIT_TEST( Rythmos_InterpolationBuffer, add_to_empty ) {
  RCP<InterpolationBuffer<double> > ibSource = interpolationBuffer<double>();
  int dim = 1;
  {
    Teuchos::Array<double> time_vec;
    Array<RCP<const VectorBase<double> > > x_vec;
    Array<RCP<const VectorBase<double> > > xdot_vec;
    time_vec.push_back(0.0);
    x_vec.push_back(createDefaultVector(dim,1.0)); 
    xdot_vec.push_back(createDefaultVector(dim,2.0)); 
    ibSource->addPoints(time_vec,x_vec,xdot_vec);
  }
  {
    Array<double> time_vec_out;
    time_vec_out.push_back(0.0);
    Array<RCP<const Thyra::VectorBase<double> > > x_vec_out;
    Array<RCP<const Thyra::VectorBase<double> > > xdot_vec_out;
    Array<double> accuracy_vec_out;
    ibSource->getPoints(time_vec_out, &x_vec_out, &xdot_vec_out, &accuracy_vec_out);
    TEST_ASSERT( x_vec_out.length() == 1 );
    RCP<VectorBase<double> > xnew = x_vec_out[0]->clone_v();
    TEST_ASSERT(!is_null(xnew));
  }
}
#endif // RYTHMOS_BROKEN_TEST

} // namespace Rythmos



