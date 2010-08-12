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

#include "Rythmos_DataStore.hpp"
#include "Rythmos_UnitTestHelpers.hpp"

namespace Rythmos {

TEUCHOS_UNIT_TEST( Rythmos_DataStore, newDataStore ) {
  DataStore<double> ds;
  TEST_EQUALITY_CONST( ds.time, -1 );
  TEST_EQUALITY_CONST( Teuchos::is_null(ds.x), true );
  TEST_EQUALITY_CONST( Teuchos::is_null(ds.xdot), true );
  TEST_EQUALITY_CONST( ds.accuracy, -1 );
}

TEUCHOS_UNIT_TEST( Rythmos_DataStore, newDataStore2 ) {
  double time = 1.0234;
  double accuracy = 0.0012;
  RCP<Thyra::VectorBase<double> > x = createDefaultVector<double>(5,1.0);
  RCP<Thyra::VectorBase<double> > xdot = createDefaultVector<double>(5,2.0);
  DataStore<double> ds(time,x,xdot,accuracy);
  TEST_EQUALITY( time, ds.time );
  TEST_EQUALITY( accuracy, ds.accuracy );
  TEST_EQUALITY( x, ds.x ); // check that pointers are equal
  TEST_EQUALITY( xdot, ds.xdot );
}

TEUCHOS_UNIT_TEST( Rythmos_DataStore, constructor ) {
  double time = 0.025;
  double accuracy = .0000025;
  DataStore<double> ds(time,Teuchos::null,Teuchos::null,accuracy);
  DataStore<double> ds2(ds);
  TEST_EQUALITY( ds2.time, time );
  TEST_EQUALITY( ds2.accuracy, accuracy );
  // This is a shallow copy constructor, use clone for a deep copy.
}

TEUCHOS_UNIT_TEST( Rythmos_DataStore, lessthan ) {
  double timeA = 1.0;
  double timeB = 2.0;
  double accuracy = 0.0;
  DataStore<double> dsA(timeA,Teuchos::null,Teuchos::null,accuracy);
  DataStore<double> dsB(timeB,Teuchos::null,Teuchos::null,accuracy);
  TEST_COMPARE( dsA, < , dsB );
  TEST_COMPARE( dsA, <=, dsB );
  TEST_COMPARE( dsB, > , dsA );
  TEST_COMPARE( dsB, >=, dsA );
  TEST_COMPARE( dsA, < , 5.0 );
  TEST_COMPARE( dsA, <=, 5.0 );
  TEST_COMPARE( dsA, > , 0.5 );
  TEST_COMPARE( dsA, >=, 0.5 );
}

TEUCHOS_UNIT_TEST( Rythmos_DataStore, lessthanequal ) {
  double timeA = 1.0;
  double timeB = 1.0;
  double accuracy = 0.0;
  DataStore<double> dsA(timeA,Teuchos::null,Teuchos::null,accuracy);
  DataStore<double> dsB(timeB,Teuchos::null,Teuchos::null,accuracy);
  TEST_COMPARE( dsA, <=, dsB );
  TEST_COMPARE( dsA, ==, dsB );
  TEST_COMPARE( dsB, >=, dsA );
  TEST_COMPARE( dsB, ==, dsA );
  TEST_COMPARE( dsA, ==, 1.0 );
}

TEUCHOS_UNIT_TEST( Rythmos_DataStore, dataStoreVectorToVector ) {
  DataStore<double>::DataStoreVector_t ds_v;
  Array<double> time_vec_gs;
  Array<double> accuracy_vec_gs;
  int N = 3;
  for (int i=0 ; i<N ; ++i) {
    double time = 1.0*i + 1.0;
    double accuracy = 1.0*i + 2.0;
    DataStore<double> ds(time,Teuchos::null,Teuchos::null,accuracy);
    ds_v.push_back(ds);
    time_vec_gs.push_back(time);
    accuracy_vec_gs.push_back(accuracy);
  }

  Array<double> time_vec;
  Array<RCP<const Thyra::VectorBase<double> > > x_vec;
  Array<RCP<const Thyra::VectorBase<double> > > xdot_vec;
  Array<double> accuracy_vec;
  // Make sure the output is cleared
  for (int i=0 ; i<2*N ; ++i) {
    time_vec.push_back(0.0);
    x_vec.push_back(Teuchos::null);
    xdot_vec.push_back(Teuchos::null);
    accuracy_vec.push_back(0.0);
  }
  dataStoreVectorToVector(ds_v,&time_vec,&x_vec,&xdot_vec,&accuracy_vec);

  using Teuchos::as;
  TEST_EQUALITY( as<int>(time_vec.size()), N );
  TEST_EQUALITY( as<int>(x_vec.size()), N );
  TEST_EQUALITY( as<int>(xdot_vec.size()), N );
  TEST_EQUALITY( as<int>(accuracy_vec.size()), N );
  TEST_COMPARE_ARRAYS( time_vec, time_vec_gs );
  TEST_COMPARE_ARRAYS( accuracy_vec, accuracy_vec_gs );
}

TEUCHOS_UNIT_TEST( Rythmos_DataStore, vectorToDataStoreVector ) {
  Array<double> time_vec;
  Array<RCP<const Thyra::VectorBase<double> > > x_vec;
  Array<RCP<const Thyra::VectorBase<double> > > xdot_vec;
  Array<double> accuracy_vec;

  int N = 3;
  for (int i=0 ; i<N ; ++i) {
    time_vec.push_back(1.0*i+1.0);
    x_vec.push_back(Teuchos::null);
    xdot_vec.push_back(Teuchos::null);
    accuracy_vec.push_back(2.0*i+2.0);
  }

  DataStore<double>::DataStoreVector_t ds_v;
  // Make sure the output is cleared
  for (int i=0 ; i<2*N ; ++i) {
    DataStore<double> ds;
    ds_v.push_back(ds);
  }
  vectorToDataStoreVector(time_vec,x_vec,xdot_vec,accuracy_vec,&ds_v);

  using Teuchos::as;
  TEST_EQUALITY( as<int>(ds_v.size()), N );
  for (int i=0 ; i<N ; ++i) {
    TEST_EQUALITY( ds_v[i].time, time_vec[i] );
    TEST_EQUALITY( ds_v[i].accuracy, accuracy_vec[i] );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_DataStore, vectorToDataStoreList ) {
  Array<double> time_vec;
  Array<RCP<const Thyra::VectorBase<double> > > x_vec;
  Array<RCP<const Thyra::VectorBase<double> > > xdot_vec;
  Array<double> accuracy_vec;

  int N = 3;
  for (int i=0 ; i<N ; ++i) {
    time_vec.push_back(1.0*i+1.0);
    x_vec.push_back(Teuchos::null);
    xdot_vec.push_back(Teuchos::null);
    accuracy_vec.push_back(2.0*i+2.0);
  }

  DataStore<double>::DataStoreList_t ds_l;
  // Make sure the output is cleared
  for (int i=0 ; i<2*N; ++i) {
    DataStore<double> ds;
    ds_l.push_back(ds);
  }
  vectorToDataStoreList(time_vec,x_vec,xdot_vec,accuracy_vec,&ds_l);

  using Teuchos::as;
  TEST_EQUALITY( as<int>(ds_l.size()), N );
  DataStore<double>::DataStoreList_t::iterator it;
  int i=0;
  for (it = ds_l.begin() ; it != ds_l.end() ; it++) {
    TEST_EQUALITY( it->time, time_vec[i] );
    TEST_EQUALITY( it->accuracy, accuracy_vec[i] );
    ++i;
  }
}

TEUCHOS_UNIT_TEST( Rythmos_DataStore, vectorToDataStoreListNoAccuracy ) {
  Array<double> time_vec;
  Array<RCP<const Thyra::VectorBase<double> > > x_vec;
  Array<RCP<const Thyra::VectorBase<double> > > xdot_vec;

  int N = 3;
  for (int i=0 ; i<N ; ++i) {
    time_vec.push_back(1.0*i+1.0);
    x_vec.push_back(Teuchos::null);
    xdot_vec.push_back(Teuchos::null);
  }

  DataStore<double>::DataStoreList_t ds_l;
  // Make sure the output is cleared
  for (int i=0 ; i<2*N; ++i) {
    DataStore<double> ds;
    ds_l.push_back(ds);
  }
  vectorToDataStoreList(time_vec,x_vec,xdot_vec,&ds_l);

  using Teuchos::as;
  TEST_EQUALITY( as<int>(ds_l.size()), N );
  DataStore<double>::DataStoreList_t::iterator it;
  int i=0;
  for (it = ds_l.begin() ; it != ds_l.end() ; it++) {
    TEST_EQUALITY( it->time, time_vec[i] );
    ++i;
  }
}

TEUCHOS_UNIT_TEST( Rythmos_DataStore, clone ) {
  {
    // Normal clone
    double t = 2.0;
    RCP<VectorBase<double> > v = createDefaultVector(2,3.0);
    RCP<VectorBase<double> > v_dot = createDefaultVector(2,4.0);
    double accuracy = 0.5;
    const DataStore<double> ds(t,v,v_dot,accuracy);

    RCP<DataStore<double> > ds_clone_ptr = ds.clone();
    DataStore<double>& ds_clone = *ds_clone_ptr;

    Thyra::V_S(v.ptr(),5.0);
    Thyra::V_S(v_dot.ptr(),6.0);
    TEST_EQUALITY_CONST(Thyra::get_ele(*(ds_clone.x),0), 3.0);
    TEST_EQUALITY_CONST(Thyra::get_ele(*(ds_clone.xdot),0), 4.0);
  }
  {
    // Clone with missing xdot
    double t = 2.0;
    RCP<VectorBase<double> > v = createDefaultVector(2,3.0);
    RCP<VectorBase<double> > v_dot; 
    double accuracy = 0.5;
    const DataStore<double> ds(t,v,v_dot,accuracy);

    RCP<DataStore<double> > ds_clone_ptr = ds.clone();
    DataStore<double>& ds_clone = *ds_clone_ptr;

    Thyra::V_S(v.ptr(),5.0);
    TEST_EQUALITY_CONST(Thyra::get_ele(*(ds_clone.x),0), 3.0);
    TEST_EQUALITY_CONST( Teuchos::is_null(ds_clone.xdot), true );
  }
  {
    // Clone with missing x
    double t = 2.0;
    RCP<VectorBase<double> > v;
    RCP<VectorBase<double> > v_dot = createDefaultVector(2,4.0);
    double accuracy = 0.5;
    const DataStore<double> ds(t,v,v_dot,accuracy);

    RCP<DataStore<double> > ds_clone_ptr = ds.clone();
    DataStore<double>& ds_clone = *ds_clone_ptr;

    Thyra::V_S(v_dot.ptr(),6.0);
    TEST_EQUALITY_CONST(Teuchos::is_null(ds_clone.x), true );
    TEST_EQUALITY_CONST(Thyra::get_ele(*(ds_clone.xdot),0), 4.0);
  }


}

TEUCHOS_UNIT_TEST( Rythmos_DataStore, shallowCopy ) {
  double t = 1.0;
  RCP<VectorBase<double> > x = createDefaultVector(2,1.0);
  RCP<VectorBase<double> > xdot = createDefaultVector(2,2.0);
  double accuracy = 0.0;
  DataStore<double> DS(t,x,xdot,accuracy);
  DataStore<double> newDS(DS);
  {
    RCP<const VectorBase<double> > newX = newDS.x;
    RCP<VectorBase<double> > newXclone = newX->clone_v();
    TEST_ASSERT( !is_null(newXclone) );
    TEST_EQUALITY( DS.x.get(), newDS.x.get() );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_DataStore, dataStoreVectorToVector_more ) {
  double time = 1.0;
  RCP<VectorBase<double> > x = createDefaultVector(2,1.0);
  RCP<VectorBase<double> > xdot = createDefaultVector(2,2.0);
  double accuracy = 0.0;
  DataStore<double> DS(time,x,xdot,accuracy);
  DataStore<double>::DataStoreVector_t data_out;
  data_out.push_back(DS);
  Array<double> time_vec;
  Array<RCP<const VectorBase<double> > > x_vec;
  Array<RCP<const VectorBase<double> > > xdot_vec;
  Array<double> accuracy_vec;
  dataStoreVectorToVector<double>(data_out,&time_vec,&x_vec,&xdot_vec,&accuracy_vec);
  TEST_EQUALITY( time_vec[0], time );
  TEST_EQUALITY( x_vec[0].get(), x.get() );
  TEST_EQUALITY( xdot_vec[0].get(), xdot.get() );
  TEST_EQUALITY( accuracy_vec[0], accuracy );
  RCP<VectorBase<double> > newX = x_vec[0]->clone_v();
  TEST_ASSERT( !is_null(newX) );
}


} // namespace Rythmos


