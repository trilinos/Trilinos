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

#ifndef Rythmos_DATA_STORE_UNITTEST_H
#define Rythmos_DATA_STORE_UNITTEST_H

#include "../test/CppUnitLite/TestHarness.h"

#include "Rythmos_DataStore.hpp"

namespace Rythmos {

TEST( DataStore, newDataStore ) {
  DataStore<double> ds;
  CHECK( -1 == ds.time );
  CHECK( Teuchos::is_null(ds.x) );
  CHECK( Teuchos::is_null(ds.xdot) );
  CHECK( -1 == ds.accuracy );
}

TEST( DataStore, newDataStore2 ) {
  double time = 1.0234;
  double accuracy = 0.0012;
  DataStore<double> ds(time,Teuchos::null,Teuchos::null,accuracy);
  CHECK( time == ds.time );
  CHECK( accuracy == ds.accuracy );
}

TEST( DataStore, clone ) {
  double time = 0.025;
  double accuracy = .0000025;
  DataStore<double> ds(time,Teuchos::null,Teuchos::null,accuracy);
  DataStore<double> ds2(ds);
  CHECK( ds2.time == time );
  CHECK( ds2.accuracy == accuracy );
}

TEST( DataStore, lessthan ) {
  double timeA = 1.0;
  double timeB = 2.0;
  double accuracy = 0.0;
  DataStore<double> dsA(timeA,Teuchos::null,Teuchos::null,accuracy);
  DataStore<double> dsB(timeB,Teuchos::null,Teuchos::null,accuracy);
  CHECK( dsA < dsB );
  CHECK( dsA <= dsB );
  CHECK( dsB > dsA );
  CHECK( dsB >= dsA );
  CHECK( dsA < 5.0 );
  CHECK( dsA <= 5.0 );
  CHECK( dsA > 0.5 );
  CHECK( dsA >= 0.5 );
}

TEST( DataStore, lessthanequal ) {
  double timeA = 1.0;
  double timeB = 1.0;
  double accuracy = 0.0;
  DataStore<double> dsA(timeA,Teuchos::null,Teuchos::null,accuracy);
  DataStore<double> dsB(timeB,Teuchos::null,Teuchos::null,accuracy);
  CHECK( dsA <= dsB );
  CHECK( dsA == dsB );
  CHECK( dsB >= dsA );
  CHECK( dsB == dsA );
  CHECK( dsA == 1.0 );
}

TEST( DataStore, dataStoreVectorToVector ) {
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
  Array<Teuchos::RCP<const Thyra::VectorBase<double> > > x_vec;
  Array<Teuchos::RCP<const Thyra::VectorBase<double> > > xdot_vec;
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
  CHECK( as<int>(time_vec.size()) == N );
  CHECK( as<int>(x_vec.size()) == N );
  CHECK( as<int>(xdot_vec.size()) == N );
  CHECK( as<int>(accuracy_vec.size()) == N );
  for (int i=0 ; i<N ; ++i) {
    CHECK( time_vec[i] == time_vec_gs[i] );
    CHECK( accuracy_vec[i] == accuracy_vec_gs[i] );
  }
}

TEST( DataStore, vectorToDataStoreVector ) {
  Array<double> time_vec;
  Array<Teuchos::RCP<const Thyra::VectorBase<double> > > x_vec;
  Array<Teuchos::RCP<const Thyra::VectorBase<double> > > xdot_vec;
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
  CHECK( as<int>(ds_v.size()) == N );
  for (int i=0 ; i<N ; ++i) {
    CHECK( ds_v[i].time == time_vec[i] );
    CHECK( ds_v[i].accuracy == accuracy_vec[i] );
  }
}

TEST( DataStore, vectorToDataStoreList ) {
  Array<double> time_vec;
  Array<Teuchos::RCP<const Thyra::VectorBase<double> > > x_vec;
  Array<Teuchos::RCP<const Thyra::VectorBase<double> > > xdot_vec;
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
  CHECK( as<int>(ds_l.size()) == N );
  DataStore<double>::DataStoreList_t::iterator it;
  int i=0;
  for (it = ds_l.begin() ; it != ds_l.end() ; it++) {
    CHECK( it->time == time_vec[i] );
    CHECK( it->accuracy == accuracy_vec[i] );
    ++i;
  }
}

TEST( DataStore, vectorToDataStoreListNoAccuracy ) {
  Array<double> time_vec;
  Array<Teuchos::RCP<const Thyra::VectorBase<double> > > x_vec;
  Array<Teuchos::RCP<const Thyra::VectorBase<double> > > xdot_vec;

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
  CHECK( as<int>(ds_l.size()) == N );
  DataStore<double>::DataStoreList_t::iterator it;
  int i=0;
  for (it = ds_l.begin() ; it != ds_l.end() ; it++) {
    CHECK( it->time == time_vec[i] );
    ++i;
  }
}

} // namespace Rythmos

#endif // Rythmos_DATA_STORE_UNITTEST_H

