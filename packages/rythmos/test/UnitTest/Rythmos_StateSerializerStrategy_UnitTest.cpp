//@HEADER
// ***********************************************************************
//
//                           Rythmos Package
//                 Copyright (2009) Sandia Corporation
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
#include "Rythmos_StateSerializerStrategy.hpp"
#include "Rythmos_UnitTestHelpers.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "../SinCos/SinCosModel.hpp"

namespace Rythmos {

using Teuchos::outArg;

TEUCHOS_UNIT_TEST( Rythmos_StateSerializerStrategy, Scalars ) {
  std::string data;
  {
    std::ostringstream oss;
    double d = 4.2;
    double e = 4.3;
    double f = 4.4;
    XMLStateSerializerStrategy<double> ss;
    ss.serializeScalar(d,oss);
    ss.serializeScalar(e,oss);
    ss.serializeScalar(f,oss);
    data = oss.str();
  }
  {
    std::istringstream iss(data);
    double d = 0.0;
    XMLStateSerializerStrategy<double> ss;
    ss.deSerializeScalar(outArg(d),iss);
    TEST_EQUALITY_CONST(d,4.2);
    ss.deSerializeScalar(outArg(d),iss);
    TEST_EQUALITY_CONST(d,4.3);
    ss.deSerializeScalar(outArg(d),iss);
    TEST_EQUALITY_CONST(d,4.4);
  }
}

TEUCHOS_UNIT_TEST( Rythmos_StateSerializerStrategy, ints ) {
  std::string data;
  {
    std::ostringstream oss;
    int i = -1;
    int j = 3;
    int k = 600;
    XMLStateSerializerStrategy<double> ss;
    ss.serializeInt(i,oss);
    ss.serializeInt(j,oss);
    ss.serializeInt(k,oss);
    data = oss.str();
  }
  {
    std::istringstream iss(data);
    int i = 0;
    XMLStateSerializerStrategy<double> ss;
    ss.deSerializeInt(outArg(i),iss);
    TEST_EQUALITY_CONST(i,-1);
    ss.deSerializeInt(outArg(i),iss);
    TEST_EQUALITY_CONST(i,3);
    ss.deSerializeInt(outArg(i),iss);
    TEST_EQUALITY_CONST(i,600);
  }
}

TEUCHOS_UNIT_TEST( Rythmos_StateSerializerStrategy, ScalarsAndInts ) {
  std::string data;
  {
    std::ostringstream oss;
    int i = -1;
    double d = 31.4159e-1;
    int k = 600;
    XMLStateSerializerStrategy<double> ss;
    ss.serializeInt(i,oss);
    ss.serializeScalar(d,oss);
    ss.serializeInt(k,oss);
    data = oss.str();
  }
  {
    std::istringstream iss(data);
    int i = 0;
    double d = 0.0;
    XMLStateSerializerStrategy<double> ss;
    ss.deSerializeInt(outArg(i),iss);
    TEST_EQUALITY_CONST(i,-1);
    ss.deSerializeScalar(outArg(d),iss);
    TEST_EQUALITY_CONST(d,3.14159);
    ss.deSerializeInt(outArg(i),iss);
    TEST_EQUALITY_CONST(i,600);
  }
}

TEUCHOS_UNIT_TEST( Rythmos_StateSerializerStrategy, bools ) {
  std::string data;
  {
    std::ostringstream oss;
    bool a = true;
    XMLStateSerializerStrategy<double> ss;
    ss.serializeBool(a,oss);
    a = false;
    ss.serializeBool(a,oss);
    ss.serializeBool(a,oss);
    data = oss.str();
  }
  {
    std::istringstream iss(data);
    bool a = false;
    XMLStateSerializerStrategy<double> ss;
    ss.deSerializeBool(outArg(a),iss);
    TEST_EQUALITY_CONST(a,true);
    ss.deSerializeBool(outArg(a),iss);
    TEST_EQUALITY_CONST(a,false);
    ss.deSerializeBool(outArg(a),iss);
    TEST_EQUALITY_CONST(a,false);
  }
}


TEUCHOS_UNIT_TEST( Rythmos_StateSerializerStrategy, paramList ) {
  std::string data;
  Teuchos::ParameterList plGold;
  {
    plGold.set("Hello","World");
    plGold.set("alpha", 1.2345 );
    plGold.set("isTrue", false );
    plGold.set("NumTimes", 35 );
  }
  {
    std::ostringstream oss;
    XMLStateSerializerStrategy<double> ss;
    ss.serializeParameterList(plGold,oss);
    data = oss.str();
  }
  {
    std::istringstream iss(data);
    Teuchos::ParameterList pl;
    XMLStateSerializerStrategy<double> ss;
    ss.deSerializeParameterList(outArg(pl),iss);

    TEST_ASSERT( plGold == pl ); // BROKEN at the moment
    TEST_EQUALITY_CONST( pl.get<std::string>("Hello"), "World" );
    TEST_EQUALITY_CONST( pl.get<double>("alpha"), 1.2345 );
    TEST_EQUALITY_CONST( pl.get<bool>("isTrue"), false );
    TEST_EQUALITY_CONST( pl.get<int>("NumTimes"), 35 );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_StateSerializerStrategy, paramLists ) {
  std::string data;
  Teuchos::ParameterList plGold_A;
  {
    plGold_A.set("Hello","World");
    plGold_A.set("alpha", 1.2345 );
  }
  Teuchos::ParameterList plGold_B;
  {
    plGold_B.set("isTrue", false );
    plGold_B.set("NumTimes", 35 );
  }
  {
    std::ostringstream oss;
    XMLStateSerializerStrategy<double> ss;
    ss.serializeParameterList(plGold_A,oss);
    ss.serializeParameterList(plGold_B,oss);
    data = oss.str();
  }
  {
    std::istringstream iss(data);
    Teuchos::ParameterList pl_A;
    Teuchos::ParameterList pl_B;
    XMLStateSerializerStrategy<double> ss;
    ss.deSerializeParameterList(outArg(pl_A),iss);
    ss.deSerializeParameterList(outArg(pl_B),iss);

    TEST_ASSERT( plGold_A == pl_A ); // BROKEN at the moment
    TEST_ASSERT( plGold_B == pl_B ); // BROKEN at the moment
    TEST_EQUALITY_CONST( pl_A.get<std::string>("Hello"), "World" );
    TEST_EQUALITY_CONST( pl_A.get<double>("alpha"), 1.2345 );
    TEST_EQUALITY_CONST( pl_B.get<bool>("isTrue"), false );
    TEST_EQUALITY_CONST( pl_B.get<int>("NumTimes"), 35 );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_StateSerializerStrategy, vector ) {
  std::string data;
  int N = 15;
  {
    RCP<VectorBase<double> > vectorGold = createDefaultVector<double>(N,0.0);
    Thyra::DetachedVectorView<double> vector_view( *vectorGold );
    for (int i=0 ; i<N ; ++i ) {
      vector_view[i] = i*1.0;
    }
    std::ostringstream oss;
    XMLStateSerializerStrategy<double> ss;
    ss.serializeVectorBase(*vectorGold,oss);
    data = oss.str();
  }
  {
    std::istringstream iss(data);
    RCP<VectorBase<double> > vector = createDefaultVector<double>(N,0.0);
    XMLStateSerializerStrategy<double> ss;
    ss.deSerializeVectorBase(outArg(*vector),iss);
    {
      Thyra::ConstDetachedVectorView<double> vector_view( *vector );
      for (int i=0 ; i<N ; ++i) {
        TEST_EQUALITY_CONST( vector_view[i], i*1.0 );
      }
    }
  }
}

TEUCHOS_UNIT_TEST( Rythmos_StateSerializerStrategy, vectors ) {
  std::string data;
  int N = 15;
  {
    RCP<VectorBase<double> > vectorGold_A = createDefaultVector<double>(N,0.0);
    RCP<VectorBase<double> > vectorGold_B = createDefaultVector<double>(N,0.0);
    Thyra::DetachedVectorView<double> vector_A_view( *vectorGold_A );
    Thyra::DetachedVectorView<double> vector_B_view( *vectorGold_B );
    for (int i=0 ; i<N ; ++i ) {
      vector_A_view[i] = i*1.0;
      vector_B_view[i] = i*2.0;
    }
    std::ostringstream oss;
    XMLStateSerializerStrategy<double> ss;
    ss.serializeVectorBase(*vectorGold_A,oss);
    ss.serializeVectorBase(*vectorGold_B,oss);
    data = oss.str();
  }
  {
    std::istringstream iss(data);
    RCP<VectorBase<double> > vector_A = createDefaultVector<double>(N,0.0);
    RCP<VectorBase<double> > vector_B = createDefaultVector<double>(N,0.0);
    XMLStateSerializerStrategy<double> ss;
    ss.deSerializeVectorBase(outArg(*vector_A),iss);
    ss.deSerializeVectorBase(outArg(*vector_B),iss);
    {
      Thyra::ConstDetachedVectorView<double> vector_A_view( *vector_A );
      Thyra::ConstDetachedVectorView<double> vector_B_view( *vector_B );
      for (int i=0 ; i<N ; ++i) {
        TEST_EQUALITY_CONST( vector_A_view[i], i*1.0 );
        TEST_EQUALITY_CONST( vector_B_view[i], i*2.0 );
      }
    }
  }
}

TEUCHOS_UNIT_TEST( Rythmos_StateSerializerStrategy, vectorsAndBools ) {
  std::string data;
  int N = 15;
  {
    RCP<VectorBase<double> > vectorGold_A = createDefaultVector<double>(N,0.0);
    RCP<VectorBase<double> > vectorGold_B = createDefaultVector<double>(N,0.0);
    Thyra::DetachedVectorView<double> vector_A_view( *vectorGold_A );
    Thyra::DetachedVectorView<double> vector_B_view( *vectorGold_B );
    for (int i=0 ; i<N ; ++i ) {
      vector_A_view[i] = i*1.0;
      vector_B_view[i] = i*2.0;
    }
    std::ostringstream oss;
    XMLStateSerializerStrategy<double> ss;
    bool flag1 = true;
    bool flag2 = false;
    bool flag3 = true;
    bool flag4 = false;
    ss.serializeBool(flag1,oss);
    ss.serializeVectorBase(*vectorGold_A,oss);
    ss.serializeBool(flag2,oss);
    ss.serializeBool(flag3,oss);
    ss.serializeVectorBase(*vectorGold_B,oss);
    ss.serializeBool(flag4,oss);
    data = oss.str();
  }
  {
    std::istringstream iss(data);
    bool flag = false;
    RCP<VectorBase<double> > vector_A = createDefaultVector<double>(N,0.0);
    RCP<VectorBase<double> > vector_B = createDefaultVector<double>(N,0.0);
    XMLStateSerializerStrategy<double> ss;
    ss.deSerializeBool(outArg(flag),iss);
    TEST_ASSERT( flag == true );
    ss.deSerializeVectorBase(outArg(*vector_A),iss);
    ss.deSerializeBool(outArg(flag),iss);
    TEST_ASSERT( flag == false );
    ss.deSerializeBool(outArg(flag),iss);
    TEST_ASSERT( flag == true );
    ss.deSerializeVectorBase(outArg(*vector_B),iss);
    {
      Thyra::ConstDetachedVectorView<double> vector_A_view( *vector_A );
      Thyra::ConstDetachedVectorView<double> vector_B_view( *vector_B );
      for (int i=0 ; i<N ; ++i) {
        TEST_EQUALITY_CONST( vector_A_view[i], i*1.0 );
        TEST_EQUALITY_CONST( vector_B_view[i], i*2.0 );
      }
    }
    ss.deSerializeBool(outArg(flag),iss);
    TEST_ASSERT( flag == false );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_StateSerializerStrategy, all ) {
  std::string data;
  Teuchos::ParameterList plGold_A;
  {
    plGold_A.set("Hello","World");
    plGold_A.set("alpha", 1.2345 );
  }
  Teuchos::ParameterList plGold_B;
  {
    plGold_B.set("isTrue", false );
    plGold_B.set("NumTimes", 35 );
  }
  {
    std::ostringstream oss;
    int i = -1;
    int j = 386;
    double d = 31.4159e-1;
    double e = 1.23456;
    bool flag1 = false;
    bool flag2 = true;
    RCP<VectorBase<double> > vec1 = createDefaultVector<double>(12,0.0);
    RCP<VectorBase<double> > vec2 = createDefaultVector<double>(9,0.0);
    RCP<VectorBase<double> > vec3 = createDefaultVector<double>(9,0.0);
    {
      Thyra::DetachedVectorView<double> vec1_view( *vec1 );
      Thyra::DetachedVectorView<double> vec2_view( *vec2 );
      Thyra::DetachedVectorView<double> vec3_view( *vec3 );
      for (int k=0 ; k<12 ; ++k) {
        vec1_view[k] = k*1.0;
      }
      for (int k=0 ; k<9 ; ++k) {
        vec2_view[k] = k*2.0;
        vec3_view[k] = k*3.0;
      }
    }

    XMLStateSerializerStrategy<double> ss;
    ss.serializeInt(i,oss);
    ss.serializeInt(j,oss);
    ss.serializeParameterList(plGold_A,oss);
    ss.serializeScalar(d,oss);
    ss.serializeScalar(e,oss);
    ss.serializeVectorBase(*vec1,oss);
    ss.serializeVectorBase(*vec2,oss);
    ss.serializeParameterList(plGold_B,oss);
    ss.serializeVectorBase(*vec3,oss);
    ss.serializeBool(flag1,oss);
    ss.serializeBool(flag2,oss);
    data = oss.str();
  }
  out << "data = >>" << data << "<<" << std::endl;
  {
    std::istringstream iss(data);
    int i = 0;
    double d = 0.0;
    bool flag = true;
    XMLStateSerializerStrategy<double> ss;
    ss.deSerializeInt(outArg(i),iss);
    TEST_EQUALITY_CONST(i,-1);
    ss.deSerializeInt(outArg(i),iss);
    TEST_EQUALITY_CONST(i,386);
    {
      ParameterList pl;
      ss.deSerializeParameterList(outArg(pl),iss);
      TEST_ASSERT( pl == plGold_A ); // BROKEN
      TEST_EQUALITY_CONST( pl.get<std::string>("Hello"), "World" );
      TEST_EQUALITY_CONST( pl.get<double>("alpha"), 1.2345 );
    }
    ss.deSerializeScalar(outArg(d),iss);
    TEST_EQUALITY_CONST(d,3.14159);
    ss.deSerializeScalar(outArg(d),iss);
    TEST_EQUALITY_CONST(d,1.23456);
    {
      RCP<VectorBase<double> > vec = createDefaultVector<double>(12,0.0);
      ss.deSerializeVectorBase(outArg(*vec),iss);
      {
        Thyra::DetachedVectorView<double> vec_view( *vec );
        for (int k=0 ; k<12 ; ++k) {
          TEST_EQUALITY_CONST( vec_view[k], k*1.0 );
        }
      }
    }
    {
      RCP<VectorBase<double> > vec = createDefaultVector<double>(9,0.0);
      ss.deSerializeVectorBase(outArg(*vec),iss);
      {
        Thyra::DetachedVectorView<double> vec_view( *vec );
        for (int k=0 ; k<9 ; ++k) {
          TEST_EQUALITY_CONST( vec_view[k], k*2.0 );
        }
      }
    }
    {
      ParameterList pl;
      ss.deSerializeParameterList(outArg(pl),iss);
      TEST_ASSERT( pl == plGold_B ); // BROKEN
      TEST_EQUALITY_CONST( pl.get<bool>("isTrue"), false );
      TEST_EQUALITY_CONST( pl.get<int>("NumTimes"), 35 );
    }
    {
      RCP<VectorBase<double> > vec = createDefaultVector<double>(9,0.0);
      ss.deSerializeVectorBase(outArg(*vec),iss);
      {
        Thyra::DetachedVectorView<double> vec_view( *vec );
        for (int k=0 ; k<9 ; ++k) {
          TEST_EQUALITY_CONST( vec_view[k], k*3.0 );
        }
      }
    }
    ss.deSerializeBool(outArg(flag),iss);
    TEST_EQUALITY_CONST(flag,false);
    ss.deSerializeBool(outArg(flag),iss);
    TEST_EQUALITY_CONST(flag,true);
  }
}

} // namespace Rythmos

