// ************************************************************************
// 
//          Kokkos: A Fast Kernel Package
//              Copyright (2004) Sandia Corporation
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
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>
#include <Teuchos_TypeNameTraits.hpp>

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_DefaultNode.hpp"
#include "Kokkos_MultiVector.hpp"
#include "Kokkos_Vector.hpp"
#include "Kokkos_Version.hpp"

namespace {

  using Kokkos::DefaultNode;
  using Kokkos::Vector;
  using Kokkos::MultiVector;

  int N = 1000;

  typedef Kokkos::DefaultNode::DefaultNodeType Node;

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
    clp.setOption("test-size",&N,"Vector length for tests.");
  }

  //
  // UNIT TESTS
  // 

  // check that default constructor zeros out, for both V and MV
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, DefaultConstructor, Scalar, Ordinal )
  {
    MultiVector<Scalar,Ordinal,Node> A;
    TEST_EQUALITY_CONST(A.getNumRows(), 0);
    TEST_EQUALITY_CONST(A.getNumCols(), 0);
    TEST_EQUALITY_CONST(A.getStride(), 0);
  }

  // check copy constructor
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, CopyConstructor, Scalar, Ordinal )
  {
    MultiVector<Scalar,Ordinal,Node> A;
    typename Node::template buffer<Scalar>::buffer_t 
      buf = A.getNode().template allocBuffer<Scalar>(2*N);
    A.initializeValues(N,2,buf,N);
    {
      MultiVector<Scalar,Ordinal,Node> Acopy(A);
      TEST_EQUALITY_CONST(Acopy.getNumRows(), N);
      TEST_EQUALITY_CONST(Acopy.getNumCols(), 2);
      TEST_EQUALITY_CONST(Acopy.getStride(), N);
      TEST_EQUALITY(Acopy.getValues(0), buf);
      TEST_INEQUALITY(Acopy.getValues(1), buf);
    }
  }

  // check that non-default constructor honors given parameters
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, InitializeAndAccess, Scalar, Ordinal )
  {
    MultiVector<Scalar,Ordinal,Node> A;
    typename Node::template buffer<Scalar>::buffer_t 
      buf = A.getNode().template allocBuffer<Scalar>(2*N);
    A.initializeValues(N,2,buf,N);
    TEST_EQUALITY_CONST(A.getNumRows(), N);
    TEST_EQUALITY_CONST(A.getNumCols(), 2);
    TEST_EQUALITY_CONST(A.getStride(), N);
    TEST_EQUALITY(A.getValues(0), buf);
    TEST_INEQUALITY(A.getValues(1), buf);
     A.getNode().template freeBuffer<Scalar>(buf);
   }
 
   // check that default constructor zeros out, for both V and MV
   TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Vector, DefaultConstructor, Scalar, Ordinal )
   {
     Vector<Scalar,Ordinal,Node> a;
     TEST_EQUALITY_CONST(a.getLength(), 0);
     TEST_EQUALITY_CONST(a.getInc(), 0);
   }
 
   // check copy constructor
   TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Vector, CopyConstructor, Scalar, Ordinal )
   {
     Vector<Scalar,Ordinal,Node> a;
     typename Node::template buffer<Scalar>::buffer_t 
       buf = a.getNode().template allocBuffer<Scalar>(N);
     a.initializeValues(N,buf,1);
     {
       Vector<Scalar,Ordinal,Node> acopy(a);
       TEST_EQUALITY_CONST(acopy.getLength(), N);
       TEST_EQUALITY_CONST(acopy.getInc(), 1);
       TEST_EQUALITY(acopy.getValues(), buf);
     }
     a.getNode().template freeBuffer<Scalar>(buf);
   }
 
   TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Vector, InitializeAndAccess, Scalar, Ordinal )
   {
     Vector<Scalar,Ordinal,Node> a;
     typename Node::template buffer<Scalar>::buffer_t 
       buf = a.getNode().template allocBuffer<Scalar>(N);
     a.initializeValues(N,buf,1);
     TEST_EQUALITY_CONST(a.getLength(), N);
     TEST_EQUALITY_CONST(a.getInc(), 1);
     TEST_EQUALITY(a.getValues(), buf);
     a.getNode().template freeBuffer<Scalar>(buf);
   }
 
 #define UNIT_TEST_GROUP_ORDINAL_SCALAR( ORDINAL, SCALAR ) \
       TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, DefaultConstructor, SCALAR, ORDINAL ) \
       TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, CopyConstructor   , SCALAR, ORDINAL ) \
       TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, InitializeAndAccess, SCALAR, ORDINAL ) \
       TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Vector,      DefaultConstructor, SCALAR, ORDINAL ) \
       TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Vector,      CopyConstructor   , SCALAR, ORDINAL ) \
       TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Vector,      InitializeAndAccess, SCALAR, ORDINAL )
 
 #define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
          UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, int) \
          UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, float)
      UNIT_TEST_GROUP_ORDINAL(int)
      typedef short int ShortInt; UNIT_TEST_GROUP_ORDINAL(ShortInt)
 
}
 
