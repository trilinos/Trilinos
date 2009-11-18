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

#include "Kokkos_MultiVector.hpp"
#include "Kokkos_DefaultArithmetic.hpp"
#include "Kokkos_Version.hpp"
#include "Kokkos_SerialNode.hpp"
#ifdef HAVE_KOKKOS_TBB
#include "Kokkos_TBBNode.hpp"
#endif
#ifdef HAVE_KOKKOS_THREADPOOL
#include "Kokkos_TPINode.hpp"
#endif
#ifdef HAVE_KOKKOS_THRUST
#include "Kokkos_ThrustGPUNode.hpp"
#endif

namespace {

  using Kokkos::MultiVector;
  using Kokkos::DefaultArithmetic;
  using Kokkos::SerialNode;
  using Teuchos::ScalarTraits;
  using Teuchos::ArrayRCP;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::null;

  RCP<SerialNode> snode;
#ifdef HAVE_KOKKOS_TBB
  using Kokkos::TBBNode;
  RCP<TBBNode> tbbnode;
#endif
#ifdef HAVE_KOKKOS_THREADPOOL
  using Kokkos::TPINode;
  RCP<TPINode> tpinode;
#endif
#ifdef HAVE_KOKKOS_THRUST
  using Kokkos::ThrustGPUNode;
  RCP<ThrustGPUNode> thrustnode;
#endif

  int N = 1000;

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
    clp.setOption("test-size",&N,"Vector length for tests.");
  }

  template <class Node>
  RCP<Node> getNode() {
    assert(false);
  }

  template <>
  RCP<SerialNode> getNode<SerialNode>() {
    if (snode == null) {
      Teuchos::ParameterList pl;
      snode = rcp(new SerialNode(pl));
    }
    return snode;
  }

#ifdef HAVE_KOKKOS_TBB
  template <>
  RCP<TBBNode> getNode<TBBNode>() {
    if (tbbnode == null) {
      Teuchos::ParameterList pl;
      pl.set<int>("Num Threads",0);
      tbbnode = rcp(new TBBNode(pl));
    }
    return tbbnode;
  }
#endif

#ifdef HAVE_KOKKOS_THREADPOOL
  template <>
  RCP<TPINode> getNode<TPINode>() {
    if (tpinode == null) {
      Teuchos::ParameterList pl;
      pl.set<int>("Num Threads",0);
      tpinode = rcp(new TPINode(pl));
    }
    return tpinode;
  }
#endif

#ifdef HAVE_KOKKOS_THRUST
  template <>
  RCP<ThrustGPUNode> getNode<ThrustGPUNode>() {
    if (thrustnode == null) {
      Teuchos::ParameterList pl;
      pl.set<int>("Num Threads",0);
      thrustnode = rcp(new ThrustGPUNode(pl));
    }
    return thrustnode;
  }
#endif

  //
  // UNIT TESTS
  // 

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( MultiVector, Scale, Scalar, Ordinal, Node )
  {
    RCP<Node> node = getNode<Node>();
    typedef MultiVector<Scalar,Node> MV;
    const int numVecs = 5;
    MV MV1(node), MV2(node), vec(node);
    ArrayRCP<Scalar> buf = node->template allocBuffer<Scalar>(2*numVecs*N);
    MV1.initializeValues(N,numVecs,buf          ,N);
    MV2.initializeValues(N,numVecs,buf+numVecs*N,N);
    vec.initializeValues(N,1,buf,N);                    // MV1 collocated with vec
    DefaultArithmetic<MV>::Init(MV1,2);                 // MV1 = twos()
    DefaultArithmetic<MV>::Recip( MV2,(const MV&)MV1);  // MV2 = 1./MV1 = ones()/two()
    ArrayRCP<const Scalar> viewMV2 = node->template viewBuffer<Scalar>(numVecs*N, buf+numVecs*N);
    ArrayRCP<Scalar> halfs(numVecs*N, static_cast<Scalar>(1.0/2.0) );
    TEST_COMPARE_FLOATING_ARRAYS(viewMV2, halfs(), ScalarTraits<Scalar>::eps());
    viewMV2 = Teuchos::null;
    buf = Teuchos::null;
  }


#define UNIT_TEST_SERIALNODE(SCALAR, ORDINAL) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( MultiVector, Scale, SCALAR, ORDINAL, SerialNode )

#ifdef HAVE_KOKKOS_TBB
#define UNIT_TEST_TBBNODE(SCALAR, ORDINAL) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( MultiVector, Scale, SCALAR, ORDINAL, TBBNode )
#else
#define UNIT_TEST_TBBNODE(SCALAR, ORDINAL)
#endif

#ifdef HAVE_KOKKOS_THREADPOOL
#define UNIT_TEST_TPINODE(SCALAR, ORDINAL) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( MultiVector, Scale, SCALAR, ORDINAL, TPINode )
#else
#define UNIT_TEST_TPINODE(SCALAR, ORDINAL)
#endif

#ifdef HAVE_KOKKOS_THRUST
#define UNIT_TEST_THRUSTGPUNODE(SCALAR, ORDINAL) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( MultiVector, Scale, SCALAR, ORDINAL, ThrustGPUNode )
#else
#define UNIT_TEST_THRUSTGPUNODE(SCALAR, ORDINAL)
#endif

#define UNIT_TEST_GROUP_ORDINAL_SCALAR( ORDINAL, SCALAR ) \
      UNIT_TEST_SERIALNODE( SCALAR, ORDINAL ) \
      UNIT_TEST_TBBNODE( SCALAR, ORDINAL ) \
      UNIT_TEST_TPINODE( SCALAR, ORDINAL ) \
      UNIT_TEST_THRUSTGPUNODE( SCALAR, ORDINAL )

#    define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
         UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, float)
     UNIT_TEST_GROUP_ORDINAL(int)
     typedef short int ShortInt; UNIT_TEST_GROUP_ORDINAL(ShortInt)

}
