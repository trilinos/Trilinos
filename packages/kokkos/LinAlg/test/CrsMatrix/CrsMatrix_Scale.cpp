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
#include <Teuchos_ScalarTraits.hpp>

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_DefaultArithmetic.hpp"
#include "Kokkos_DefaultKernels.hpp"
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
  using Kokkos::CrsMatrix;
  using Kokkos::CrsGraph;
  using Kokkos::DefaultArithmetic;
  using Kokkos::DefaultKernels;
  using Kokkos::SerialNode;
  using Teuchos::ArrayRCP;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::null;
  using std::endl;

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
      pl.set<int>("Verbose",1);
      thrustnode = rcp(new ThrustGPUNode(pl));
    }
    return thrustnode;
  }
#endif

  //
  // UNIT TESTS
  // 

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsMatrix, Scale, Ordinal, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    typedef typename DefaultKernels<Scalar,Ordinal,Node>::SparseOps DSM;
    typedef CrsGraph<Ordinal,Node,DSM>                             GRPH;
    typedef CrsMatrix<Scalar,Ordinal,Node,DSM>                      MAT;
    typedef MultiVector<Scalar,Node>                                 MV;
    typedef Teuchos::ScalarTraits<Scalar>                            ST;
    const Scalar ONE = ST::one(),
                ZERO = ST::zero();
    // generate tridiagonal matrix:
    // [ 2 -1                   ]
    // [-1  3  -1               ]
    // [   -1   3  -1           ]
    // [                        ]
    // [                -1  3 -1]
    // [                   -1  2]
    if (N<2) return;
    GRPH G(N,node);
    MAT  A(G);
    // allocate buffers for offsets, indices and values
    const size_t totalNNZ = 3*N - 2;
    ArrayRCP<size_t> offsets(N+1);
    ArrayRCP<Ordinal>   inds(totalNNZ);
    ArrayRCP<Scalar>    vals(totalNNZ);
    // fill the buffers on the host
    {
      size_t NNZsofar = 0;
      offsets[0] = NNZsofar;
      inds[NNZsofar] = 0; inds[NNZsofar+1] =  1;
      vals[NNZsofar] = 2; vals[NNZsofar+1] = -1;
      NNZsofar += 2;
      for (int i=1; i != N-1; ++i) {
        offsets[i] = NNZsofar;
        inds[NNZsofar] = i-1; inds[NNZsofar+1] = i; inds[NNZsofar+2] = i+1;
        vals[NNZsofar] =  -1; vals[NNZsofar+1] = 3; vals[NNZsofar+2] =  -1;
        NNZsofar += 3;
      }
      offsets[N-1] = NNZsofar;
      inds[NNZsofar] = N-2; inds[NNZsofar+1] = N-1;
      vals[NNZsofar] =  -1; vals[NNZsofar+1] = 2;
      NNZsofar += 2;
      offsets[N]   = NNZsofar;
      TEST_FOR_EXCEPT(NNZsofar != totalNNZ);
    }
    G.set1DStructure(inds, offsets, offsets.persistingView(1,N));
    offsets = Teuchos::null;
    inds    = Teuchos::null;
    A.set1DValues(vals);
    vals    = Teuchos::null;
    A.finalize(true);
    typename DSM::template rebind<Scalar>::other dsm(node);
    out << "Testing with sparse ops: " << Teuchos::typeName(dsm) << std::endl;
    dsm.initializeStructure(G);
    dsm.initializeValues(A);

    ArrayRCP<Scalar> xdat, axdat;
    xdat  = node->template allocBuffer<Scalar>(N);
    axdat = node->template allocBuffer<Scalar>(N);
    MV X(node), AX(node);
    X.initializeValues( N,1, xdat,N);
    AX.initializeValues(N,1,axdat,N);

    // Scale by x, unscale by x and test the multiply.
    DefaultArithmetic<MV>::Init( X,2);
    dsm.leftScale(X);
    DefaultArithmetic<MV>::Init( X,0.5);
    dsm.leftScale(X);

    DefaultArithmetic<MV>::Init( X,1.0);
    dsm.multiply(Teuchos::NO_TRANS,ONE,X,AX);

    // AX should be all ones
    {
      ArrayRCP<const Scalar> axview = node->template viewBuffer<Scalar>(N,axdat);
      Scalar err = ZERO;
      for (int i=0; i<N; ++i) {
        err += ST::magnitude(ONE - axview[i]);
      }
      TEST_EQUALITY_CONST(err, ZERO);
    }

    xdat = null;
    axdat = null;
  }

#include "CrsMatrix_DefaultMultiplyTests.hpp"

#define ALL_UNIT_TESTS_ORDINAL_SCALAR_NODE( ORDINAL, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsMatrix,        Scale, ORDINAL, SCALAR, NODE ) 

#define UNIT_TEST_SERIALNODE(ORDINAL, SCALAR) \
      ALL_UNIT_TESTS_ORDINAL_SCALAR_NODE( ORDINAL, SCALAR, SerialNode )

#ifdef HAVE_KOKKOS_TBB
#define UNIT_TEST_TBBNODE(ORDINAL, SCALAR) \
      ALL_UNIT_TESTS_ORDINAL_SCALAR_NODE( ORDINAL, SCALAR, TBBNode )
#else
#define UNIT_TEST_TBBNODE(ORDINAL, SCALAR)
#endif

#ifdef HAVE_KOKKOS_THREADPOOL
#define UNIT_TEST_TPINODE(ORDINAL, SCALAR) \
      ALL_UNIT_TESTS_ORDINAL_SCALAR_NODE( ORDINAL, SCALAR, TPINode )
#else
#define UNIT_TEST_TPINODE(ORDINAL, SCALAR)
#endif

#ifdef HAVE_KOKKOS_THRUST
#define UNIT_TEST_THRUSTGPUNODE(ORDINAL, SCALAR) \
      ALL_UNIT_TESTS_ORDINAL_SCALAR_NODE( ORDINAL, SCALAR, ThrustGPUNode )
#else
#define UNIT_TEST_THRUSTGPUNODE(ORDINAL, SCALAR)
#endif

#define UNIT_TEST_GROUP_ORDINAL_SCALAR( ORDINAL, SCALAR ) \
        UNIT_TEST_SERIALNODE( ORDINAL, SCALAR ) \
        UNIT_TEST_TBBNODE( ORDINAL, SCALAR ) \
        UNIT_TEST_TPINODE( ORDINAL, SCALAR ) \

#define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
        UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, float)

     UNIT_TEST_GROUP_ORDINAL(int)
     typedef short int ShortInt; UNIT_TEST_GROUP_ORDINAL(ShortInt)

}

