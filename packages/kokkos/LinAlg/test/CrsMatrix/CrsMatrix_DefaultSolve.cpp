//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
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

#define TEST_DATA_FOR_ONES(dat) {                                               \
    ArrayRCP<const Scalar> view = node->template viewBuffer<Scalar>(N,dat);     \
    Scalar err = ZERO;                                                          \
    for (int i=0; i<N; ++i) {                                                   \
      err += ST::magnitude(ONE - view[i]);                                    \
    }                                                                           \
    TEST_EQUALITY_CONST(err, ZERO);                                             \
  }                                                                             

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsMatrix, SparseSolveIdentity, Ordinal, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    typedef typename DefaultKernels<Scalar,Ordinal,Node>::SparseOps          DSM;
    typedef typename DSM::template bind_scalar<Scalar>::other_type           OPS;
    typedef typename OPS::template matrix<Scalar,Ordinal,Node>::matrix_type  MAT;
    typedef typename OPS::template graph<Ordinal,Node>::graph_type          GRPH;
    typedef MultiVector<Scalar,Node>                                          MV;
    typedef Teuchos::ScalarTraits<Scalar>                                     ST;
    const Scalar ONE = ST::one(),
                ZERO = ST::zero();
    // generate tridiagonal matrix:
    // [ 1 0           ]
    // [   1 0         ]
    // [     1 0       ]
    // [         . .   ]
    // [            1 0]
    // [              1]
    if (N<2) return;
    RCP<GRPH> G = rcp(new GRPH (N,node,null) );
    RCP<MAT>  A = rcp(new MAT  (G,null) );
    // allocate buffers for offsets, indices and values
    const size_t totalNNZ = 2*N-1;
    ArrayRCP<size_t> offsets(N+1);
    ArrayRCP<Ordinal>   inds(totalNNZ);
    ArrayRCP<Scalar>    vals(totalNNZ);
    // fill the buffers on the host
    {
      size_t num = 0;
      for (int i=0; i < N-1; ++i) {
        offsets[i] = num;
        inds[num] = i; inds[num+1] = i+1;
        vals[num] = 1; vals[num+1] = 0;
        num += 2;
      }
      inds[num] = N-1;
      vals[num] = 1;
      offsets[N] = num+1;
    }
    G->setStructure(offsets, inds);
    offsets = Teuchos::null;
    inds    = Teuchos::null;
    A->setValues(vals);
    vals    = Teuchos::null;
    OPS::finalizeGraphAndMatrix(*G,*A,null);
    OPS dsm(node);
    out << "Testing with sparse ops: " << Teuchos::typeName(dsm) << std::endl;
    dsm.setGraphAndMatrix(G,A);

    // A is the identity, which allows us to easily test transpose and non-transpose, upper and lower tri
    ArrayRCP<Scalar> ydat, x1dat, x2dat, x3dat, x4dat;
    ydat  = node->template allocBuffer<Scalar>(N);
    x1dat = node->template allocBuffer<Scalar>(N);
    x2dat = node->template allocBuffer<Scalar>(N);
    x3dat = node->template allocBuffer<Scalar>(N);
    x4dat = node->template allocBuffer<Scalar>(N);
    MV Y(node), X1(node), X2(node), X3(node), X4(node); 
    Y.initializeValues(N,1, ydat,N);
    X1.initializeValues(N,1,x1dat,N);
    X2.initializeValues(N,1,x2dat,N);
    X3.initializeValues(N,1,x3dat,N);
    X4.initializeValues(N,1,x4dat,N);
    // solve A*X=Y
    DefaultArithmetic<MV>::Init(Y,1);
    dsm.solve(Teuchos::NO_TRANS,Teuchos::UPPER_TRI,Teuchos::NON_UNIT_DIAG,Y,X1);
    TEST_DATA_FOR_ONES(x1dat)
    dsm.solve(Teuchos::NO_TRANS,Teuchos::LOWER_TRI,Teuchos::NON_UNIT_DIAG,Y,X2);
    TEST_DATA_FOR_ONES(x2dat)
    dsm.solve(Teuchos::CONJ_TRANS,Teuchos::UPPER_TRI,Teuchos::NON_UNIT_DIAG,Y,X3);
    TEST_DATA_FOR_ONES(x3dat)
    dsm.solve(Teuchos::CONJ_TRANS,Teuchos::LOWER_TRI,Teuchos::NON_UNIT_DIAG,Y,X4);
    TEST_DATA_FOR_ONES(x4dat)
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsMatrix, SparseSolveImplicitIdentity, Ordinal, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    typedef typename DefaultKernels<Scalar,Ordinal,Node>::SparseOps          DSM;
    typedef typename DSM::template bind_scalar<Scalar>::other_type           OPS;
    typedef typename OPS::template matrix<Scalar,Ordinal,Node>::matrix_type  MAT;
    typedef typename OPS::template graph<Ordinal,Node>::graph_type          GRPH;
    typedef MultiVector<Scalar,Node>                                          MV;
    typedef Teuchos::ScalarTraits<Scalar>                                     ST;
    const Scalar ONE = ST::one(),
                ZERO = ST::zero();
    // generate tridiagonal matrix:
    // [ 1 0           ]
    // [   1 0         ]
    // [     1 0       ]
    // [         . .   ]
    // [            1 0]
    // [              1]
    // but don't store the diagonal
    if (N<2) return;
    RCP<GRPH> G = rcp(new GRPH (N,node,null) );
    RCP<MAT>  A = rcp(new MAT  (G,null) );
    // allocate buffers for offsets, indices and values
    const size_t totalNNZ = N-1;
    ArrayRCP<size_t> offsets(N+1);
    ArrayRCP<Ordinal>   inds(totalNNZ);
    ArrayRCP<Scalar>    vals(totalNNZ);
    // fill the buffers on the host
    {
      size_t num = 0;
      for (int i=0; i < N-1; ++i) {
        offsets[i] = i;
        vals[i] = 0;
        inds[i] = i;
      }
      offsets[N] = N-1;
    }
    G->setStructure(offsets, inds);
    offsets = Teuchos::null;
    inds    = Teuchos::null;
    A->setValues(vals);
    vals    = Teuchos::null;
    OPS::finalizeGraphAndMatrix(*G,*A,null);
    OPS dsm(node);
    out << "Testing with sparse ops: " << Teuchos::typeName(dsm) << std::endl;
    dsm.setGraphAndMatrix(G,A);

    // A is the identity, which allows us to easily test transpose and non-transpose, upper and lower tri
    ArrayRCP<Scalar> ydat, x1dat, x2dat, x3dat, x4dat;
    ydat  = node->template allocBuffer<Scalar>(N);
    x1dat = node->template allocBuffer<Scalar>(N);
    x2dat = node->template allocBuffer<Scalar>(N);
    x3dat = node->template allocBuffer<Scalar>(N);
    x4dat = node->template allocBuffer<Scalar>(N);
    MV Y(node), X1(node), X2(node), X3(node), X4(node); 
    Y.initializeValues(N,1, ydat,N);
    X1.initializeValues(N,1,x1dat,N);
    X2.initializeValues(N,1,x2dat,N);
    X3.initializeValues(N,1,x3dat,N);
    X4.initializeValues(N,1,x4dat,N);
    // solve A*X=Y
    DefaultArithmetic<MV>::Init(Y,1);
    dsm.solve(Teuchos::NO_TRANS,Teuchos::UPPER_TRI,Teuchos::UNIT_DIAG,Y,X1);
    TEST_DATA_FOR_ONES(x1dat)
    dsm.solve(Teuchos::NO_TRANS,Teuchos::LOWER_TRI,Teuchos::UNIT_DIAG,Y,X2);
    TEST_DATA_FOR_ONES(x2dat)
    dsm.solve(Teuchos::CONJ_TRANS,Teuchos::UPPER_TRI,Teuchos::UNIT_DIAG,Y,X3);
    TEST_DATA_FOR_ONES(x3dat)
    dsm.solve(Teuchos::CONJ_TRANS,Teuchos::LOWER_TRI,Teuchos::UNIT_DIAG,Y,X4);
    TEST_DATA_FOR_ONES(x4dat)
  }

#define ALL_UNIT_TESTS_ORDINAL_SCALAR_NODE( ORDINAL, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsMatrix,        SparseSolveIdentity, ORDINAL, SCALAR, NODE )

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
        UNIT_TEST_THRUSTGPUNODE( ORDINAL, SCALAR )

#define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
        UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, int) \
        UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, float)

     UNIT_TEST_GROUP_ORDINAL(int)
     typedef short int ShortInt; UNIT_TEST_GROUP_ORDINAL(ShortInt)

}

