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
#ifdef HAVE_KOKKOSCLASSIC_TBB
#include "Kokkos_TBBNode.hpp"
#endif
#ifdef HAVE_KOKKOSCLASSIC_THREADPOOL
#include "Kokkos_TPINode.hpp"
#endif
#ifdef HAVE_KOKKOSCLASSIC_OPENMP
#include "Kokkos_OpenMPNode.hpp"
#endif

#if (defined(HAVE_KOKKOSCLASSIC_CUSPARSE) || defined(HAVE_KOKKOSCLASSIC_CUSP)) && defined(HAVE_KOKKOSCLASSIC_THRUST)
  #include "Kokkos_ThrustGPUNode.hpp"
  #ifdef HAVE_KOKKOSCLASSIC_CUSPARSE
    #include "Kokkos_CUSPARSEOps.hpp"
  #endif
  #ifdef HAVE_KOKKOSCLASSIC_CUSP
    #include "Kokkos_CuspOps.hpp"
  #endif
  #define TEST_CUDA
  #ifdef HAVE_KOKKOSCLASSIC_CUDA_FLOAT
    #define TEST_CUDA_FLOAT
  #endif
  #ifdef HAVE_KOKKOSCLASSIC_CUDA_DOUBLE
    #define TEST_CUDA_DOUBLE
  #endif
  #ifdef HAVE_KOKKOSCLASSIC_CUDA_COMPLEX_FLOAT
    #define TEST_CUDA_COMPLEX_FLOAT
  #endif
  #ifdef HAVE_KOKKOSCLASSIC_CUDA_COMPLEX_DOUBLE
    #define TEST_CUDA_COMPLEX_DOUBLE
  #endif
#endif

namespace {

  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Kokkos::MultiVector;
  using Kokkos::DefaultArithmetic;
  using Kokkos::DefaultKernels;
  using Teuchos::ArrayRCP;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::null;
  using Teuchos::tuple;
  using std::endl;

  using Kokkos::SerialNode;
  RCP<SerialNode> snode;
#ifdef HAVE_KOKKOSCLASSIC_TBB
  using Kokkos::TBBNode;
  RCP<TBBNode> tbbnode;
#endif
#ifdef HAVE_KOKKOSCLASSIC_THREADPOOL
  using Kokkos::TPINode;
  RCP<TPINode> tpinode;
#endif
#ifdef HAVE_KOKKOSCLASSIC_OPENMP
  using Kokkos::OpenMPNode;
  RCP<OpenMPNode> ompnode;
#endif
#ifdef TEST_CUDA
  using Kokkos::ThrustGPUNode;
  RCP<ThrustGPUNode> gpunode;
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
    using Teuchos::TypeNameTraits;
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "getNode() has not been "
      "specialized for the Kokkos Node type \"" << TypeNameTraits<Node>::name ()
      << "\".  Please report this bug to the Kokkos developers.");
  }

  template <>
  RCP<SerialNode> getNode<SerialNode>() {
    if (snode == null) {
      Teuchos::ParameterList pl;
      snode = rcp(new SerialNode(pl));
    }
    return snode;
  }

#ifdef HAVE_KOKKOSCLASSIC_TBB
  template <>
  RCP<TBBNode> getNode<TBBNode>() {
    if (tbbnode == null) {
      Teuchos::ParameterList pl;
      tbbnode = rcp (new TBBNode (pl));
    }
    return tbbnode;
  }
#endif

#ifdef HAVE_KOKKOSCLASSIC_OPENMP
  template <>
  RCP<OpenMPNode> getNode<OpenMPNode>() {
    if (ompnode == null) {
      Teuchos::ParameterList pl;
      pl.set<int>("Num Threads",0);
      ompnode = rcp (new OpenMPNode (pl));
    }
    return ompnode;
  }
#endif

#ifdef HAVE_KOKKOSCLASSIC_THREADPOOL
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

#ifdef TEST_CUDA
  template <>
  RCP<ThrustGPUNode> getNode<ThrustGPUNode>() {
    if (gpunode == null) {
      Teuchos::ParameterList pl;
      pl.set<int>("Num Threads",0);
      gpunode = rcp(new ThrustGPUNode(pl));
    }
    return gpunode;
  }
#endif

  //
  // UNIT TESTS
  //

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsMatrix, TransposeMultiply, Ordinal, Scalar, Node )
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
    // generate rectangular matrix:
    // [-1 1]
    // [-1 1]
    // [-1 1]
    if (N<2) return;
    RCP<GRPH> G = rcp(new GRPH (N,2,node,null) );
    RCP<MAT>  A = rcp(new MAT  (G,null) );
    // allocate buffers for ptrs, indices and values
    const Ordinal totalNNZ = 2*N;
    ArrayRCP<size_t>  ptrs(N+1);
    ArrayRCP<Ordinal> inds(totalNNZ);
    ArrayRCP<Scalar>  vals(totalNNZ);
    // fill the buffers on the host
    {
      ptrs[0] = 0;
      for (int i=0; i != N; ++i) {
        ptrs[i] = 2*i;
        inds[2*i  ] = 0;
        inds[2*i+1] = 1;
        vals[2*i  ] = -1;
        vals[2*i+1] =  1;
      }
      ptrs[N] = totalNNZ;
    }
    G->setStructure(ptrs, inds);
    ptrs = Teuchos::null;
    inds = Teuchos::null;
    A->setValues(vals);
    vals = Teuchos::null;
    RCP<ParameterList> plist = parameterList();
    plist->set("Prepare Transpose Multiply",true);
    OPS::finalizeGraphAndMatrix(Teuchos::UNDEF_TRI,Teuchos::NON_UNIT_DIAG,*G,*A,plist);
    OPS dsm(node);
    out << "Testing with sparse ops: " << Teuchos::typeName(dsm) << std::endl;
    dsm.setGraphAndMatrix(G,A);

    ArrayRCP<Scalar> xdat, axdat, atxdat;
    xdat  = node->template allocBuffer<Scalar>(2);
    node->template copyToBuffer<Scalar>(2, tuple<Scalar>(1.0,2.0), xdat);
    axdat = node->template allocBuffer<Scalar>(N);
    atxdat = node->template allocBuffer<Scalar>(2);
    MV X(node), AX(node), ATX(node);
    X.initializeValues(  2,1,  xdat,2);
    AX.initializeValues( N,1, axdat,N);
    ATX.initializeValues(2,1,atxdat,2);
#ifdef HAVE_KOKKOSCLASSIC_DEBUG
    node->sync();
#endif
    // AX = A*X
    dsm.multiply(Teuchos::NO_TRANS,ONE,X,AX);     //  AX = ones()
    dsm.multiply(Teuchos::TRANS,   ONE,AX,ATX);   // ATX = [-N;N]
    dsm.multiply(Teuchos::TRANS,  -ONE,AX,ONE,X); // X = [1;2] - [-N;N] = [1+N;2-N]
#ifdef HAVE_KOKKOSCLASSIC_DEBUG
    node->sync();
#endif
    // AX should be all ones
    {
      ArrayRCP<const Scalar> axview = node->template viewBuffer<Scalar>(N,axdat);
      Scalar err = ZERO;
      for (int i=0; i<N; ++i) {
        err += ST::magnitude(ONE - axview[i]);
      }
      TEST_EQUALITY_CONST(err, ZERO);
    }
    // ATX should be [-N;N]
    {
      ArrayRCP<const Scalar> atxview = node->template viewBuffer<Scalar>(2,atxdat);
      TEST_COMPARE_FLOATING_ARRAYS( atxview, tuple<Scalar>(-N,N), ZERO );
    }
    // X should be [1+N;2-N]
    {
      ArrayRCP<const Scalar> xview = node->template viewBuffer<Scalar>(2,xdat);
      TEST_COMPARE_FLOATING_ARRAYS( xview, tuple<Scalar>(1+N,2-N), ZERO );
    }
  }


  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsMatrix, SparseMultiply, Ordinal, Scalar, Node )
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
    // [ 1 -1                   ]
    // [-1  3  -1               ]
    // [   -1   3  -1           ]
    // [                        ]
    // [                -1  3 -1]
    // [                   -1  2]
    if (N<2) return;
    RCP<GRPH> G = rcp(new GRPH (N,N,node,null) );
    RCP<MAT>  A = rcp(new MAT  (G,null) );
    // allocate buffers for ptrs, indices and values
    const Ordinal totalNNZ = 3*N - 2;
    ArrayRCP<size_t>  ptrs(N+1);
    ArrayRCP<Ordinal> inds(totalNNZ);
    ArrayRCP<Scalar>  vals(totalNNZ);
    // fill the buffers on the host
    {
      Ordinal NNZsofar = 0;
      ptrs[0] = NNZsofar;
      inds[NNZsofar] = 0; inds[NNZsofar+1] =  1;
      vals[NNZsofar] = 2; vals[NNZsofar+1] = -1;
      NNZsofar += 2;
      for (int i=1; i != N-1; ++i) {
        ptrs[i] = NNZsofar;
        inds[NNZsofar] = i-1; inds[NNZsofar+1] = i; inds[NNZsofar+2] = i+1;
        vals[NNZsofar] =  -1; vals[NNZsofar+1] = 3; vals[NNZsofar+2] =  -1;
        NNZsofar += 3;
      }
      ptrs[N-1] = NNZsofar;
      inds[NNZsofar] = N-2; inds[NNZsofar+1] = N-1;
      vals[NNZsofar] =  -1; vals[NNZsofar+1] = 2;
      NNZsofar += 2;
      ptrs[N]   = NNZsofar;
      TEUCHOS_TEST_FOR_EXCEPT(NNZsofar != totalNNZ);
    }
    G->setStructure(ptrs, inds);
    ptrs = Teuchos::null;
    inds = Teuchos::null;
    A->setValues(vals);
    vals = Teuchos::null;
    OPS::finalizeGraphAndMatrix(Teuchos::UNDEF_TRI,Teuchos::NON_UNIT_DIAG,*G,*A,null);
    Teuchos::EDiag diag;
    Teuchos::EUplo uplo;
    G->getMatDesc(uplo,diag);
    TEST_EQUALITY_CONST( uplo, Teuchos::UNDEF_TRI );
    TEST_EQUALITY_CONST( diag, Teuchos::NON_UNIT_DIAG );
    OPS dsm(node);
    out << "Testing with sparse ops: " << Teuchos::typeName(dsm) << std::endl;
    dsm.setGraphAndMatrix(G,A);

    ArrayRCP<Scalar> xdat, axdat;
    xdat  = node->template allocBuffer<Scalar>(N);
    axdat = node->template allocBuffer<Scalar>(N);
    MV X(node), AX(node);
    X.initializeValues( N,1, xdat,N);
    AX.initializeValues(N,1,axdat,N);
    DefaultArithmetic<MV>::Init( X,1);
#ifdef HAVE_KOKKOSCLASSIC_DEBUG
    node->sync();
#endif
    // AX = A*X
    dsm.multiply(Teuchos::NO_TRANS,ONE,X,AX);
#ifdef HAVE_KOKKOSCLASSIC_DEBUG
    node->sync();
#endif
    // AX should be all ones
    {
      ArrayRCP<const Scalar> axview = node->template viewBuffer<Scalar>(N,axdat);
      Scalar err = ZERO;
      for (int i=0; i<N; ++i) {
        err += ST::magnitude(ONE - axview[i]);
      }
      TEST_EQUALITY_CONST(err, ZERO);
    }
    // do that same multiplication, testing alpha=-1 and beta=1 accumulation
    // AX = AX - A*X = A*X - A*X = 0
    dsm.multiply(Teuchos::NO_TRANS,-ONE,X,ONE,AX);
    // AX should be zero
    {
      ArrayRCP<const Scalar> axview = node->template viewBuffer<Scalar>(N,axdat);
      ArrayRCP<const Scalar> xview = node->template viewBuffer<Scalar>(N,xdat);
      Scalar err = ZERO, nrm = ZERO;
      for (int i=0; i<N; ++i) {
        err += ST::magnitude(axview[i]);
        nrm += ST::magnitude(xview[i]);
      }
      TEST_INEQUALITY_CONST(nrm, ZERO);
      TEST_EQUALITY_CONST(err, ZERO);
    }
    xdat = null;
    axdat = null;
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsMatrix, ZeroMultiply, Ordinal, Scalar, Node )
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
    // generate zero matrix diag(zeros()), one "entry" per row
    if (N<2) return;
    RCP<GRPH> G = rcp(new GRPH (N,N,node,null) );
    RCP<MAT>  A = rcp(new MAT  (G,null) );
    // allocate buffers for ptrs, indices and values
    const Ordinal totalNNZ = N;
    ArrayRCP<size_t>  ptrs(N+1);
    ArrayRCP<Ordinal> inds(totalNNZ);
    ArrayRCP<Scalar>  vals(totalNNZ);
    // fill the buffers on the host
    {
      std::fill(vals.begin(), vals.end(), ZERO);
      for (int i=0; i < N; ++i) {
        inds[i] = i;
        ptrs[i] = i;
      }
      ptrs[N] = N;
    }
    G->setStructure(ptrs, inds);
    A->setValues(vals);
    ptrs = Teuchos::null;
    inds = Teuchos::null;
    vals = Teuchos::null;
    OPS::finalizeGraphAndMatrix(Teuchos::UNDEF_TRI,Teuchos::NON_UNIT_DIAG,*G,*A,null);
    OPS dsm(node);
    out << "Testing with sparse ops: " << Teuchos::typeName(dsm) << std::endl;
    dsm.setGraphAndMatrix(G,A);

    ArrayRCP<Scalar> xdat, axdat;
    xdat  = node->template allocBuffer<Scalar>(N);
    axdat = node->template allocBuffer<Scalar>(N);
    MV X(node), AX(node);
    X.initializeValues( N,1, xdat,N);
    AX.initializeValues(N,1,axdat,N);
    DefaultArithmetic<MV>::Init( X,1);
    DefaultArithmetic<MV>::Init(AX,1);
    // AX = A*X
    dsm.multiply(Teuchos::NO_TRANS,ONE,X,AX);
    // AX should be zero
    {
      ArrayRCP<const Scalar> axview = node->template viewBuffer<Scalar>(N,axdat);
      Scalar err = ZERO;
      for (int i=0; i<N; ++i) {
        err += ST::magnitude(axview[i]);
      }
      TEST_EQUALITY_CONST(err, ZERO);
    }
    // do that same multiplication, testing alpha=1 and beta=0 "accumulation"
    // AX = 0*AX + A*X = 0*ones() + 1*Zero*ones()
    DefaultArithmetic<MV>::Init( X,1);
    DefaultArithmetic<MV>::Init(AX,1);
    dsm.multiply(Teuchos::NO_TRANS,ONE,X,ZERO,AX);
    // AX should again be zero
    {
      ArrayRCP<const Scalar> axview = node->template viewBuffer<Scalar>(N,axdat);
      Scalar err = ZERO;
      for (int i=0; i<N; ++i) {
        err += ST::magnitude(axview[i]);
      }
      TEST_EQUALITY_CONST(err, ZERO);
    }
    xdat = null;
    axdat = null;
  }

#include "CrsMatrix_DefaultMultiplyTests.hpp"

#define ALL_UNIT_TESTS_ORDINAL_SCALAR_NODE( ORDINAL, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsMatrix,        SparseMultiply, ORDINAL, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsMatrix,        ZeroMultiply,   ORDINAL, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsMatrix,     TransposeMultiply, ORDINAL, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( DefaultSparseOps, NodeTest,       ORDINAL, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( DefaultSparseOps, ResubmitMatrix, ORDINAL, SCALAR, NODE )

#define UNIT_TEST_SERIALNODE(ORDINAL, SCALAR) \
      ALL_UNIT_TESTS_ORDINAL_SCALAR_NODE( ORDINAL, SCALAR, SerialNode )

#ifdef HAVE_KOKKOSCLASSIC_TBB
#define UNIT_TEST_TBBNODE(ORDINAL, SCALAR) \
      ALL_UNIT_TESTS_ORDINAL_SCALAR_NODE( ORDINAL, SCALAR, TBBNode )
#else
#define UNIT_TEST_TBBNODE(ORDINAL, SCALAR)
#endif

#ifdef HAVE_KOKKOSCLASSIC_OPENMP
#define UNIT_TEST_OPENMPNODE(ORDINAL, SCALAR) \
      ALL_UNIT_TESTS_ORDINAL_SCALAR_NODE( ORDINAL, SCALAR, OpenMPNode )
#else
#define UNIT_TEST_OPENMPNODE(ORDINAL, SCALAR)
#endif

#ifdef HAVE_KOKKOSCLASSIC_THREADPOOL
#define UNIT_TEST_TPINODE(ORDINAL, SCALAR) \
      ALL_UNIT_TESTS_ORDINAL_SCALAR_NODE( ORDINAL, SCALAR, TPINode )
#else
#define UNIT_TEST_TPINODE(ORDINAL, SCALAR)
#endif

#define UNIT_TEST_GROUP_ORDINAL_SCALAR( ORDINAL, SCALAR ) \
        UNIT_TEST_SERIALNODE( ORDINAL, SCALAR ) \
        UNIT_TEST_TBBNODE( ORDINAL, SCALAR ) \
        UNIT_TEST_OPENMPNODE( ORDINAL, SCALAR ) \
        UNIT_TEST_TPINODE( ORDINAL, SCALAR )

#define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
        UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, int) \
        UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, float)

typedef std::complex<float>  ComplexFloat;
typedef std::complex<double> ComplexDouble;
#ifdef TEST_CUDA_FLOAT
ALL_UNIT_TESTS_ORDINAL_SCALAR_NODE( int, float, ThrustGPUNode )
#endif
#ifdef TEST_CUDA_DOUBLE
ALL_UNIT_TESTS_ORDINAL_SCALAR_NODE( int, double, ThrustGPUNode )
#endif
#ifdef TEST_CUDA_COMPLEX_FLOAT
ALL_UNIT_TESTS_ORDINAL_SCALAR_NODE( int, ComplexFloat, ThrustGPUNode )
#endif
#ifdef TEST_CUDA_COMPLEX_DOUBLE
ALL_UNIT_TESTS_ORDINAL_SCALAR_NODE( int, ComplexDouble, ThrustGPUNode )
#endif

#ifdef HAVE_KOKKOSCLASSIC_CUSP 
ALL_UNIT_TESTS_ORDINAL_SCALAR_NODE( short, float, ThrustGPUNode )
#endif

UNIT_TEST_GROUP_ORDINAL(int)
typedef short int ShortInt; UNIT_TEST_GROUP_ORDINAL(ShortInt)

}

