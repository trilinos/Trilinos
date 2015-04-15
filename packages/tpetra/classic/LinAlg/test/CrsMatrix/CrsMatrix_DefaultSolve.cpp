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

#include "Kokkos_DefaultNode.hpp"

namespace {

  using KokkosClassic::MultiVector;
  using KokkosClassic::DefaultArithmetic;
  using KokkosClassic::DefaultKernels;
  using Teuchos::ArrayRCP;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::null;
  using std::endl;

  int N = 1000;

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
    clp.setOption("test-size",&N,"Vector length for tests.");
  }

  template <class Node>
  RCP<Node> getNode() {
    Teuchos::ParameterList pl;
    return rcp (new Node (pl));
  }

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

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsMatrix, SparseSolveIdentityLower, Ordinal, Scalar, Node )
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
    // generate lower triangular identity matrix:
    // [ 1           ]
    // [ 0 1         ]
    // [   0 1       ]
    // [     . .     ]
    // [        0 1  ]
    // [          0 1]
    if (N<2) return;
    RCP<GRPH> G = rcp(new GRPH (N,N,node,null) );
    RCP<MAT>  A = rcp(new MAT  (G,null) );
    // allocate buffers for ptrs, indices and values
    const Ordinal totalNNZ = 2*N-1;
    ArrayRCP<size_t>   ptrs(N+1);
    ArrayRCP<Ordinal>  inds(totalNNZ);
    ArrayRCP<Scalar>   vals(totalNNZ);
    // fill the buffers on the host
    {
      Ordinal num = 0;
      ptrs[0] = num;
      inds[num] = 0;
      vals[num] = 1;
      num += 1;
      for (int i=1; i < N; ++i) {
        ptrs[i] = num;
        inds[num] = i-1; inds[num+1] = i;
        vals[num] = 0; vals[num+1] = 1;
        num += 2;
      }
      ptrs[N] = num;
    }
    G->setStructure(ptrs, inds);
    ptrs = Teuchos::null;
    inds = Teuchos::null;
    A->setValues(vals);
    vals = Teuchos::null;
    OPS::finalizeGraphAndMatrix(Teuchos::LOWER_TRI,Teuchos::NON_UNIT_DIAG,*G,*A,null);
    Teuchos::EDiag diag;
    Teuchos::EUplo uplo;
    G->getMatDesc(uplo,diag);
    TEST_EQUALITY_CONST( uplo, Teuchos::LOWER_TRI );
    TEST_EQUALITY_CONST( diag, Teuchos::NON_UNIT_DIAG );
    OPS dsm(node);
    out << "Testing with sparse ops: " << Teuchos::typeName(dsm) << std::endl;
    dsm.setGraphAndMatrix(G,A);

    // A is the identity, which allows us to easily test transpose and non-transpose, upper and lower tri
    ArrayRCP<Scalar> ydat, x1dat, x3dat;
    ydat  = node->template allocBuffer<Scalar>(N);
    x1dat = node->template allocBuffer<Scalar>(N);
    x3dat = node->template allocBuffer<Scalar>(N);
    MV Y(node), X1(node), X3(node);
    Y.initializeValues(N,1, ydat,N);
    X1.initializeValues(N,1,x1dat,N);
    X3.initializeValues(N,1,x3dat,N);
    // solve A*X=Y
    DefaultArithmetic<MV>::Init(Y,1);
    dsm.solve(Teuchos::NO_TRANS,Y,X1);
    TEST_DATA_FOR_ONES(x1dat)
    dsm.solve(Teuchos::CONJ_TRANS,Y,X3);
    TEST_DATA_FOR_ONES(x3dat)
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsMatrix, SparseSolveImplicitIdentityLower, Ordinal, Scalar, Node )
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
    // generate lower triangular identity matrix:
    // [ 1           ]
    // [ 0 1         ]
    // [   0 1       ]
    // [     . .     ]
    // [        0 1  ]
    // [          0 1]
    // but don't store the diagonal
    if (N<2) return;
    RCP<GRPH> G = rcp(new GRPH (N,N,node,null) );
    RCP<MAT>  A = rcp(new MAT  (G,null) );
    // allocate buffers for ptrs, indices and values
    const Ordinal totalNNZ = N-1;
    ArrayRCP<size_t>   ptrs(N+1);
    ArrayRCP<Ordinal>  inds(totalNNZ);
    ArrayRCP<Scalar>   vals(totalNNZ);
    // fill the buffers on the host
    {
      Ordinal num = 0;
      ptrs[0] = num;
      for (int i=1; i < N; ++i) {
        ptrs[i] = num;
        inds[num] = i-1;
        vals[num] = 0;
        num += 1;
      }
      ptrs[N] = num;
    }
    G->setStructure(ptrs, inds);
    ptrs = Teuchos::null;
    inds = Teuchos::null;
    A->setValues(vals);
    vals = Teuchos::null;
    OPS::finalizeGraphAndMatrix(Teuchos::LOWER_TRI,Teuchos::UNIT_DIAG,*G,*A,null);
    Teuchos::EDiag diag;
    Teuchos::EUplo uplo;
    G->getMatDesc(uplo,diag);
    TEST_EQUALITY_CONST( uplo, Teuchos::LOWER_TRI );
    TEST_EQUALITY_CONST( diag, Teuchos::UNIT_DIAG );
    OPS dsm(node);
    out << "Testing with sparse ops: " << Teuchos::typeName(dsm) << std::endl;
    dsm.setGraphAndMatrix(G,A);

    // A is the identity, which allows us to easily test transpose and non-transpose, upper and lower tri
    ArrayRCP<Scalar> ydat, x1dat, x3dat;
    ydat  = node->template allocBuffer<Scalar>(N);
    x1dat = node->template allocBuffer<Scalar>(N);
    x3dat = node->template allocBuffer<Scalar>(N);
    MV Y(node), X1(node), X3(node);
    Y.initializeValues(N,1, ydat,N);
    X1.initializeValues(N,1,x1dat,N);
    X3.initializeValues(N,1,x3dat,N);
    // solve A*X=Y
    DefaultArithmetic<MV>::Init(Y,1);
    dsm.solve(Teuchos::NO_TRANS,Y,X1);
    TEST_DATA_FOR_ONES(x1dat)
    dsm.solve(Teuchos::CONJ_TRANS,Y,X3);
    TEST_DATA_FOR_ONES(x3dat)
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsMatrix, SparseSolveIdentityUpper, Ordinal, Scalar, Node )
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
    // generate upper triangular identity matrix:
    // [ 1 0           ]
    // [   1 0         ]
    // [     1 0       ]
    // [         . .   ]
    // [            1 0]
    // [              1]
    if (N<2) return;
    RCP<GRPH> G = rcp(new GRPH (N,N,node,null) );
    RCP<MAT>  A = rcp(new MAT  (G,null) );
    // allocate buffers for ptrs, indices and values
    const Ordinal totalNNZ = 2*N-1;
    ArrayRCP<size_t>  ptrs(N+1);
    ArrayRCP<Ordinal> inds(totalNNZ);
    ArrayRCP<Scalar>  vals(totalNNZ);
    // fill the buffers on the host
    {
      Ordinal num = 0;
      for (int i=0; i < N-1; ++i) {
        ptrs[i] = num;
        inds[num] = i; inds[num+1] = i+1;
        vals[num] = 1; vals[num+1] = 0;
        num += 2;
      }
      ptrs[N-1] = num;
      inds[num] = N-1;
      vals[num] = 1;
      num += 1;
      ptrs[N] = num;
    }
    G->setStructure(ptrs, inds);
    ptrs = Teuchos::null;
    inds = Teuchos::null;
    A->setValues(vals);
    vals = Teuchos::null;
    OPS::finalizeGraphAndMatrix(Teuchos::UPPER_TRI,Teuchos::NON_UNIT_DIAG,*G,*A,null);
    Teuchos::EDiag diag;
    Teuchos::EUplo uplo;
    G->getMatDesc(uplo,diag);
    TEST_EQUALITY_CONST( uplo, Teuchos::UPPER_TRI );
    TEST_EQUALITY_CONST( diag, Teuchos::NON_UNIT_DIAG );
    OPS dsm(node);
    out << "Testing with sparse ops: " << Teuchos::typeName(dsm) << std::endl;
    dsm.setGraphAndMatrix(G,A);

    // A is the identity, which allows us to easily test transpose and non-transpose, upper and lower tri
    ArrayRCP<Scalar> ydat, x1dat, x3dat;
    ydat  = node->template allocBuffer<Scalar>(N);
    x1dat = node->template allocBuffer<Scalar>(N);
    x3dat = node->template allocBuffer<Scalar>(N);
    MV Y(node), X1(node), X3(node);
    Y.initializeValues(N,1, ydat,N);
    X1.initializeValues(N,1,x1dat,N);
    X3.initializeValues(N,1,x3dat,N);
    // solve A*X=Y
    DefaultArithmetic<MV>::Init(Y,1);
    dsm.solve(Teuchos::NO_TRANS,Y,X1);
    TEST_DATA_FOR_ONES(x1dat)
    dsm.solve(Teuchos::CONJ_TRANS,Y,X3);
    TEST_DATA_FOR_ONES(x3dat)
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsMatrix, SparseSolveImplicitIdentityUpper, Ordinal, Scalar, Node )
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
    // generate upper triangular identity matrix:
    // [ 1 0           ]
    // [   1 0         ]
    // [     1 0       ]
    // [         . .   ]
    // [            1 0]
    // [              1]
    // but don't store the diagonal
    if (N<2) return;
    RCP<GRPH> G = rcp(new GRPH (N,N,node,null) );
    RCP<MAT>  A = rcp(new MAT  (G,null) );
    // allocate buffers for ptrs, indices and values
    const Ordinal totalNNZ = N-1;
    ArrayRCP<size_t>  ptrs(N+1);
    ArrayRCP<Ordinal> inds(totalNNZ);
    ArrayRCP<Scalar>  vals(totalNNZ);
    // fill the buffers on the host
    {
      Ordinal num = 0;
      for (int i=0; i < N-1; ++i) {
        ptrs[i] = num;
        inds[i] = i+1;
        vals[i] = 0;
        num += 1;
      }
      ptrs[N-1] = num;
      ptrs[N]   = num;
    }
    G->setStructure(ptrs, inds);
    ptrs = Teuchos::null;
    inds = Teuchos::null;
    A->setValues(vals);
    vals = Teuchos::null;
    OPS::finalizeGraphAndMatrix(Teuchos::UPPER_TRI,Teuchos::UNIT_DIAG,*G,*A,null);
    Teuchos::EDiag diag;
    Teuchos::EUplo uplo;
    G->getMatDesc(uplo,diag);
    TEST_EQUALITY_CONST( uplo, Teuchos::UPPER_TRI );
    TEST_EQUALITY_CONST( diag, Teuchos::UNIT_DIAG );
    OPS dsm(node);
    out << "Testing with sparse ops: " << Teuchos::typeName(dsm) << std::endl;
    dsm.setGraphAndMatrix(G,A);

    // A is the identity, which allows us to easily test transpose and non-transpose, upper and lower tri
    ArrayRCP<Scalar> ydat, x1dat, x3dat;
    ydat  = node->template allocBuffer<Scalar>(N);
    x1dat = node->template allocBuffer<Scalar>(N);
    x3dat = node->template allocBuffer<Scalar>(N);
    MV Y(node), X1(node), X3(node);
    Y.initializeValues(N,1, ydat,N);
    X1.initializeValues(N,1,x1dat,N);
    X3.initializeValues(N,1,x3dat,N);
    // solve A*X=Y
    DefaultArithmetic<MV>::Init(Y,1);
    dsm.solve(Teuchos::NO_TRANS,Y,X1);
    TEST_DATA_FOR_ONES(x1dat)
    dsm.solve(Teuchos::CONJ_TRANS,Y,X3);
    TEST_DATA_FOR_ONES(x3dat)
  }

#define ALL_UNIT_TESTS_ORDINAL_SCALAR_NODE( ORDINAL, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsMatrix,          SparseSolveIdentityLower, ORDINAL, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsMatrix,          SparseSolveIdentityUpper, ORDINAL, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsMatrix,  SparseSolveImplicitIdentityLower, ORDINAL, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsMatrix,  SparseSolveImplicitIdentityUpper, ORDINAL, SCALAR, NODE )

#ifdef HAVE_TPETRACLASSIC_SERIAL
  typedef KokkosClassic::DoNotUse::SerialNode KokkosClassic_SerialNode;
#define UNIT_TEST_SERIALNODE(ORDINAL, SCALAR) \
  ALL_UNIT_TESTS_ORDINAL_SCALAR_NODE( ORDINAL, SCALAR, KokkosClassic_SerialNode )
#else
#define UNIT_TEST_SERIALNODE(ORDINAL, SCALAR)
#endif // HAVE_TPETRACLASSIC_SERIAL

#ifdef HAVE_TPETRACLASSIC_TBB
  typedef KokkosClassic::TBBNode KokkosClassic_TBBNode;
#define UNIT_TEST_TBBNODE(ORDINAL, SCALAR) \
  ALL_UNIT_TESTS_ORDINAL_SCALAR_NODE( ORDINAL, SCALAR, KokkosClassic_TBBNode )
#else
#define UNIT_TEST_TBBNODE(ORDINAL, SCALAR)
#endif

#ifdef HAVE_TPETRACLASSIC_OPENMP
  typedef KokkosClassic::OpenMPNode KokkosClassic_OpenMPNode;
#define UNIT_TEST_OPENMPNODE(ORDINAL, SCALAR) \
  ALL_UNIT_TESTS_ORDINAL_SCALAR_NODE( ORDINAL, SCALAR, KokkosClassic_OpenMPNode )
#else
#define UNIT_TEST_OPENMPNODE(ORDINAL, SCALAR)
#endif

#ifdef HAVE_TPETRACLASSIC_THREADPOOL
  typedef KokkosClassic::DoNotUse::TPINode KokkosClassic_TPINode;
#define UNIT_TEST_TPINODE(ORDINAL, SCALAR) \
  ALL_UNIT_TESTS_ORDINAL_SCALAR_NODE( ORDINAL, SCALAR, KokkosClassic_TPINode )
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

  UNIT_TEST_GROUP_ORDINAL(int)

  typedef short int ShortInt;
  UNIT_TEST_GROUP_ORDINAL(ShortInt)

} // namespace (anonymous)

