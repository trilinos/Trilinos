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
#include <Teuchos_ScalarTraits.hpp>

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_DefaultArithmetic.hpp"
#include "Kokkos_DefaultKernels.hpp"
#include "Kokkos_DefaultNode.hpp"
#include "Kokkos_DefaultRelaxation.hpp"
#include "Kokkos_MultiVector.hpp"
#include "Kokkos_Version.hpp"

#define TPI_CHUNKS 4

namespace {

  using KokkosClassic::MultiVector;
  using KokkosClassic::DefaultArithmetic;
  using KokkosClassic::DefaultKernels;
  using KokkosClassic::DefaultRelaxation;
  using Teuchos::ArrayRCP;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::null;
  using std::endl;

  int N = 10;

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

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsMatrix_1D, Jacobi, Ordinal, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    typedef typename DefaultKernels<Scalar,Ordinal,Node>::SparseOps          BASESOPS;
    typedef typename BASESOPS::template bind_scalar<Scalar>::other_type           DSM;
    typedef typename DSM::template graph<Ordinal,Node>::graph_type               GRPH;
    typedef typename DSM::template matrix<Scalar,Ordinal,Node>::matrix_type       MAT;
    typedef MultiVector<Scalar,Node>                                               MV;
    // generate tridiagonal matrix:
    // [ 2 -1                   ]
    // [-1  2  -1               ]
    // [   -1   2  -1           ]
    // [                        ]
    // [                -1  2 -1]
    // [                   -1  2]
    if (N<2) return;
    RCP<GRPH> G = Teuchos::rcp( new GRPH(N,N,node,null) );
    RCP<MAT > A = Teuchos::rcp( new MAT(G,null) );
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
        vals[NNZsofar] =  -1; vals[NNZsofar+1] = 2; vals[NNZsofar+2] =  -1;
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
    DSM::finalizeGraphAndMatrix(Teuchos::UNDEF_TRI,Teuchos::NON_UNIT_DIAG,*G,*A,null);

    int its=10;

    // Allocate Relaxation Object
    DefaultRelaxation<Scalar,Ordinal,Node> dj(node);
    dj.initializeData(G,A);

    // Allocate Vectors
    MV X0(node),RHS(node), AX(node);
    ArrayRCP<Scalar> x0dat  = node->template allocBuffer<Scalar>(N);
    for (int i=0; i<N; ++i) x0dat[i] = i;
    ArrayRCP<Scalar> rhsdat  = node->template allocBuffer<Scalar>(N);
    X0.initializeValues(N,1,x0dat,N);
    RHS.initializeValues(N,1,rhsdat,N);
    Scalar norms,norm0;

    // Set starting vector & run Jacobi
    DefaultArithmetic<MV>::Init(RHS,0);
    norm0=DefaultArithmetic<MV>::Norm2Squared(X0);
    for(int i=0;i<its;i++){
      dj.sweep_jacobi((Scalar)1.0,X0,RHS);
    }
    norms=DefaultArithmetic<MV>::Norm2Squared(X0);
    Scalar goodNorm = 0.280595;
    //Note that this is (||x_10|| / ||x_0|| )^2.
    TEST_FLOATING_EQUALITY(norms/norm0, goodNorm, 1e-6);

  } //Jacobi unit test

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsMatrix_1D, Chebyshev, Ordinal, Scalar, Node )
  {
    out << "Tests Chebyshev relaxation on matrix stored in \"1D\" format." << std::endl;
    RCP<Node> node = getNode<Node>();
    typedef typename DefaultKernels<Scalar,Ordinal,Node>::SparseOps          BASESOPS;
    typedef typename BASESOPS::template bind_scalar<Scalar>::other_type           DSM;
    typedef typename DSM::template graph<Ordinal,Node>::graph_type               GRPH;
    typedef typename DSM::template matrix<Scalar,Ordinal,Node>::matrix_type       MAT;
    typedef MultiVector<Scalar,Node>                                               MV;
    // generate tridiagonal matrix:
    // [ 2 -1                   ]
    // [-1  2  -1               ]
    // [   -1   2  -1           ]
    // [                        ]
    // [                -1  2 -1]
    // [                   -1  2]
    if (N<2) return;
    RCP<GRPH> G = Teuchos::rcp( new GRPH(N,N,node,null) );
    RCP<MAT > A = Teuchos::rcp( new MAT(G,null) );
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
        vals[NNZsofar] =  -1; vals[NNZsofar+1] = 2; vals[NNZsofar+2] =  -1;
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
    DSM::finalizeGraphAndMatrix(Teuchos::UNDEF_TRI,Teuchos::NON_UNIT_DIAG,*G,*A,null);

    int its=10;

    // Allocate Relaxation Object
    DefaultRelaxation<Scalar,Ordinal,Node> dj(node);
    dj.initializeData(G,A);

    // Allocate Vectors
    MV X0(node),RHS(node), AX(node);
    ArrayRCP<Scalar> x0dat  = node->template allocBuffer<Scalar>(N);
    ArrayRCP<Scalar> rhsdat  = node->template allocBuffer<Scalar>(N);
    X0.initializeValues(N,1,x0dat,N);
    RHS.initializeValues(N,1,rhsdat,N);
    typedef  Teuchos::ScalarTraits<Scalar> SCT;
    typedef  typename SCT::magnitudeType   Magnitude;
    Magnitude norms;
    // FIXME (mfh 29 Jun 2012) norm0 gets assigned (below), but never used.
    // gcc 4.6 warns about this; gcc 4.5 doesn't.
    //Magnitude norm0;

    // Set starting vector & run Chebyshev
    DefaultArithmetic<MV>::Init(X0,0);
    DefaultArithmetic<MV>::Init(RHS,1);
    norms=-666;
    // FIXME (mfh 29 Jun 2012) norm0 gets assigned, but never used.
    // gcc 4.6 warns about this; gcc 4.5 doesn't.
    //norm0=DefaultArithmetic<MV>::Norm2Squared(X0);
    norms=DefaultArithmetic<MV>::Norm2Squared(RHS);
    dj.setup_chebyshev((Scalar)1.9594929736,(Scalar) 0.097974648681);
    for(int i=0;i<its;i++){
      dj.sweep_chebyshev(X0,RHS);
    }
    norms=DefaultArithmetic<MV>::Norm2Squared(X0);
    Magnitude goodNorm = 23.209151403;
    TEST_FLOATING_EQUALITY(norms, goodNorm*goodNorm, 1e-6);

    x0dat = null;
    rhsdat= null;
  } //Chebyshev unit test

#ifdef ENABLE_ALL_OTHER_RELAXATION
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsMatrix, HybridRelaxation, Ordinal, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    typedef typename DefaultKernels<Scalar,Ordinal,Node>::SparseOps          BASESOPS;
    typedef typename BASESOPS::template bind_scalar<Scalar>::other_type           DSM;
    typedef typename DSM::template graph<Ordinal,Node>::graph_type               GRPH;
    typedef typename DSM::template matrix<Scalar,Ordinal,Node>::matrix_type       MAT;
    typedef MultiVector<Scalar,Node>                                               MV;
    // generate tridiagonal matrix:
    // [ 2 -1                   ]
    // [-1  2  -1               ]
    // [   -1   2  -1           ]
    // [                        ]
    // [                -1  2 -1]
    // [                   -1  2]
    if (N<2) return;
    RCP<GRPH> G = Teuchos::rcp( new GRPH(N,N,node,null) );
    RCP<MAT > A = Teuchos::rcp( new MAT(G,null) );
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
        vals[NNZsofar] =  -1; vals[NNZsofar+1] = 2; vals[NNZsofar+2] =  -1;
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
    DSM::finalizeGraphAndMatrix(Teuchos::UNDEF_TRI,Teuchos::NON_UNIT_DIAG,*G,*A,null);

    printf("\n");

    int its=10;

    // Allocate Relaxation Object
    DefaultRelaxation<Scalar,Ordinal,Node> dj(node);
    dj.initializeData(G,A);

    // Allocate Vectors
    MV X0(node),RHS(node), AX(node);
    ArrayRCP<Scalar> x0dat  = node->template allocBuffer<Scalar>(N);
    ArrayRCP<Scalar> rhsdat  = node->template allocBuffer<Scalar>(N);
    X0.initializeValues(N,1,x0dat,N);
    RHS.initializeValues(N,1,rhsdat,N);
    Scalar norms,norm0;

    // Set starting vector & run Fine Hybrid
    Teuchos::ScalarTraits<Scalar>::seedrandom(24601);
    DefaultArithmetic<MV>::Random(X0);
    DefaultArithmetic<MV>::Init(RHS,0);
    norms=-666;norm0=DefaultArithmetic<MV>::Norm2Squared(X0);
    printf("*** Fine Hybrid ***\n");
    for(int i=0;i<its;i++){
      dj.sweep_fine_hybrid((Scalar)1.0,X0,RHS);
    }
    norms=DefaultArithmetic<MV>::Norm2Squared(X0);
    printf("[%3d] ||x0|| = %22.16e\n",its,(double)norms/norm0);


    // Set starting vector & run Coarse Hybrid (2)
    Teuchos::ScalarTraits<Scalar>::seedrandom(24601);
    DefaultArithmetic<MV>::Random(X0);
    DefaultArithmetic<MV>::Init(RHS,0);
    norms=-666;norm0=DefaultArithmetic<MV>::Norm2Squared(X0);
    printf("*** Coarse Hybrid (%d chunks) ***\n",2);
    for(int i=0;i<its;i++){
      dj.sweep_coarse_hybrid((Scalar)1.0,2,X0,RHS);
    }
    norms=DefaultArithmetic<MV>::Norm2Squared(X0);
    printf("[%3d] ||x0|| = %22.16e\n",its,(double)norms/norm0);

    // Set starting vector & run Coarse Hybrid (N)
    Teuchos::ScalarTraits<Scalar>::seedrandom(24601);
    DefaultArithmetic<MV>::Random(X0);
    DefaultArithmetic<MV>::Init(RHS,0);
    norms=-666;norm0=DefaultArithmetic<MV>::Norm2Squared(X0);
    int num_chunks=TPI_CHUNKS;
    printf("*** Coarse Hybrid (%d chunks) ***\n",num_chunks);
    for(int i=0;i<its;i++){
      dj.sweep_coarse_hybrid((Scalar)1.0,(size_t)num_chunks,X0,RHS);
    }
    norms=DefaultArithmetic<MV>::Norm2Squared(X0);
    printf("[%3d] ||x0|| = %22.16e\n",its,(double)norms/norm0);

    x0dat = null;
    rhsdat= null;
  } //HybridRelaxation unit test
#endif //ifdef ENABLE_ALL_OTHER_RELAXATION


#define ALL_UNIT_TESTS_ORDINAL_SCALAR_NODE( ORDINAL, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsMatrix_1D, Jacobi, ORDINAL, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsMatrix_1D, Chebyshev, ORDINAL, SCALAR, NODE )

#ifdef HAVE_TPETRACLASSIC_SERIAL
#define UNIT_TEST_SERIALNODE(ORDINAL, SCALAR) \
      ALL_UNIT_TESTS_ORDINAL_SCALAR_NODE( ORDINAL, SCALAR, SerialNode )
#else
#define UNIT_TEST_SERIALNODE(ORDINAL, SCALAR)
#endif // HAVE_TPETRACLASSIC_SERIAL

#ifdef HAVE_TPETRACLASSIC_TBB
#define UNIT_TEST_TBBNODE(ORDINAL, SCALAR) \
      ALL_UNIT_TESTS_ORDINAL_SCALAR_NODE( ORDINAL, SCALAR, TBBNode )
#else
#define UNIT_TEST_TBBNODE(ORDINAL, SCALAR)
#endif

#ifdef HAVE_TPETRACLASSIC_OPENMP
#define UNIT_TEST_OPENMPNODE(ORDINAL, SCALAR) \
      ALL_UNIT_TESTS_ORDINAL_SCALAR_NODE( ORDINAL, SCALAR, OpenMPNode )
#else
#define UNIT_TEST_OPENMPNODE(ORDINAL, SCALAR)
#endif

#ifdef HAVE_TPETRACLASSIC_THREADPOOL
#define UNIT_TEST_TPINODE(ORDINAL, SCALAR) \
      ALL_UNIT_TESTS_ORDINAL_SCALAR_NODE( ORDINAL, SCALAR, TPINode )
#else
#define UNIT_TEST_TPINODE(ORDINAL, SCALAR)
#endif

#define UNIT_TEST_GROUP_ORDINAL_SCALAR( ORDINAL, SCALAR ) \
        UNIT_TEST_SERIALNODE( ORDINAL, SCALAR ) \
        UNIT_TEST_OPENMPNODE( ORDINAL, SCALAR ) \
        UNIT_TEST_TBBNODE( ORDINAL, SCALAR ) \
        UNIT_TEST_TPINODE( ORDINAL, SCALAR )

#define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
        UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, double) \
        UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, float)

  /* Macro to run the actual tests */
     UNIT_TEST_GROUP_ORDINAL(int)

}

