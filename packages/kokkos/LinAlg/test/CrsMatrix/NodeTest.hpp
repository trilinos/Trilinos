/*
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
*/

#ifndef NODE_TEST_HPP_
#define NODE_TEST_HPP_ 

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_ParameterList.hpp>

#include <functional>
#include <algorithm>

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_DefaultSparseOps.hpp"

namespace NodeTest {
  extern int N;
  extern int numIters;
  extern int numThreads;
  extern int verbose;
}

namespace {

  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::Time;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::null;
  using Teuchos::TimeMonitor;
  using Teuchos::ParameterList;
  using Kokkos::MultiVector;
  using Kokkos::DefaultArithmetic;
  using Kokkos::DefaultHostSparseOps;

  template <class NODE>
  RCP<NODE> getNode() {
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"Node type not defined.");
  }

  template <class NODE>
  void initNode() {
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"Node type not defined.");
  }

  //
  // UNIT TESTS
  // 

  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( AAAAA_Is_First, InitNode, NODE )
  {
    out << "Initializing " << Teuchos::TypeNameTraits<NODE>::name() << std::endl;
    initNode<NODE>();
    TEST_EQUALITY(0,0);
  }

  /////////////////////////////////
  /////////////////////////////////
  //  tests
  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( CrsMatrix, TimingTest, BaseOPS )
  {
    typedef double                             Scalar;
    typedef typename BaseOPS::ordinal_type     Ordinal;
    typedef typename BaseOPS::node_type        Node;
    typedef typename BaseOPS::template bind_scalar<Scalar>::other_type       OPS;
    typedef typename OPS::template matrix<Scalar,Ordinal,Node>::matrix_type  MAT;
    typedef typename OPS::template graph<Ordinal,Node>::graph_type          GRPH;
    typedef MultiVector<Scalar,Node>                                          MV;
    typedef Teuchos::ScalarTraits<Scalar>                                     ST;
    const Scalar ONE = ST::one();
    RCP<Node> node = getNode<Node>();
    // generate tridiagonal matrix:
    // [ 2 -1                   ]
    // [-1  3  -1               ]
    // [   -1   3  -1           ]
    // [                        ]
    // [                -1  3 -1]
    // [                   -1  2]
    const int N = NodeTest::N;
    if (N<2) return;
    RCP<GRPH> G = rcp(new GRPH (N,N,node,null) );
    RCP<MAT>  A = rcp(new MAT  (G,null) );
    // allocate data for ptrs, indices and values
    const Ordinal totalNNZ = 3*N - 2;
    ArrayRCP<size_t>  ptrs(N+1);
    ArrayRCP<Ordinal> inds(totalNNZ);
    ArrayRCP<Scalar>  vals(totalNNZ);
    // fill the data
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
    OPS dsm(node);
    out << "Testing with sparse ops: " << Teuchos::typeName(dsm) << std::endl;
    TEST_NOTHROW( dsm.setGraphAndMatrix(G,A) )

    ArrayRCP<Scalar> xdat, axdat;
    xdat  = node->template allocBuffer<Scalar>(N);
    axdat = node->template allocBuffer<Scalar>(N);
    MV X(node), Y(node);
    X.initializeValues( N,1, xdat,N);
    Y.initializeValues(N,1,axdat,N);
    DefaultArithmetic<MV>::Init( X,1);
    Teuchos::RCP<Teuchos::Time> matvectime = Teuchos::TimeMonitor::getNewTimer("LocalTimer");
    {
      Teuchos::TimeMonitor lcltimer(*matvectime);
      for (int i=0; i<NodeTest::numIters; ++i) {
        // Y = A*X
        dsm.multiply(Teuchos::NO_TRANS,ONE,X,Y);
      }
    }
    out << "Time is " << matvectime->totalElapsedTime() << std::endl;
  }


  // 
  // INSTANTIATIONS
  //

  #define TEST_NODE( NODE ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( AAAAA_Is_First, InitNode, NODE ) \
    typedef DefaultHostSparseOps<void,int,NODE,Kokkos::details::DefaultCRSAllocator   >   NoFirstTouch; \
    typedef DefaultHostSparseOps<void,int,NODE,Kokkos::details::FirstTouchCRSAllocator>  YesFirstTouch; \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( CrsMatrix, TimingTest, NoFirstTouch) \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( CrsMatrix, TimingTest, YesFirstTouch)

}

#endif
