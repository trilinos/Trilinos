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

#include "Kokkos_MultiVector.hpp"
#include "Kokkos_DefaultArithmetic.hpp"
#include "Kokkos_VbrMatrix.hpp"
#include "Kokkos_DefaultBlockSparseOps.hpp"
#include "Kokkos_Version.hpp"

#include "Kokkos_SerialNode.hpp"
#ifdef HAVE_KOKKOSCLASSIC_TBB
#include "Kokkos_TBBNode.hpp"
#endif
#ifdef HAVE_KOKKOSCLASSIC_THREADPOOL
#include "Kokkos_TPINode.hpp"
#endif
//#ifdef HAVE_KOKKOSCLASSIC_THRUST
//#include "Kokkos_ThrustGPUNode.hpp"
//#endif

namespace {

  using Kokkos::MultiVector;
  using Kokkos::VbrMatrix;
  using Kokkos::DefaultArithmetic;
  using Kokkos::DefaultBlockSparseOps;
  using Kokkos::SerialNode;
  using Teuchos::ArrayRCP;
  using Teuchos::arcp;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::null;
  using std::endl;

  RCP<SerialNode> snode;
#ifdef HAVE_KOKKOSCLASSIC_TBB
  using Kokkos::TBBNode;
  RCP<TBBNode> tbbnode;
#endif
#ifdef HAVE_KOKKOSCLASSIC_THREADPOOL
  using Kokkos::TPINode;
  RCP<TPINode> tpinode;
#endif
//#ifdef HAVE_KOKKOSCLASSIC_THRUST
//  using Kokkos::ThrustGPUNode;
//  RCP<ThrustGPUNode> thrustnode;
//#endif

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

#ifdef HAVE_KOKKOSCLASSIC_TBB
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

/*#ifdef HAVE_KOKKOSCLASSIC_THRUST
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
  #endif*/

  //
  // UNIT TESTS
  // 

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( VbrMatrix, SparseMultiply1, Ordinal, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    typedef VbrMatrix<Scalar,Ordinal,Node>  VBR;
    typedef MultiVector<Scalar,Node> MV;
    // generate small 2x2 block matrix:
    // [ 1  1    2  2 ]
    // [ 1  1    2  2 ]
    //
    // [ 3  3    4  4 ]
    // [ 3  3    4  4 ]

    // allocate buffers
    const size_t num_point_rows = 4;
    const size_t num_block_rows = 2;
    const size_t num_block_cols = 2;
    const size_t num_block_nz = 4;
    const size_t totalNNZ = 16;
    ArrayRCP<Ordinal> rptr = node->template allocBuffer<Ordinal> (num_block_rows+1);
    ArrayRCP<Ordinal> cptr = node->template allocBuffer<Ordinal> (num_block_cols+1);
    ArrayRCP<size_t> bptr = node->template allocBuffer<size_t> (num_block_rows+1);
    ArrayRCP<Ordinal> bindx = node->template allocBuffer<Ordinal>(num_block_nz);
    ArrayRCP<Ordinal> indx = node->template allocBuffer<Ordinal>(num_block_nz+1);
    ArrayRCP<Scalar>  vals = node->template allocBuffer<Scalar >(totalNNZ);
    // fill the buffers on the host
    {
      ArrayRCP<Ordinal>  rptr_h = node->template viewBufferNonConst<Ordinal>(Kokkos::WriteOnly,num_block_rows+1,rptr);
      ArrayRCP<Ordinal>  cptr_h = node->template viewBufferNonConst<Ordinal>(Kokkos::WriteOnly,num_block_cols+1,cptr);
      ArrayRCP<size_t>  bptr_h = node->template viewBufferNonConst<size_t>(Kokkos::WriteOnly,num_block_rows+1,bptr);
      ArrayRCP<Ordinal>  bindx_h = node->template viewBufferNonConst<Ordinal>(Kokkos::WriteOnly,num_block_nz,bindx);
      ArrayRCP<Ordinal>  indx_h = node->template viewBufferNonConst<Ordinal>(Kokkos::WriteOnly,num_block_nz+1,indx);
      ArrayRCP<Scalar>   vals_h = node->template viewBufferNonConst<Scalar >(Kokkos::WriteOnly,totalNNZ,vals);

      rptr_h[ 0] = 0; rptr_h[ 1] = 2; rptr_h[ 2] = 4;
      cptr_h[ 0] = 0; cptr_h[ 1] = 2; cptr_h[ 2] = 4;
      bptr_h[ 0] = 0; bptr_h[ 1] = 2; bptr_h[ 2] = 4;
      bindx_h[0] = 0; bindx_h[1] = 1; bindx_h[2] = 0; bindx_h[3] = 1;
      indx_h[ 0] = 0; indx_h[ 1] = 4; indx_h[ 2] = 8; indx_h[ 3] = 12; indx_h[4] = 16;

      vals_h[ 0] = 1; vals_h[ 1] = 1; vals_h[ 2] = 1; vals_h[ 3] = 1;
      vals_h[ 4] = 2; vals_h[ 5] = 2; vals_h[ 6] = 2; vals_h[ 7] = 2;
      vals_h[ 8] = 3; vals_h[ 9] = 3; vals_h[10] = 3; vals_h[11] = 3;
      vals_h[12] = 4; vals_h[13] = 4; vals_h[14] = 4; vals_h[15] = 4;
    }
    VBR A(num_block_rows,node);
    A.setPackedValues(vals,rptr,cptr,bptr,bindx,indx);
    DefaultBlockSparseOps<Scalar,Ordinal,Node> dbsm(node);
    dbsm.initializeValues(A);

    ArrayRCP<Scalar> xdat, axdat, ax_check;
    xdat  = node->template allocBuffer<Scalar>(num_point_rows);
    axdat = node->template allocBuffer<Scalar>(num_point_rows);
    ax_check = arcp<Scalar>(num_point_rows);
    ax_check[0] = 6;
    ax_check[1] = 6;
    ax_check[2] = 14;
    ax_check[3] = 14;
    MV X(node), AX(node);
    X.initializeValues( num_point_rows,1, xdat,num_point_rows);
    AX.initializeValues(num_point_rows,1,axdat,num_point_rows);
    DefaultArithmetic<MV>::Init(X,1);
    dbsm.multiply(Teuchos::NO_TRANS,Teuchos::ScalarTraits<Scalar>::one(),X,Teuchos::ScalarTraits<Scalar>::zero(),AX);
    ArrayRCP<const Scalar> axview = node->template viewBuffer<Scalar>(num_point_rows,axdat);
    TEST_COMPARE_FLOATING_ARRAYS(axview, ax_check, Teuchos::ScalarTraits<Scalar>::zero());
    xdat = null;
    axdat = null;
    ax_check = null;

    //now check the case where x and ax have 3 columns instead of 1.
    size_t numVecs = 3;
    xdat  = node->template allocBuffer<Scalar>(num_point_rows*numVecs);
    axdat = node->template allocBuffer<Scalar>(num_point_rows*numVecs);
    ax_check = arcp<Scalar>(num_point_rows*numVecs);
    ax_check[0] = 6;
    ax_check[1] = 6;
    ax_check[2] = 14;
    ax_check[3] = 14;
    ax_check[4] = 6;
    ax_check[5] = 6;
    ax_check[6] = 14;
    ax_check[7] = 14;
    ax_check[8] = 6;
    ax_check[9] = 6;
    ax_check[10] = 14;
    ax_check[11] = 14;
    X.initializeValues( num_point_rows,numVecs, xdat,num_point_rows);
    AX.initializeValues(num_point_rows,numVecs,axdat,num_point_rows);
    DefaultArithmetic<MV>::Init(X,1);
    dbsm.multiply(Teuchos::NO_TRANS,Teuchos::ScalarTraits<Scalar>::one(),X,Teuchos::ScalarTraits<Scalar>::zero(),AX);
    axview = node->template viewBuffer<Scalar>(num_point_rows*numVecs,axdat);
    TEST_COMPARE_FLOATING_ARRAYS(axview, ax_check, Teuchos::ScalarTraits<Scalar>::zero());
    xdat = null;
    axdat = null;
    ax_check = null;
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( VbrMatrix, SparseMultiply2, Ordinal, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    typedef VbrMatrix<Scalar,Ordinal,Node>  VBR;
    typedef MultiVector<Scalar,Node> MV;
    // generate small 2x2 block matrix:
    // [ 1   1  2  2 ]
    // 
    // [ 1   1  2  2 ]
    // [ 3   3  4  4 ]
    // [ 3   3  4  4 ]

    // allocate buffers
    const size_t num_point_rows = 4;
    const size_t num_block_rows = 2;
    const size_t num_block_cols = 2;
    const size_t num_block_nz = 4;
    const size_t totalNNZ = 16;
    ArrayRCP<Ordinal> rptr = node->template allocBuffer<Ordinal> (num_block_rows+1);
    ArrayRCP<Ordinal> cptr = node->template allocBuffer<Ordinal> (num_block_cols+1);
    ArrayRCP<size_t> bptr = node->template allocBuffer<size_t> (num_block_rows+1);
    ArrayRCP<Ordinal> bindx = node->template allocBuffer<Ordinal>(num_block_nz);
    ArrayRCP<Ordinal> indx = node->template allocBuffer<Ordinal>(num_block_nz+1);
    ArrayRCP<Scalar>  vals = node->template allocBuffer<Scalar >(totalNNZ);
    // fill the buffers on the host
    {
      ArrayRCP<Ordinal>  rptr_h = node->template viewBufferNonConst<Ordinal>(Kokkos::WriteOnly,num_block_rows+1,rptr);
      ArrayRCP<Ordinal>  cptr_h = node->template viewBufferNonConst<Ordinal>(Kokkos::WriteOnly,num_block_cols+1,cptr);
      ArrayRCP<size_t>  bptr_h = node->template viewBufferNonConst<size_t>(Kokkos::WriteOnly,num_block_rows+1,bptr);
      ArrayRCP<Ordinal>  bindx_h = node->template viewBufferNonConst<Ordinal>(Kokkos::WriteOnly,num_block_nz,bindx);
      ArrayRCP<Ordinal>  indx_h = node->template viewBufferNonConst<Ordinal>(Kokkos::WriteOnly,num_block_nz+1,indx);
      ArrayRCP<Scalar>   vals_h = node->template viewBufferNonConst<Scalar >(Kokkos::WriteOnly,totalNNZ,vals);

      rptr_h[0] = 0; rptr_h[1] = 1; rptr_h[2] = 4;
      cptr_h[0] = 0; cptr_h[1] = 1; cptr_h[2] = 4;
      bptr_h[0] = 0; bptr_h[1] = 2; bptr_h[2] = 4;
      bindx_h[0] = 0; bindx_h[1] = 1; bindx_h[2] = 0; bindx_h[3] = 1;
      indx_h[0] = 0; indx_h[1] = 1; indx_h[2] = 4; indx_h[3] = 7; indx_h[4] = 16;

      vals_h[0] = 1; vals_h[1] = 1; vals_h[2] = 2; vals_h[3] = 2;
      vals_h[4] = 1; vals_h[5] = 3; vals_h[6] = 3; vals_h[7] = 1;
      vals_h[8] = 3; vals_h[9] = 3; vals_h[10] = 2; vals_h[11] = 4;
      vals_h[12] = 4; vals_h[13] = 2; vals_h[14] = 4; vals_h[15] = 4;

    }
    VBR  A(num_block_rows,node);
    A.setPackedValues(vals,rptr,cptr,bptr,bindx,indx);
    DefaultBlockSparseOps<Scalar,Ordinal,Node> dbsm(node);
    dbsm.initializeValues(A);

    ArrayRCP<Scalar> xdat, axdat, ax_check;
    xdat  = node->template allocBuffer<Scalar>(num_point_rows);
    axdat = node->template allocBuffer<Scalar>(num_point_rows);
    ax_check = arcp<Scalar>(num_point_rows);
    ax_check[0] = 6;
    ax_check[1] = 6;
    ax_check[2] = 14;
    ax_check[3] = 14;
    MV X(node), AX(node);
    X.initializeValues( num_point_rows,1, xdat,num_point_rows);
    AX.initializeValues(num_point_rows,1,axdat,num_point_rows);
    DefaultArithmetic<MV>::Init(X,1);
    dbsm.multiply(Teuchos::NO_TRANS,Teuchos::ScalarTraits<Scalar>::one(),X,Teuchos::ScalarTraits<Scalar>::zero(),AX);
    ArrayRCP<const Scalar> axview = node->template viewBuffer<Scalar>(num_point_rows,axdat);
    TEST_COMPARE_FLOATING_ARRAYS(axview, ax_check, Teuchos::ScalarTraits<Scalar>::zero());
    xdat = null;
    axdat = null;
    ax_check = null;

    //now check the case where x and ax have 3 columns instead of 1.
    size_t numVecs = 3;
    xdat  = node->template allocBuffer<Scalar>(num_point_rows*numVecs);
    axdat = node->template allocBuffer<Scalar>(num_point_rows*numVecs);
    ax_check = arcp<Scalar>(num_point_rows*numVecs);
    ax_check[0] = 6;
    ax_check[1] = 6;
    ax_check[2] = 14;
    ax_check[3] = 14;
    ax_check[4] = 6;
    ax_check[5] = 6;
    ax_check[6] = 14;
    ax_check[7] = 14;
    ax_check[8] = 6;
    ax_check[9] = 6;
    ax_check[10] = 14;
    ax_check[11] = 14;
    X.initializeValues( num_point_rows,numVecs, xdat,num_point_rows);
    AX.initializeValues(num_point_rows,numVecs,axdat,num_point_rows);
    DefaultArithmetic<MV>::Init(X,1);
    dbsm.multiply(Teuchos::NO_TRANS,Teuchos::ScalarTraits<Scalar>::one(),X,Teuchos::ScalarTraits<Scalar>::zero(),AX);
    axview = node->template viewBuffer<Scalar>(num_point_rows*numVecs,axdat);
    TEST_COMPARE_FLOATING_ARRAYS(axview, ax_check, Teuchos::ScalarTraits<Scalar>::zero());
    xdat = null;
    axdat = null;
    ax_check = null;
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( VbrMatrix, SparseMultiply1Transpose, Ordinal, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    typedef VbrMatrix<Scalar,Ordinal,Node>  VBR;
    typedef MultiVector<Scalar,Node> MV;
    // generate small 2x2 block matrix:
    // [ 1  2    2  3 ]
    // [ 1  1    2  2 ]
    //
    // [ 3  4    4  5 ]
    // [ 3  3    4  4 ]

    // allocate buffers
    const size_t num_point_rows = 4;
    const size_t num_block_rows = 2;
    const size_t num_block_cols = 2;
    const size_t num_block_nz = 4;
    const size_t totalNNZ = 16;
    ArrayRCP<Ordinal> rptr = node->template allocBuffer<Ordinal> (num_block_rows+1);
    ArrayRCP<Ordinal> cptr = node->template allocBuffer<Ordinal> (num_block_cols+1);
    ArrayRCP<size_t> bptr = node->template allocBuffer<size_t> (num_block_rows+1);
    ArrayRCP<Ordinal> bindx = node->template allocBuffer<Ordinal>(num_block_nz);
    ArrayRCP<Ordinal> indx = node->template allocBuffer<Ordinal>(num_block_nz+1);
    ArrayRCP<Scalar>  vals = node->template allocBuffer<Scalar >(totalNNZ);
    // fill the buffers on the host
    {
      ArrayRCP<Ordinal>  rptr_h = node->template viewBufferNonConst<Ordinal>(Kokkos::WriteOnly,num_block_rows+1,rptr);
      ArrayRCP<Ordinal>  cptr_h = node->template viewBufferNonConst<Ordinal>(Kokkos::WriteOnly,num_block_cols+1,cptr);
      ArrayRCP<size_t>  bptr_h = node->template viewBufferNonConst<size_t>(Kokkos::WriteOnly,num_block_rows+1,bptr);
      ArrayRCP<Ordinal>  bindx_h = node->template viewBufferNonConst<Ordinal>(Kokkos::WriteOnly,num_block_nz,bindx);
      ArrayRCP<Ordinal>  indx_h = node->template viewBufferNonConst<Ordinal>(Kokkos::WriteOnly,num_block_nz+1,indx);
      ArrayRCP<Scalar>   vals_h = node->template viewBufferNonConst<Scalar >(Kokkos::WriteOnly,totalNNZ,vals);

      rptr_h[0] = 0; rptr_h[1] = 2; rptr_h[2] = 4;
      cptr_h[0] = 0; cptr_h[1] = 2; cptr_h[2] = 4;
      bptr_h[0] = 0; bptr_h[1] = 2; bptr_h[2] = 4;
      bindx_h[0] = 0; bindx_h[1] = 1; bindx_h[2] = 0; bindx_h[3] = 1;
      indx_h[0] = 0; indx_h[1] = 4; indx_h[2] = 8; indx_h[3] = 12; indx_h[4] = 16;

      vals_h[0] = 1; vals_h[1] = 1; vals_h[2] = 2; vals_h[3] = 1;
      vals_h[4] = 2; vals_h[5] = 2; vals_h[6] = 3; vals_h[7] = 2;
      vals_h[8] = 3; vals_h[9] = 3; vals_h[10] = 4; vals_h[11] = 3;
      vals_h[12] = 4; vals_h[13] = 4; vals_h[14] = 5; vals_h[15] = 4;

    }
    VBR  A(num_block_rows,node);
    A.setPackedValues(vals,rptr,cptr,bptr,bindx,indx);
    DefaultBlockSparseOps<Scalar,Ordinal,Node> dbsm(node);
    dbsm.initializeValues(A);

    {
      ArrayRCP<Scalar> ax_check = arcp<Scalar>(num_point_rows);
      ax_check[0] = 8;
      ax_check[1] = 10;
      ax_check[2] = 12;
      ax_check[3] = 14;
      MV X(node), AX(node);
      ArrayRCP<Scalar> xdat  = node->template allocBuffer<Scalar>(num_point_rows);
      X.initializeValues( num_point_rows,1, xdat,num_point_rows);
      ArrayRCP<Scalar> axdat = node->template allocBuffer<Scalar>(num_point_rows);
      AX.initializeValues(num_point_rows,1,axdat,num_point_rows);
      DefaultArithmetic<MV>::Init(X,1);
      dbsm.multiply(Teuchos::TRANS,Teuchos::ScalarTraits<Scalar>::one(),X,Teuchos::ScalarTraits<Scalar>::zero(),AX);
      ArrayRCP<const Scalar> axview = node->template viewBuffer<Scalar>(num_point_rows,axdat);
      TEST_COMPARE_FLOATING_ARRAYS(axview, ax_check, Teuchos::ScalarTraits<Scalar>::zero());
      xdat = null;
      axdat = null;
      ax_check = null;
    }

    //now check the case where x and ax have 3 columns instead of 1.
    {
      const int num_vecs = 3;
      ArrayRCP<Scalar> ax_check = arcp<Scalar>(num_vecs*num_point_rows);
      ax_check[0] = 8;
      ax_check[1] = 10;
      ax_check[2] = 12;
      ax_check[3] = 14;
      std::copy( ax_check.begin(), ax_check.begin()+4, ax_check.begin()+4 );
      std::copy( ax_check.begin(), ax_check.begin()+4, ax_check.begin()+8 );
      //
      MV X(node), AX(node);
      ArrayRCP<Scalar> xdat  = node->template allocBuffer<Scalar>(num_point_rows*num_point_rows);
      X.initializeValues( num_point_rows,num_point_rows, xdat,num_point_rows);
      DefaultArithmetic<MV>::Init(X,1);
      ArrayRCP<Scalar> axdat = node->template allocBuffer<Scalar>(num_point_rows*num_point_rows);
      AX.initializeValues(num_point_rows,num_point_rows,axdat,num_point_rows);
      //
      dbsm.multiply(Teuchos::TRANS,Teuchos::ScalarTraits<Scalar>::one(),X,Teuchos::ScalarTraits<Scalar>::zero(),AX);
      ArrayRCP<const Scalar> axview = node->template viewBuffer<Scalar>(num_vecs*num_point_rows,axdat);
      TEST_COMPARE_FLOATING_ARRAYS(axview, ax_check, Teuchos::ScalarTraits<Scalar>::zero());
      xdat = null;
      axdat = null;
      ax_check = null;
    }
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( VbrMatrix, SparseMultiply2Transpose, Ordinal, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    typedef VbrMatrix<Scalar,Ordinal,Node>  VBR;
    typedef MultiVector<Scalar,Node> MV;
    // generate small 2x2 block matrix:
    // [ 1    2  2  3 ]
    //
    // [ 1    1  2  2 ]
    // [ 3    4  4  5 ]
    // [ 3    3  4  4 ]

    // allocate buffers
    const size_t num_point_rows = 4;
    const size_t num_block_rows = 2;
    const size_t num_block_cols = 2;
    const size_t num_block_nz = 4;
    const size_t totalNNZ = 16;
    ArrayRCP<Ordinal> rptr = node->template allocBuffer<Ordinal> (num_block_rows+1);
    ArrayRCP<Ordinal> cptr = node->template allocBuffer<Ordinal> (num_block_cols+1);
    ArrayRCP<size_t> bptr = node->template allocBuffer<size_t> (num_block_rows+1);
    ArrayRCP<Ordinal> bindx = node->template allocBuffer<Ordinal>(num_block_nz);
    ArrayRCP<Ordinal> indx = node->template allocBuffer<Ordinal>(num_block_nz+1);
    ArrayRCP<Scalar>  vals = node->template allocBuffer<Scalar >(totalNNZ);
    // fill the buffers on the host
    {
      ArrayRCP<Ordinal>  rptr_h = node->template viewBufferNonConst<Ordinal>(Kokkos::WriteOnly,num_block_rows+1,rptr);
      ArrayRCP<Ordinal>  cptr_h = node->template viewBufferNonConst<Ordinal>(Kokkos::WriteOnly,num_block_cols+1,cptr);
      ArrayRCP<size_t>  bptr_h = node->template viewBufferNonConst<size_t>(Kokkos::WriteOnly,num_block_rows+1,bptr);
      ArrayRCP<Ordinal>  bindx_h = node->template viewBufferNonConst<Ordinal>(Kokkos::WriteOnly,num_block_nz,bindx);
      ArrayRCP<Ordinal>  indx_h = node->template viewBufferNonConst<Ordinal>(Kokkos::WriteOnly,num_block_nz+1,indx);
      ArrayRCP<Scalar>   vals_h = node->template viewBufferNonConst<Scalar >(Kokkos::WriteOnly,totalNNZ,vals);

      rptr_h[0] = 0; rptr_h[1] = 1; rptr_h[2] = 4;
      cptr_h[0] = 0; cptr_h[1] = 1; cptr_h[2] = 4;
      bptr_h[0] = 0; bptr_h[1] = 2; bptr_h[2] = 4;
      bindx_h[0] = 0; bindx_h[1] = 1; bindx_h[2] = 0; bindx_h[3] = 1;
      indx_h[0] = 0; indx_h[1] = 1; indx_h[2] = 4; indx_h[3] = 7; indx_h[4] = 16;

      vals_h[0] = 1; vals_h[1] = 2; vals_h[2] = 2; vals_h[3] = 3;
      vals_h[4] = 1; vals_h[5] = 3; vals_h[6] = 3; vals_h[7] = 1;
      vals_h[8] = 4; vals_h[9] = 3; vals_h[10] = 2; vals_h[11] = 4;
      vals_h[12] = 4; vals_h[13] = 2; vals_h[14] = 5; vals_h[15] = 4;

    }
    VBR  A(num_block_rows,node);
    A.setPackedValues(vals,rptr,cptr,bptr,bindx,indx);
    DefaultBlockSparseOps<Scalar,Ordinal,Node> dbsm(node);
    dbsm.initializeValues(A);

    ArrayRCP<Scalar> xdat, axdat, ax_check;
    xdat  = node->template allocBuffer<Scalar>(num_point_rows);
    axdat = node->template allocBuffer<Scalar>(num_point_rows);
    ax_check = arcp<Scalar>(num_point_rows);
    ax_check[0] = 8;
    ax_check[1] = 10;
    ax_check[2] = 12;
    ax_check[3] = 14;
    MV X(node), AX(node);
    X.initializeValues( num_point_rows,1, xdat,num_point_rows);
    AX.initializeValues(num_point_rows,1,axdat,num_point_rows);
    DefaultArithmetic<MV>::Init(X,1);
    dbsm.multiply(Teuchos::TRANS,Teuchos::ScalarTraits<Scalar>::one(),X,Teuchos::ScalarTraits<Scalar>::zero(),AX);
    ArrayRCP<const Scalar> axview = node->template viewBuffer<Scalar>(num_point_rows,axdat);
    TEST_COMPARE_FLOATING_ARRAYS(axview, ax_check, Teuchos::ScalarTraits<Scalar>::zero());
    xdat = null;
    axdat = null;
    ax_check = null;

    //now check the case where x and ax have 3 columns instead of 1.
    size_t numVecs = 3;
    xdat  = node->template allocBuffer<Scalar>(num_point_rows*numVecs);
    axdat = node->template allocBuffer<Scalar>(num_point_rows*numVecs);
    ax_check = arcp<Scalar>(num_point_rows*numVecs);
    ax_check[0] = 8;
    ax_check[1] = 10;
    ax_check[2] = 12;
    ax_check[3] = 14;
    ax_check[4] = 8;
    ax_check[5] = 10;
    ax_check[6] = 12;
    ax_check[7] = 14;
    ax_check[8] = 8;
    ax_check[9] = 10;
    ax_check[10] = 12;
    ax_check[11] = 14;
    X.initializeValues( num_point_rows,numVecs, xdat,num_point_rows);
    AX.initializeValues(num_point_rows,numVecs,axdat,num_point_rows);
    DefaultArithmetic<MV>::Init(X,1);
    dbsm.multiply(Teuchos::TRANS,Teuchos::ScalarTraits<Scalar>::one(),X,Teuchos::ScalarTraits<Scalar>::zero(),AX);
    axview = node->template viewBuffer<Scalar>(num_point_rows*numVecs,axdat);
    TEST_COMPARE_FLOATING_ARRAYS(axview, ax_check, Teuchos::ScalarTraits<Scalar>::zero());
    xdat = null;
    axdat = null;
    ax_check = null;
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( VbrMatrix, SparseMultiply3Transpose, Ordinal, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    typedef VbrMatrix<Scalar,Ordinal,Node>  VBR;
    typedef MultiVector<Scalar,Node> MV;
    // generate small 2x3 block matrix:
    // [ 1    2  2  3   4 ]
    //
    // [ 1    1  2  2   3 ]
    // [ 3    4  4  5   6 ]
    // [ 3    3  4  4   5 ]

    // allocate buffers
    const size_t num_point_rows = 4;
    const size_t num_point_cols = 5;
    const size_t num_block_rows = 2;
    const size_t num_block_cols = 3;
    const size_t num_block_nz = 6;
    const size_t totalNNZ = 20;
    ArrayRCP<Ordinal> rptr = node->template allocBuffer<Ordinal> (num_block_rows+1);
    ArrayRCP<Ordinal> cptr = node->template allocBuffer<Ordinal> (num_block_cols+1);
    ArrayRCP<size_t> bptr = node->template allocBuffer<size_t> (num_block_rows+1);
    ArrayRCP<Ordinal> bindx = node->template allocBuffer<Ordinal>(num_block_nz);
    ArrayRCP<Ordinal> indx = node->template allocBuffer<Ordinal>(num_block_nz+1);
    ArrayRCP<Scalar>  vals = node->template allocBuffer<Scalar >(totalNNZ);
    // fill the buffers on the host
    {
      ArrayRCP<Ordinal>  rptr_h = node->template viewBufferNonConst<Ordinal>(Kokkos::WriteOnly,num_block_rows+1,rptr);
      ArrayRCP<Ordinal>  cptr_h = node->template viewBufferNonConst<Ordinal>(Kokkos::WriteOnly,num_block_cols+1,cptr);
      ArrayRCP<size_t>  bptr_h = node->template viewBufferNonConst<size_t>(Kokkos::WriteOnly,num_block_rows+1,bptr);
      ArrayRCP<Ordinal>  bindx_h = node->template viewBufferNonConst<Ordinal>(Kokkos::WriteOnly,num_block_nz,bindx);
      ArrayRCP<Ordinal>  indx_h = node->template viewBufferNonConst<Ordinal>(Kokkos::WriteOnly,num_block_nz+1,indx);
      ArrayRCP<Scalar>   vals_h = node->template viewBufferNonConst<Scalar >(Kokkos::WriteOnly,totalNNZ,vals);

      rptr_h[0] = 0; rptr_h[1] = 1; rptr_h[2] = 4;
      cptr_h[0] = 0; cptr_h[1] = 1; cptr_h[2] = 4; cptr_h[3] = 5;
      bptr_h[0] = 0; bptr_h[1] = 3; bptr_h[2] = 6;
      bindx_h[0] = 0; bindx_h[1] = 1; bindx_h[2] = 2; bindx_h[3] = 0; bindx_h[4] = 1; bindx_h[5] = 2;
      indx_h[0] = 0; indx_h[1] = 1; indx_h[2] = 4; indx_h[3] = 5; indx_h[4] = 8; indx_h[5] = 17; indx_h[6] = 20;

      vals_h[0] = 1; vals_h[1] = 2; vals_h[2] = 2; vals_h[3] = 3; vals_h[4] = 4;
      vals_h[5] = 1; vals_h[6] = 3; vals_h[7] = 3; vals_h[8] = 1; vals_h[9] = 4;
      vals_h[10] = 3; vals_h[11] = 2; vals_h[12] = 4; vals_h[13] = 4; vals_h[14] = 2;
      vals_h[15] = 5; vals_h[16] = 4; vals_h[17] = 3; vals_h[18] = 6; vals_h[19] = 5;

    }
    VBR  A(num_block_rows,node);
    A.setPackedValues(vals,rptr,cptr,bptr,bindx,indx);
    DefaultBlockSparseOps<Scalar,Ordinal,Node> dbsm(node);
    dbsm.initializeValues(A);

    ArrayRCP<Scalar> xdat, axdat, ax_check;
    xdat  = node->template allocBuffer<Scalar>(num_point_rows);
    axdat = node->template allocBuffer<Scalar>(num_point_cols);
    ax_check = arcp<Scalar>(num_point_cols);
    ax_check[0] = 8;
    ax_check[1] = 10;
    ax_check[2] = 12;
    ax_check[3] = 14;
    ax_check[4] = 18;
    MV X(node), AX(node);
    X.initializeValues( num_point_rows,1, xdat,num_point_rows);
    AX.initializeValues(num_point_cols,1,axdat,num_point_cols);
    DefaultArithmetic<MV>::Init(X,1);
    dbsm.multiply(Teuchos::TRANS,Teuchos::ScalarTraits<Scalar>::one(),X,Teuchos::ScalarTraits<Scalar>::zero(),AX);
    ArrayRCP<const Scalar> axview = node->template viewBuffer<Scalar>(num_point_cols,axdat);
    TEST_COMPARE_FLOATING_ARRAYS(axview, ax_check, Teuchos::ScalarTraits<Scalar>::zero());
    xdat = null;
    axdat = null;
    ax_check = null;

    //now check the case where x and ax have 3 columns instead of 1.
    size_t numVecs = 3;
    xdat  = node->template allocBuffer<Scalar>(num_point_rows*numVecs);
    axdat = node->template allocBuffer<Scalar>(num_point_cols*numVecs);
    ax_check = arcp<Scalar>(num_point_cols*numVecs);
    ax_check[0] = 8;
    ax_check[1] = 10;
    ax_check[2] = 12;
    ax_check[3] = 14;
    ax_check[4] = 18;
    ax_check[5] = 8;
    ax_check[6] = 10;
    ax_check[7] = 12;
    ax_check[8] = 14;
    ax_check[9] = 18;
    ax_check[10] = 8;
    ax_check[11] = 10;
    ax_check[12] = 12;
    ax_check[13] = 14;
    ax_check[14] = 18;
    X.initializeValues( num_point_rows,numVecs, xdat,num_point_rows);
    AX.initializeValues(num_point_cols,numVecs,axdat,num_point_cols);
    DefaultArithmetic<MV>::Init(X,1);
    dbsm.multiply(Teuchos::TRANS,Teuchos::ScalarTraits<Scalar>::one(),X,Teuchos::ScalarTraits<Scalar>::zero(),AX);
    axview = node->template viewBuffer<Scalar>(num_point_cols*numVecs,axdat);
    TEST_COMPARE_FLOATING_ARRAYS(axview, ax_check, Teuchos::ScalarTraits<Scalar>::zero());
    xdat = null;
    axdat = null;
    ax_check = null;
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( VbrMatrix, SolveUpperNoDiag, Ordinal, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    typedef VbrMatrix<Scalar,Ordinal,Node>  VBR;
    typedef MultiVector<Scalar,Node> MV;
    // generate small 2x2 block matrix:
    // [ 0  1    2  2 ]
    // [ 0  0    2  2 ]
    //
    // [         0  4 ]
    // [         0  0 ]

    // allocate buffers
    const size_t num_point_rows = 4;
    const size_t num_block_rows = 2;
    const size_t num_block_cols = 2;
    const size_t num_block_nz = 3;
    const size_t totalNNZ = 12;
    ArrayRCP<Ordinal> rptr = node->template allocBuffer<Ordinal> (num_block_rows+1);
    ArrayRCP<Ordinal> cptr = node->template allocBuffer<Ordinal> (num_block_cols+1);
    ArrayRCP<size_t> bptr = node->template allocBuffer<size_t> (num_block_rows+1);
    ArrayRCP<Ordinal> bindx = node->template allocBuffer<Ordinal>(num_block_nz);
    ArrayRCP<Ordinal> indx = node->template allocBuffer<Ordinal>(num_block_nz+1);
    ArrayRCP<Scalar>  vals = node->template allocBuffer<Scalar >(totalNNZ);
    // fill the buffers on the host
    {
      ArrayRCP<Ordinal>  rptr_h = node->template viewBufferNonConst<Ordinal>(Kokkos::WriteOnly,num_block_rows+1,rptr);
      ArrayRCP<Ordinal>  cptr_h = node->template viewBufferNonConst<Ordinal>(Kokkos::WriteOnly,num_block_cols+1,cptr);
      ArrayRCP<size_t>  bptr_h = node->template viewBufferNonConst<size_t>(Kokkos::WriteOnly,num_block_rows+1,bptr);
      ArrayRCP<Ordinal>  bindx_h = node->template viewBufferNonConst<Ordinal>(Kokkos::WriteOnly,num_block_nz,bindx);
      ArrayRCP<Ordinal>  indx_h = node->template viewBufferNonConst<Ordinal>(Kokkos::WriteOnly,num_block_nz+1,indx);
      ArrayRCP<Scalar>   vals_h = node->template viewBufferNonConst<Scalar >(Kokkos::WriteOnly,totalNNZ,vals);

      rptr_h[ 0] = 0; rptr_h[ 1] = 2; rptr_h[ 2] = 4;
      cptr_h[ 0] = 0; cptr_h[ 1] = 2; cptr_h[ 2] = 4;
      bptr_h[ 0] = 0; bptr_h[ 1] = 2; bptr_h[ 2] = 3;
      bindx_h[0] = 0; bindx_h[1] = 1; bindx_h[2] = 1;
      indx_h[ 0] = 0; indx_h[ 1] = 4; indx_h[ 2] = 8; indx_h[ 3] = 12;

      vals_h[ 0] = 0; vals_h[ 1] = 0; vals_h[ 2] = 1; vals_h[ 3] = 0;
      vals_h[ 4] = 2; vals_h[ 5] = 2; vals_h[ 6] = 2; vals_h[ 7] = 2;
      vals_h[ 8] = 0; vals_h[ 9] = 0; vals_h[10] = 4; vals_h[11] = 0;
    }
    VBR A(num_block_rows,node);
    A.setPackedValues(vals,rptr,cptr,bptr,bindx,indx);
    DefaultBlockSparseOps<Scalar,Ordinal,Node> dbsm(node);
    dbsm.initializeValues(A);

    ArrayRCP<Scalar> xdat, ydat, x_check;
    xdat = node->template allocBuffer<Scalar>(num_point_rows);
    ydat = node->template allocBuffer<Scalar>(num_point_rows);
    {
      ArrayRCP<Scalar> yview = node->template viewBufferNonConst<Scalar>(Kokkos::WriteOnly,ydat.size(),ydat);
      yview[0] = 6;
      yview[1] = 5;
      yview[2] = 5;
      yview[3] = 1;
    }
    x_check = arcp<Scalar>(num_point_rows);
    x_check[0] = 1;
    x_check[1] = 1;
    x_check[2] = 1;
    x_check[3] = 1;
    MV X(node), Y(node);
    X.initializeValues( num_point_rows,1, xdat,num_point_rows);
    Y.initializeValues( num_point_rows,1, ydat,num_point_rows);
    dbsm.solve(Teuchos::NO_TRANS, Teuchos::UPPER_TRI, Teuchos::UNIT_DIAG, Y,X);
    ArrayRCP<const Scalar> xview = node->template viewBuffer<Scalar>(num_point_rows,xdat);
    TEST_COMPARE_FLOATING_ARRAYS(xview, x_check, Teuchos::ScalarTraits<Scalar>::zero());
    xdat = null;
    ydat = null;
    x_check = null;
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( VbrMatrix, SolveUpperUnitDiag, Ordinal, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    typedef VbrMatrix<Scalar,Ordinal,Node>  VBR;
    typedef MultiVector<Scalar,Node> MV;
    // generate small 2x2 block matrix:
    // [ 1  1    2  2 ]
    // [ 0  1    2  2 ]
    //
    // [         1  4 ]
    // [         0  1 ]

    // allocate buffers
    const size_t num_point_rows = 4;
    const size_t num_block_rows = 2;
    const size_t num_block_cols = 2;
    const size_t num_block_nz = 3;
    const size_t totalNNZ = 12;
    ArrayRCP<Ordinal> rptr = node->template allocBuffer<Ordinal> (num_block_rows+1);
    ArrayRCP<Ordinal> cptr = node->template allocBuffer<Ordinal> (num_block_cols+1);
    ArrayRCP<size_t> bptr = node->template allocBuffer<size_t> (num_block_rows+1);
    ArrayRCP<Ordinal> bindx = node->template allocBuffer<Ordinal>(num_block_nz);
    ArrayRCP<Ordinal> indx = node->template allocBuffer<Ordinal>(num_block_nz+1);
    ArrayRCP<Scalar>  vals = node->template allocBuffer<Scalar >(totalNNZ);
    // fill the buffers on the host
    {
      ArrayRCP<Ordinal>  rptr_h = node->template viewBufferNonConst<Ordinal>(Kokkos::WriteOnly,num_block_rows+1,rptr);
      ArrayRCP<Ordinal>  cptr_h = node->template viewBufferNonConst<Ordinal>(Kokkos::WriteOnly,num_block_cols+1,cptr);
      ArrayRCP<size_t>  bptr_h = node->template viewBufferNonConst<size_t>(Kokkos::WriteOnly,num_block_rows+1,bptr);
      ArrayRCP<Ordinal>  bindx_h = node->template viewBufferNonConst<Ordinal>(Kokkos::WriteOnly,num_block_nz,bindx);
      ArrayRCP<Ordinal>  indx_h = node->template viewBufferNonConst<Ordinal>(Kokkos::WriteOnly,num_block_nz+1,indx);
      ArrayRCP<Scalar>   vals_h = node->template viewBufferNonConst<Scalar >(Kokkos::WriteOnly,totalNNZ,vals);

      rptr_h[ 0] = 0; rptr_h[ 1] = 2; rptr_h[ 2] = 4;
      cptr_h[ 0] = 0; cptr_h[ 1] = 2; cptr_h[ 2] = 4;
      bptr_h[ 0] = 0; bptr_h[ 1] = 2; bptr_h[ 2] = 3;
      bindx_h[0] = 0; bindx_h[1] = 1; bindx_h[2] = 1;
      indx_h[ 0] = 0; indx_h[ 1] = 4; indx_h[ 2] = 8; indx_h[ 3] = 12;

      vals_h[ 0] = 1; vals_h[ 1] = 0; vals_h[ 2] = 1; vals_h[ 3] = 1;
      vals_h[ 4] = 2; vals_h[ 5] = 2; vals_h[ 6] = 2; vals_h[ 7] = 2;
      vals_h[ 8] = 1; vals_h[ 9] = 0; vals_h[10] = 4; vals_h[11] = 1;
    }
    VBR A(num_block_rows,node);
    A.setPackedValues(vals,rptr,cptr,bptr,bindx,indx);
    DefaultBlockSparseOps<Scalar,Ordinal,Node> dbsm(node);
    dbsm.initializeValues(A);

    ArrayRCP<Scalar> xdat, ydat, x_check;
    xdat = node->template allocBuffer<Scalar>(num_point_rows);
    ydat = node->template allocBuffer<Scalar>(num_point_rows);
    {
      ArrayRCP<Scalar> yview = node->template viewBufferNonConst<Scalar>(Kokkos::WriteOnly,ydat.size(),ydat);
      yview[0] = 6;
      yview[1] = 5;
      yview[2] = 5;
      yview[3] = 1;
    }
    x_check = arcp<Scalar>(num_point_rows);
    x_check[0] = 1;
    x_check[1] = 1;
    x_check[2] = 1;
    x_check[3] = 1;
    MV X(node), Y(node);
    X.initializeValues( num_point_rows,1, xdat,num_point_rows);
    Y.initializeValues( num_point_rows,1, ydat,num_point_rows);
    dbsm.solve(Teuchos::NO_TRANS, Teuchos::UPPER_TRI, Teuchos::UNIT_DIAG, Y,X);
    ArrayRCP<const Scalar> xview = node->template viewBuffer<Scalar>(num_point_rows,xdat);
    TEST_COMPARE_FLOATING_ARRAYS(xview, x_check, Teuchos::ScalarTraits<Scalar>::zero());
    xdat = null;
    ydat = null;
    x_check = null;
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( VbrMatrix, SolveUpperNonUnitDiag, Ordinal, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    typedef VbrMatrix<Scalar,Ordinal,Node>  VBR;
    typedef MultiVector<Scalar,Node> MV;
    // generate small 2x2 block matrix:
    // [ 1  1    2  2 ]
    // [ 0  2    2  2 ]
    //
    // [         3  4 ]
    // [         0  4 ]

    // allocate buffers
    const size_t num_point_rows = 4;
    const size_t num_block_rows = 2;
    const size_t num_block_cols = 2;
    const size_t num_block_nz = 3;
    const size_t totalNNZ = 12;
    ArrayRCP<Ordinal> rptr = node->template allocBuffer<Ordinal> (num_block_rows+1);
    ArrayRCP<Ordinal> cptr = node->template allocBuffer<Ordinal> (num_block_cols+1);
    ArrayRCP<size_t> bptr = node->template allocBuffer<size_t> (num_block_rows+1);
    ArrayRCP<Ordinal> bindx = node->template allocBuffer<Ordinal>(num_block_nz);
    ArrayRCP<Ordinal> indx = node->template allocBuffer<Ordinal>(num_block_nz+1);
    ArrayRCP<Scalar>  vals = node->template allocBuffer<Scalar >(totalNNZ);
    // fill the buffers on the host
    {
      ArrayRCP<Ordinal>  rptr_h = node->template viewBufferNonConst<Ordinal>(Kokkos::WriteOnly,num_block_rows+1,rptr);
      ArrayRCP<Ordinal>  cptr_h = node->template viewBufferNonConst<Ordinal>(Kokkos::WriteOnly,num_block_cols+1,cptr);
      ArrayRCP<size_t>  bptr_h = node->template viewBufferNonConst<size_t>(Kokkos::WriteOnly,num_block_rows+1,bptr);
      ArrayRCP<Ordinal>  bindx_h = node->template viewBufferNonConst<Ordinal>(Kokkos::WriteOnly,num_block_nz,bindx);
      ArrayRCP<Ordinal>  indx_h = node->template viewBufferNonConst<Ordinal>(Kokkos::WriteOnly,num_block_nz+1,indx);
      ArrayRCP<Scalar>   vals_h = node->template viewBufferNonConst<Scalar >(Kokkos::WriteOnly,totalNNZ,vals);

      rptr_h[ 0] = 0; rptr_h[ 1] = 2; rptr_h[ 2] = 4;
      cptr_h[ 0] = 0; cptr_h[ 1] = 2; cptr_h[ 2] = 4;
      bptr_h[ 0] = 0; bptr_h[ 1] = 2; bptr_h[ 2] = 3;
      bindx_h[0] = 0; bindx_h[1] = 1; bindx_h[2] = 1;
      indx_h[ 0] = 0; indx_h[ 1] = 4; indx_h[ 2] = 8; indx_h[ 3] = 12;

      vals_h[ 0] = 1; vals_h[ 1] = 0; vals_h[ 2] = 1; vals_h[ 3] = 2;
      vals_h[ 4] = 2; vals_h[ 5] = 2; vals_h[ 6] = 2; vals_h[ 7] = 2;
      vals_h[ 8] = 3; vals_h[ 9] = 0; vals_h[10] = 4; vals_h[11] = 4;
    }
    VBR A(num_block_rows,node);
    A.setPackedValues(vals,rptr,cptr,bptr,bindx,indx);
    DefaultBlockSparseOps<Scalar,Ordinal,Node> dbsm(node);
    dbsm.initializeValues(A);

    ArrayRCP<Scalar> xdat, ydat, x_check;
    xdat = node->template allocBuffer<Scalar>(num_point_rows);
    ydat = node->template allocBuffer<Scalar>(num_point_rows);
    {
      ArrayRCP<Scalar> yview = node->template viewBufferNonConst<Scalar>(Kokkos::WriteOnly,ydat.size(),ydat);
      yview[0] = 6;
      yview[1] = 6;
      yview[2] = 7;
      yview[3] = 4;
    }
    x_check = arcp<Scalar>(num_point_rows);
    x_check[0] = 1;
    x_check[1] = 1;
    x_check[2] = 1;
    x_check[3] = 1;
    MV X(node), Y(node);
    X.initializeValues( num_point_rows,1, xdat,num_point_rows);
    Y.initializeValues( num_point_rows,1, ydat,num_point_rows);
    dbsm.solve(Teuchos::NO_TRANS, Teuchos::UPPER_TRI, Teuchos::NON_UNIT_DIAG, Y,X);
    ArrayRCP<const Scalar> xview = node->template viewBuffer<Scalar>(num_point_rows,xdat);
    TEST_COMPARE_FLOATING_ARRAYS(xview, x_check, Teuchos::ScalarTraits<Scalar>::zero());
    xdat = null;
    ydat = null;
    x_check = null;
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( VbrMatrix, SolveTransposeUpperNonUnitDiag, Ordinal, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    typedef VbrMatrix<Scalar,Ordinal,Node>  VBR;
    typedef MultiVector<Scalar,Node> MV;
    // generate small 2x2 block matrix:
    // [ 1  1    2  2 ]
    // [ 0  2    2  2 ]
    //
    // [         3  4 ]
    // [         0  4 ]

    // allocate buffers
    const size_t num_point_rows = 4;
    const size_t num_block_rows = 2;
    const size_t num_block_cols = 2;
    const size_t num_block_nz = 3;
    const size_t totalNNZ = 12;
    ArrayRCP<Ordinal> rptr = node->template allocBuffer<Ordinal> (num_block_rows+1);
    ArrayRCP<Ordinal> cptr = node->template allocBuffer<Ordinal> (num_block_cols+1);
    ArrayRCP<size_t> bptr = node->template allocBuffer<size_t> (num_block_rows+1);
    ArrayRCP<Ordinal> bindx = node->template allocBuffer<Ordinal>(num_block_nz);
    ArrayRCP<Ordinal> indx = node->template allocBuffer<Ordinal>(num_block_nz+1);
    ArrayRCP<Scalar>  vals = node->template allocBuffer<Scalar >(totalNNZ);
    // fill the buffers on the host
    {
      ArrayRCP<Ordinal>  rptr_h = node->template viewBufferNonConst<Ordinal>(Kokkos::WriteOnly,num_block_rows+1,rptr);
      ArrayRCP<Ordinal>  cptr_h = node->template viewBufferNonConst<Ordinal>(Kokkos::WriteOnly,num_block_cols+1,cptr);
      ArrayRCP<size_t>  bptr_h = node->template viewBufferNonConst<size_t>(Kokkos::WriteOnly,num_block_rows+1,bptr);
      ArrayRCP<Ordinal>  bindx_h = node->template viewBufferNonConst<Ordinal>(Kokkos::WriteOnly,num_block_nz,bindx);
      ArrayRCP<Ordinal>  indx_h = node->template viewBufferNonConst<Ordinal>(Kokkos::WriteOnly,num_block_nz+1,indx);
      ArrayRCP<Scalar>   vals_h = node->template viewBufferNonConst<Scalar >(Kokkos::WriteOnly,totalNNZ,vals);

      rptr_h[ 0] = 0; rptr_h[ 1] = 2; rptr_h[ 2] = 4;
      cptr_h[ 0] = 0; cptr_h[ 1] = 2; cptr_h[ 2] = 4;
      bptr_h[ 0] = 0; bptr_h[ 1] = 2; bptr_h[ 2] = 3;
      bindx_h[0] = 0; bindx_h[1] = 1; bindx_h[2] = 1;
      indx_h[ 0] = 0; indx_h[ 1] = 4; indx_h[ 2] = 8; indx_h[ 3] = 12;

      vals_h[ 0] = 1; vals_h[ 1] = 0; vals_h[ 2] = 1; vals_h[ 3] = 2;
      vals_h[ 4] = 2; vals_h[ 5] = 2; vals_h[ 6] = 2; vals_h[ 7] = 2;
      vals_h[ 8] = 3; vals_h[ 9] = 0; vals_h[10] = 4; vals_h[11] = 4;
    }
    VBR A(num_block_rows,node);
    A.setPackedValues(vals,rptr,cptr,bptr,bindx,indx);
    DefaultBlockSparseOps<Scalar,Ordinal,Node> dbsm(node);
    dbsm.initializeValues(A);

    ArrayRCP<Scalar> xdat, ydat, x_check;
    xdat = node->template allocBuffer<Scalar>(num_point_rows);
    ydat = node->template allocBuffer<Scalar>(num_point_rows);
    {
      ArrayRCP<Scalar> yview = node->template viewBufferNonConst<Scalar>(Kokkos::WriteOnly,ydat.size(),ydat);
      yview[0] = 1;
      yview[1] = 3;
      yview[2] = 7;
      yview[3] = 12;
    }
    x_check = arcp<Scalar>(num_point_rows);
    x_check[0] = 1;
    x_check[1] = 1;
    x_check[2] = 1;
    x_check[3] = 1;
    MV X(node), Y(node);
    X.initializeValues( num_point_rows,1, xdat,num_point_rows);
    Y.initializeValues( num_point_rows,1, ydat,num_point_rows);
    dbsm.solve(Teuchos::TRANS, Teuchos::UPPER_TRI, Teuchos::NON_UNIT_DIAG, Y,X);
    ArrayRCP<const Scalar> xview = node->template viewBuffer<Scalar>(num_point_rows,xdat);
    TEST_COMPARE_FLOATING_ARRAYS(xview, x_check, Teuchos::ScalarTraits<Scalar>::zero());
    xdat = null;
    ydat = null;
    x_check = null;
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( VbrMatrix, SolveLowerNoDiag, Ordinal, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    typedef VbrMatrix<Scalar,Ordinal,Node>  VBR;
    typedef MultiVector<Scalar,Node> MV;
    // generate small 2x2 block matrix:
    // [ 0  0         ]
    // [ 1  0         ]
    //
    // [ 2  2    0  0 ]
    // [ 2  2    4  0 ]

    // allocate buffers
    const size_t num_point_rows = 4;
    const size_t num_block_rows = 2;
    const size_t num_block_cols = 2;
    const size_t num_block_nz = 3;
    const size_t totalNNZ = 12;
    ArrayRCP<Ordinal> rptr = node->template allocBuffer<Ordinal> (num_block_rows+1);
    ArrayRCP<Ordinal> cptr = node->template allocBuffer<Ordinal> (num_block_cols+1);
    ArrayRCP<size_t> bptr = node->template allocBuffer<size_t> (num_block_rows+1);
    ArrayRCP<Ordinal> bindx = node->template allocBuffer<Ordinal>(num_block_nz);
    ArrayRCP<Ordinal> indx = node->template allocBuffer<Ordinal>(num_block_nz+1);
    ArrayRCP<Scalar>  vals = node->template allocBuffer<Scalar >(totalNNZ);
    // fill the buffers on the host
    {
      ArrayRCP<Ordinal>  rptr_h = node->template viewBufferNonConst<Ordinal>(Kokkos::WriteOnly,num_block_rows+1,rptr);
      ArrayRCP<Ordinal>  cptr_h = node->template viewBufferNonConst<Ordinal>(Kokkos::WriteOnly,num_block_cols+1,cptr);
      ArrayRCP<size_t>  bptr_h = node->template viewBufferNonConst<size_t>(Kokkos::WriteOnly,num_block_rows+1,bptr);
      ArrayRCP<Ordinal>  bindx_h = node->template viewBufferNonConst<Ordinal>(Kokkos::WriteOnly,num_block_nz,bindx);
      ArrayRCP<Ordinal>  indx_h = node->template viewBufferNonConst<Ordinal>(Kokkos::WriteOnly,num_block_nz+1,indx);
      ArrayRCP<Scalar>   vals_h = node->template viewBufferNonConst<Scalar >(Kokkos::WriteOnly,totalNNZ,vals);

      rptr_h[ 0] = 0; rptr_h[ 1] = 2; rptr_h[ 2] = 4;
      cptr_h[ 0] = 0; cptr_h[ 1] = 2; cptr_h[ 2] = 4;
      bptr_h[ 0] = 0; bptr_h[ 1] = 1; bptr_h[ 2] = 3;
      bindx_h[0] = 0; bindx_h[1] = 0; bindx_h[2] = 1;
      indx_h[ 0] = 0; indx_h[ 1] = 4; indx_h[ 2] = 8; indx_h[ 3] = 12;

      vals_h[ 0] = 0; vals_h[ 1] = 1; vals_h[ 2] = 0; vals_h[ 3] = 0;
      vals_h[ 4] = 2; vals_h[ 5] = 2; vals_h[ 6] = 2; vals_h[ 7] = 2;
      vals_h[ 8] = 0; vals_h[ 9] = 4; vals_h[10] = 0; vals_h[11] = 0;
    }
    VBR A(num_block_rows,node);
    A.setPackedValues(vals,rptr,cptr,bptr,bindx,indx);
    DefaultBlockSparseOps<Scalar,Ordinal,Node> dbsm(node);
    dbsm.initializeValues(A);

    ArrayRCP<Scalar> xdat, ydat, x_check;
    xdat = node->template allocBuffer<Scalar>(num_point_rows);
    ydat = node->template allocBuffer<Scalar>(num_point_rows);
    {
      ArrayRCP<Scalar> yview = node->template viewBufferNonConst<Scalar>(Kokkos::WriteOnly,ydat.size(),ydat);
      yview[0] = 1;
      yview[1] = 2;
      yview[2] = 5;
      yview[3] = 9;
    }
    x_check = arcp<Scalar>(num_point_rows);
    x_check[0] = 1;
    x_check[1] = 1;
    x_check[2] = 1;
    x_check[3] = 1;
    MV X(node), Y(node);
    X.initializeValues( num_point_rows,1, xdat,num_point_rows);
    Y.initializeValues( num_point_rows,1, ydat,num_point_rows);
    dbsm.solve(Teuchos::NO_TRANS, Teuchos::LOWER_TRI, Teuchos::UNIT_DIAG, Y,X);
    ArrayRCP<const Scalar> xview = node->template viewBuffer<Scalar>(num_point_rows,xdat);
    TEST_COMPARE_FLOATING_ARRAYS(xview, x_check, Teuchos::ScalarTraits<Scalar>::zero());
    xdat = null;
    ydat = null;
    x_check = null;
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( VbrMatrix, SolveLowerUnitDiag, Ordinal, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    typedef VbrMatrix<Scalar,Ordinal,Node>  VBR;
    typedef MultiVector<Scalar,Node> MV;
    // generate small 2x2 block matrix:
    // [ 1  0         ]
    // [ 1  1         ]
    //
    // [ 2  2    1  0 ]
    // [ 2  2    4  1 ]

    // allocate buffers
    const size_t num_point_rows = 4;
    const size_t num_block_rows = 2;
    const size_t num_block_cols = 2;
    const size_t num_block_nz = 3;
    const size_t totalNNZ = 12;
    ArrayRCP<Ordinal> rptr = node->template allocBuffer<Ordinal> (num_block_rows+1);
    ArrayRCP<Ordinal> cptr = node->template allocBuffer<Ordinal> (num_block_cols+1);
    ArrayRCP<size_t> bptr = node->template allocBuffer<size_t> (num_block_rows+1);
    ArrayRCP<Ordinal> bindx = node->template allocBuffer<Ordinal>(num_block_nz);
    ArrayRCP<Ordinal> indx = node->template allocBuffer<Ordinal>(num_block_nz+1);
    ArrayRCP<Scalar>  vals = node->template allocBuffer<Scalar >(totalNNZ);
    // fill the buffers on the host
    {
      ArrayRCP<Ordinal>  rptr_h = node->template viewBufferNonConst<Ordinal>(Kokkos::WriteOnly,num_block_rows+1,rptr);
      ArrayRCP<Ordinal>  cptr_h = node->template viewBufferNonConst<Ordinal>(Kokkos::WriteOnly,num_block_cols+1,cptr);
      ArrayRCP<size_t>  bptr_h = node->template viewBufferNonConst<size_t>(Kokkos::WriteOnly,num_block_rows+1,bptr);
      ArrayRCP<Ordinal>  bindx_h = node->template viewBufferNonConst<Ordinal>(Kokkos::WriteOnly,num_block_nz,bindx);
      ArrayRCP<Ordinal>  indx_h = node->template viewBufferNonConst<Ordinal>(Kokkos::WriteOnly,num_block_nz+1,indx);
      ArrayRCP<Scalar>   vals_h = node->template viewBufferNonConst<Scalar >(Kokkos::WriteOnly,totalNNZ,vals);

      rptr_h[ 0] = 0; rptr_h[ 1] = 2; rptr_h[ 2] = 4;
      cptr_h[ 0] = 0; cptr_h[ 1] = 2; cptr_h[ 2] = 4;
      bptr_h[ 0] = 0; bptr_h[ 1] = 1; bptr_h[ 2] = 3;
      bindx_h[0] = 0; bindx_h[1] = 0; bindx_h[2] = 1;
      indx_h[ 0] = 0; indx_h[ 1] = 4; indx_h[ 2] = 8; indx_h[ 3] = 12;

      vals_h[ 0] = 1; vals_h[ 1] = 1; vals_h[ 2] = 0; vals_h[ 3] = 1;
      vals_h[ 4] = 2; vals_h[ 5] = 2; vals_h[ 6] = 2; vals_h[ 7] = 2;
      vals_h[ 8] = 1; vals_h[ 9] = 4; vals_h[10] = 0; vals_h[11] = 1;
    }
    VBR A(num_block_rows,node);
    A.setPackedValues(vals,rptr,cptr,bptr,bindx,indx);
    DefaultBlockSparseOps<Scalar,Ordinal,Node> dbsm(node);
    dbsm.initializeValues(A);

    ArrayRCP<Scalar> xdat, ydat, x_check;
    xdat = node->template allocBuffer<Scalar>(num_point_rows);
    ydat = node->template allocBuffer<Scalar>(num_point_rows);
    {
      ArrayRCP<Scalar> yview = node->template viewBufferNonConst<Scalar>(Kokkos::WriteOnly,ydat.size(),ydat);
      yview[0] = 1;
      yview[1] = 2;
      yview[2] = 5;
      yview[3] = 9;
    }
    x_check = arcp<Scalar>(num_point_rows);
    x_check[0] = 1;
    x_check[1] = 1;
    x_check[2] = 1;
    x_check[3] = 1;
    MV X(node), Y(node);
    X.initializeValues( num_point_rows,1, xdat,num_point_rows);
    Y.initializeValues( num_point_rows,1, ydat,num_point_rows);
    dbsm.solve(Teuchos::NO_TRANS, Teuchos::LOWER_TRI, Teuchos::UNIT_DIAG, Y,X);
    ArrayRCP<const Scalar> xview = node->template viewBuffer<Scalar>(num_point_rows,xdat);
    TEST_COMPARE_FLOATING_ARRAYS(xview, x_check, Teuchos::ScalarTraits<Scalar>::zero());
    xdat = null;
    ydat = null;
    x_check = null;
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( VbrMatrix, SolveLowerNonUnitDiag, Ordinal, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    typedef VbrMatrix<Scalar,Ordinal,Node>  VBR;
    typedef MultiVector<Scalar,Node> MV;
    // generate small 2x2 block matrix:
    // [ 1  0         ]
    // [ 1  2         ]
    //
    // [ 2  2    3  0 ]
    // [ 2  2    4  4 ]

    // allocate buffers
    const size_t num_point_rows = 4;
    const size_t num_block_rows = 2;
    const size_t num_block_cols = 2;
    const size_t num_block_nz = 3;
    const size_t totalNNZ = 12;
    ArrayRCP<Ordinal> rptr = node->template allocBuffer<Ordinal> (num_block_rows+1);
    ArrayRCP<Ordinal> cptr = node->template allocBuffer<Ordinal> (num_block_cols+1);
    ArrayRCP<size_t> bptr = node->template allocBuffer<size_t> (num_block_rows+1);
    ArrayRCP<Ordinal> bindx = node->template allocBuffer<Ordinal>(num_block_nz);
    ArrayRCP<Ordinal> indx = node->template allocBuffer<Ordinal>(num_block_nz+1);
    ArrayRCP<Scalar>  vals = node->template allocBuffer<Scalar >(totalNNZ);
    // fill the buffers on the host
    {
      ArrayRCP<Ordinal>  rptr_h = node->template viewBufferNonConst<Ordinal>(Kokkos::WriteOnly,num_block_rows+1,rptr);
      ArrayRCP<Ordinal>  cptr_h = node->template viewBufferNonConst<Ordinal>(Kokkos::WriteOnly,num_block_cols+1,cptr);
      ArrayRCP<size_t>  bptr_h = node->template viewBufferNonConst<size_t>(Kokkos::WriteOnly,num_block_rows+1,bptr);
      ArrayRCP<Ordinal>  bindx_h = node->template viewBufferNonConst<Ordinal>(Kokkos::WriteOnly,num_block_nz,bindx);
      ArrayRCP<Ordinal>  indx_h = node->template viewBufferNonConst<Ordinal>(Kokkos::WriteOnly,num_block_nz+1,indx);
      ArrayRCP<Scalar>   vals_h = node->template viewBufferNonConst<Scalar >(Kokkos::WriteOnly,totalNNZ,vals);

      rptr_h[ 0] = 0; rptr_h[ 1] = 2; rptr_h[ 2] = 4;
      cptr_h[ 0] = 0; cptr_h[ 1] = 2; cptr_h[ 2] = 4;
      bptr_h[ 0] = 0; bptr_h[ 1] = 1; bptr_h[ 2] = 3;
      bindx_h[0] = 0; bindx_h[1] = 0; bindx_h[2] = 1;
      indx_h[ 0] = 0; indx_h[ 1] = 4; indx_h[ 2] = 8; indx_h[ 3] = 12;

      vals_h[ 0] = 1; vals_h[ 1] = 1; vals_h[ 2] = 0; vals_h[ 3] = 2;
      vals_h[ 4] = 2; vals_h[ 5] = 2; vals_h[ 6] = 2; vals_h[ 7] = 2;
      vals_h[ 8] = 3; vals_h[ 9] = 4; vals_h[10] = 0; vals_h[11] = 4;
    }
    VBR A(num_block_rows,node);
    A.setPackedValues(vals,rptr,cptr,bptr,bindx,indx);
    DefaultBlockSparseOps<Scalar,Ordinal,Node> dbsm(node);
    dbsm.initializeValues(A);

    ArrayRCP<Scalar> xdat, ydat, x_check;
    xdat = node->template allocBuffer<Scalar>(num_point_rows);
    ydat = node->template allocBuffer<Scalar>(num_point_rows);
    {
      ArrayRCP<Scalar> yview = node->template viewBufferNonConst<Scalar>(Kokkos::WriteOnly,ydat.size(),ydat);
      yview[0] = 1;
      yview[1] = 3;
      yview[2] = 7;
      yview[3] = 12;
    }
    x_check = arcp<Scalar>(num_point_rows);
    x_check[0] = 1;
    x_check[1] = 1;
    x_check[2] = 1;
    x_check[3] = 1;
    MV X(node), Y(node);
    X.initializeValues( num_point_rows,1, xdat,num_point_rows);
    Y.initializeValues( num_point_rows,1, ydat,num_point_rows);
    dbsm.solve(Teuchos::NO_TRANS, Teuchos::LOWER_TRI, Teuchos::NON_UNIT_DIAG, Y,X);
    ArrayRCP<const Scalar> xview = node->template viewBuffer<Scalar>(num_point_rows,xdat);
    TEST_COMPARE_FLOATING_ARRAYS(xview, x_check, Teuchos::ScalarTraits<Scalar>::zero());
    xdat = null;
    ydat = null;
    x_check = null;
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( VbrMatrix, SolveTransposeLowerNonUnitDiag, Ordinal, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    typedef VbrMatrix<Scalar,Ordinal,Node>  VBR;
    typedef MultiVector<Scalar,Node> MV;
    // generate small 2x2 block matrix:
    // [ 1  0         ]
    // [ 1  2         ]
    //
    // [ 2  2    3  0 ]
    // [ 2  2    4  4 ]

    // allocate buffers
    const size_t num_point_rows = 4;
    const size_t num_block_rows = 2;
    const size_t num_block_cols = 2;
    const size_t num_block_nz = 3;
    const size_t totalNNZ = 12;
    ArrayRCP<Ordinal> rptr = node->template allocBuffer<Ordinal> (num_block_rows+1);
    ArrayRCP<Ordinal> cptr = node->template allocBuffer<Ordinal> (num_block_cols+1);
    ArrayRCP<size_t> bptr = node->template allocBuffer<size_t> (num_block_rows+1);
    ArrayRCP<Ordinal> bindx = node->template allocBuffer<Ordinal>(num_block_nz);
    ArrayRCP<Ordinal> indx = node->template allocBuffer<Ordinal>(num_block_nz+1);
    ArrayRCP<Scalar>  vals = node->template allocBuffer<Scalar >(totalNNZ);
    // fill the buffers on the host
    {
      ArrayRCP<Ordinal>  rptr_h = node->template viewBufferNonConst<Ordinal>(Kokkos::WriteOnly,num_block_rows+1,rptr);
      ArrayRCP<Ordinal>  cptr_h = node->template viewBufferNonConst<Ordinal>(Kokkos::WriteOnly,num_block_cols+1,cptr);
      ArrayRCP<size_t>  bptr_h = node->template viewBufferNonConst<size_t>(Kokkos::WriteOnly,num_block_rows+1,bptr);
      ArrayRCP<Ordinal>  bindx_h = node->template viewBufferNonConst<Ordinal>(Kokkos::WriteOnly,num_block_nz,bindx);
      ArrayRCP<Ordinal>  indx_h = node->template viewBufferNonConst<Ordinal>(Kokkos::WriteOnly,num_block_nz+1,indx);
      ArrayRCP<Scalar>   vals_h = node->template viewBufferNonConst<Scalar >(Kokkos::WriteOnly,totalNNZ,vals);

      rptr_h[ 0] = 0; rptr_h[ 1] = 2; rptr_h[ 2] = 4;
      cptr_h[ 0] = 0; cptr_h[ 1] = 2; cptr_h[ 2] = 4;
      bptr_h[ 0] = 0; bptr_h[ 1] = 1; bptr_h[ 2] = 3;
      bindx_h[0] = 0; bindx_h[1] = 0; bindx_h[2] = 1;
      indx_h[ 0] = 0; indx_h[ 1] = 4; indx_h[ 2] = 8; indx_h[ 3] = 12;

      vals_h[ 0] = 1; vals_h[ 1] = 1; vals_h[ 2] = 0; vals_h[ 3] = 2;
      vals_h[ 4] = 2; vals_h[ 5] = 2; vals_h[ 6] = 2; vals_h[ 7] = 2;
      vals_h[ 8] = 3; vals_h[ 9] = 4; vals_h[10] = 0; vals_h[11] = 4;
    }
    VBR A(num_block_rows,node);
    A.setPackedValues(vals,rptr,cptr,bptr,bindx,indx);
    DefaultBlockSparseOps<Scalar,Ordinal,Node> dbsm(node);
    dbsm.initializeValues(A);

    ArrayRCP<Scalar> xdat, ydat, x_check;
    xdat = node->template allocBuffer<Scalar>(num_point_rows);
    ydat = node->template allocBuffer<Scalar>(num_point_rows);
    {
      ArrayRCP<Scalar> yview = node->template viewBufferNonConst<Scalar>(Kokkos::WriteOnly,ydat.size(),ydat);
      yview[0] = 6;
      yview[1] = 6;
      yview[2] = 7;
      yview[3] = 4;
    }
    x_check = arcp<Scalar>(num_point_rows);
    x_check[0] = 1;
    x_check[1] = 1;
    x_check[2] = 1;
    x_check[3] = 1;
    MV X(node), Y(node);
    X.initializeValues( num_point_rows,1, xdat,num_point_rows);
    Y.initializeValues( num_point_rows,1, ydat,num_point_rows);
    dbsm.solve(Teuchos::TRANS, Teuchos::LOWER_TRI, Teuchos::NON_UNIT_DIAG, Y,X);
    ArrayRCP<const Scalar> xview = node->template viewBuffer<Scalar>(num_point_rows,xdat);
    TEST_COMPARE_FLOATING_ARRAYS(xview, x_check, Teuchos::ScalarTraits<Scalar>::zero());
    xdat = null;
    ydat = null;
    x_check = null;
  }

#define ALL_UNIT_TESTS_ORDINAL_SCALAR_NODE( ORDINAL, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( VbrMatrix, SparseMultiply1, ORDINAL, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( VbrMatrix, SparseMultiply2, ORDINAL, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( VbrMatrix, SparseMultiply1Transpose, ORDINAL, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( VbrMatrix, SparseMultiply2Transpose, ORDINAL, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( VbrMatrix, SparseMultiply3Transpose, ORDINAL, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( VbrMatrix, SolveUpperNoDiag, ORDINAL, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( VbrMatrix, SolveUpperUnitDiag, ORDINAL, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( VbrMatrix, SolveUpperNonUnitDiag, ORDINAL, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( VbrMatrix, SolveTransposeUpperNonUnitDiag, ORDINAL, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( VbrMatrix, SolveLowerNoDiag, ORDINAL, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( VbrMatrix, SolveLowerUnitDiag, ORDINAL, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( VbrMatrix, SolveLowerNonUnitDiag, ORDINAL, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( VbrMatrix, SolveTransposeLowerNonUnitDiag, ORDINAL, SCALAR, NODE )

#define UNIT_TEST_SERIALNODE(ORDINAL, SCALAR) \
      ALL_UNIT_TESTS_ORDINAL_SCALAR_NODE( ORDINAL, SCALAR, SerialNode )

#ifdef HAVE_KOKKOSCLASSIC_TBB
#define UNIT_TEST_TBBNODE(ORDINAL, SCALAR) \
      ALL_UNIT_TESTS_ORDINAL_SCALAR_NODE( ORDINAL, SCALAR, TBBNode )
#else
#define UNIT_TEST_TBBNODE(ORDINAL, SCALAR)
#endif

#ifdef HAVE_KOKKOSCLASSIC_THREADPOOL
#define UNIT_TEST_TPINODE(ORDINAL, SCALAR) \
      ALL_UNIT_TESTS_ORDINAL_SCALAR_NODE( ORDINAL, SCALAR, TPINode )
#else
#define UNIT_TEST_TPINODE(ORDINAL, SCALAR)
#endif

/*#ifdef HAVE_KOKKOSCLASSIC_THRUST
  #define UNIT_TEST_THRUSTGPUNODE(ORDINAL, SCALAR) \
        ALL_UNIT_TESTS_ORDINAL_SCALAR_NODE( ORDINAL, SCALAR, ThrustGPUNode )
  #else
  #define UNIT_TEST_THRUSTGPUNODE(ORDINAL, SCALAR)
  #endif*/

#define UNIT_TEST_GROUP_ORDINAL_SCALAR( ORDINAL, SCALAR ) \
        UNIT_TEST_SERIALNODE( ORDINAL, SCALAR ) \
        UNIT_TEST_TBBNODE( ORDINAL, SCALAR ) \
        UNIT_TEST_TPINODE( ORDINAL, SCALAR ) 
        //UNIT_TEST_THRUSTGPUNODE( ORDINAL, SCALAR )

#define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
        UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, int) \
        UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, float)

     UNIT_TEST_GROUP_ORDINAL(int)
     typedef short int ShortInt; UNIT_TEST_GROUP_ORDINAL(ShortInt)

}

