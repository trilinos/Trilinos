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

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_DefaultNode.hpp"
#include "Kokkos_DefaultKernels.hpp"
#include "Kokkos_Version.hpp"

namespace {

  using Kokkos::DefaultNode;
  using Teuchos::ArrayRCP;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::null;
  using Teuchos::arcp;
  using Teuchos::ParameterList;
  using Teuchos::parameterList;

  typedef Kokkos::DefaultNode::DefaultNodeType Node;

  int N = 100;

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
    clp.setOption("test-size",&N,"Vector length for tests.");
  }

  //
  // UNIT TESTS
  // 

  // intialize using static graph
  TEUCHOS_UNIT_TEST( CrsGraph, RuntimeExceptions )
  {
    typedef Kokkos::DefaultKernels<double,int,Node>::SparseOps SparseOps;
    typedef SparseOps::graph<int,Node>::graph_type             Graph;
    const int N = 10;
    RCP<Node> node = Kokkos::DefaultNode::getDefaultNode();
    RCP<Graph> G = rcp(new Graph(N,N,node,parameterList()));
    {
      ArrayRCP<size_t> ptrs_tooSmall(N), ptrs_tooBig(N+2);
      ArrayRCP<int>    inds;
      TEST_THROW( G->setStructure(ptrs_tooSmall, inds), std::runtime_error );
      TEST_THROW( G->setStructure(ptrs_tooBig,   inds), std::runtime_error );
    }
    {
      ArrayRCP<size_t> ptrs(N+1);
      for (int i=0; i<=N; ++i) ptrs[i] = i;
      ArrayRCP<int> tooFewInds(N-1);
      TEST_THROW( G->setStructure(ptrs, tooFewInds), std::runtime_error );
    }
  }

  TEUCHOS_UNIT_TEST( CrsGraph, StatusTests )
  {
    typedef Kokkos::DefaultKernels<double,int,Node>::SparseOps SparseOps;
    typedef SparseOps::graph<int,Node>::graph_type             Graph;
    const int N = 10;
    ArrayRCP<size_t> ptrs = arcp<size_t>(N+1);
    std::fill(ptrs.begin(), ptrs.end(), 0);
    RCP<Node> node = Kokkos::DefaultNode::getDefaultNode();
    RCP<Graph> G = rcp(new Graph(N,N,node,parameterList()));
    G->setStructure(ptrs,null);
    RCP<ParameterList> params = parameterList();
    SparseOps::finalizeGraph(Teuchos::UNDEF_TRI,Teuchos::NON_UNIT_DIAG,*G,params);
    TEST_EQUALITY_CONST( G->isEmpty(), true );
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( CrsMatrix, StaticGraph, Scalar, Ordinal )
  {
    typedef typename Kokkos::DefaultKernels<Scalar,Ordinal,Node>::SparseOps     BaseSparseOps;
    typedef typename BaseSparseOps::template bind_scalar<Scalar>::other_type        SparseOps;
    typedef typename SparseOps::template graph<Ordinal,Node>::graph_type                Graph;
    typedef typename SparseOps::template matrix<Scalar,Ordinal,Node>::matrix_type      Matrix;
    const Ordinal N = 10;
    RCP<Node> node = Kokkos::DefaultNode::getDefaultNode();
    // build a non-empty graph, tridiagonal
    const Ordinal testNumEntries = 3*N-2;
    ArrayRCP<Ordinal> inds(testNumEntries);
    ArrayRCP<Scalar > vals(testNumEntries);
    ArrayRCP<size_t> ptrs(N+1);
    {
      std::fill( inds.begin(), inds.end(), 0 );
      std::fill( vals.begin(), vals.end(), 0 );
      Ordinal curoffset = 0;
      for (Ordinal r=0; r < N; ++r) {
        ptrs[r] = curoffset;
        if (r > 0 && r < N-1) curoffset += 3;
        else                  curoffset += 2;
      }
      ptrs[N] = curoffset;
    }
    RCP<Graph> G = rcp(new Graph(N,N,node,parameterList()));
    TEST_EQUALITY( G->getNumRows(), N );
    TEST_EQUALITY( G->getNumCols(), N );
    G->setStructure(ptrs, inds);
    SparseOps::finalizeGraph(Teuchos::UNDEF_TRI,Teuchos::NON_UNIT_DIAG,*G,null);
    ArrayRCP<const Ordinal> chkInds;
    chkInds = G->getIndices();
    TEST_EQUALITY( inds, chkInds );
    TEST_EQUALITY_CONST( G->isEmpty(), false );
    TEST_EQUALITY( G->getNumRows(), N );
    TEST_EQUALITY( G->getNumCols(), N );
    Matrix M(G,null);
    M.setValues(vals);
    SparseOps::finalizeMatrix(*G,M,null);
  }


  // intialize using static graph
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( CrsMatrix, DynamicGraph, Scalar, Ordinal )
  {
    typedef typename Kokkos::DefaultKernels<Scalar,Ordinal,Node>::SparseOps     BaseSparseOps;
    typedef typename BaseSparseOps::template bind_scalar<Scalar>::other_type        SparseOps;
    typedef typename SparseOps::template graph<Ordinal,Node>::graph_type                Graph;
    typedef typename SparseOps::template matrix<Scalar,Ordinal,Node>::matrix_type      Matrix;
    const Ordinal N = 10;
    RCP<Node> node = Kokkos::DefaultNode::getDefaultNode();
    // build a non-empty graph, tridiagonal
    const Ordinal testNumEntries = 3*N-2;
    ArrayRCP<Ordinal> inds(testNumEntries);
    ArrayRCP<Scalar > vals(testNumEntries);
    ArrayRCP<size_t> ptrs(N+1);
    {
      std::fill( inds.begin(), inds.end(), 0 );
      std::fill( vals.begin(), vals.end(), 0 );
      Ordinal curoffset = 0;
      for (Ordinal r=0; r < N; ++r) {
        Ordinal numper;
        if (r > 0 && r < N-1) numper = 3;
        else numper = 2;
        ptrs[r] = curoffset;
        curoffset += numper;
      }
      ptrs[N] = curoffset;
    }
    RCP<Graph> G = rcp(new Graph(N,N,node,null) );
    Matrix M(G,null);
    G->setStructure(ptrs,inds);
    M.setValues(vals);
    ArrayRCP<const Ordinal> chkInds;
    ArrayRCP<const Scalar> chkVals;
    chkInds = G->getIndices();
    chkVals = M.getValues();
    TEST_EQUALITY( inds, chkInds );
    TEST_EQUALITY( vals, chkVals );
    SparseOps::finalizeGraphAndMatrix(Teuchos::UNDEF_TRI,Teuchos::NON_UNIT_DIAG,*G,M,null);
    TEST_EQUALITY_CONST( G->isEmpty(),     false );
  }

  // no rows
  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( CrsGraph, NoRows, Ordinal )
  {
    typedef typename Kokkos::DefaultKernels<void  ,Ordinal,Node>::SparseOps     BaseSparseOps;
    typedef typename BaseSparseOps::template graph<Ordinal,Node>::graph_type            Graph;
    const Ordinal N = 0;
    RCP<Node> node = Kokkos::DefaultNode::getDefaultNode();

    // build a null-size graph
    Graph G(N,N,node,null);
    {
      ArrayRCP<size_t>  ptrs;
      ArrayRCP<Ordinal> inds;
      ptrs = null;
      inds = null;
      // ptrs is null; not allowed
      TEST_THROW( G.setStructure(ptrs, inds), std::runtime_error );
      ptrs = arcp<size_t>(1);
      ptrs[0] = 0;
      G.setStructure(ptrs, inds);
      ptrs = null;
      inds = null;
      ArrayRCP<const size_t>  getptrs;
      ArrayRCP<const Ordinal> getinds;
      getptrs = G.getPointers();
      getinds = G.getIndices();
      TEST_INEQUALITY_CONST( getptrs, null );
      TEST_EQUALITY_CONST( getinds, null );
    }
  }

 #define UNIT_TEST_GROUP_ORDINAL_SCALAR( ORDINAL, SCALAR ) \
       TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsMatrix, StaticGraph,  SCALAR, ORDINAL ) \
       TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsMatrix, DynamicGraph, SCALAR, ORDINAL )

 #define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
          TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( CrsGraph, NoRows,  ORDINAL ) \
          UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, int) \
          UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, float)
      UNIT_TEST_GROUP_ORDINAL(int)
      typedef short int ShortInt; UNIT_TEST_GROUP_ORDINAL(ShortInt)

}
