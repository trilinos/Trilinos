/*
//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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

#include <Teuchos_TimeMonitor.hpp>

/** \file CrsMatrix_CrsTimingTests.hpp
    \brief A file unit-testing and demonstrating comparison of different CRS object types.
 */

  namespace Test {
    int numRows = 1000;
    int numIters = 1000;
    int numThreads = 0;
  }

  template <class Node>
  RCP<Node> getNode() {
    assert(false);
  }

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
    clp.setOption("test-size",&Test::numRows,"Vector length for tests.");
    clp.setOption("test-iters",&Test::numIters,"Number of mat-vecs to time.");
    clp.setOption("num-threads",&Test::numThreads,"Number of threads.");
  }


  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( CrsTiming, PowerTriDiag, GRPH, MAT )
  {
    typedef typename MAT::ScalarType        Scalar;
    typedef typename MAT::OrdinalType       Ordinal;
    typedef typename MAT::NodeType          Node;
    typedef typename MAT::LocalMatOpsType   SparseOps;
    typedef MultiVector<Scalar,Node>                                 MV;
    typedef Teuchos::ScalarTraits<Scalar>                            ST;
    
    using Teuchos::arcpClone;
    using Teuchos::tuple;

    out << "Testing " << Teuchos::TypeNameTraits<GRPH>::name() << std::endl 
        << "    and " << Teuchos::TypeNameTraits<MAT>::name()  << std::endl;

    RCP<Node> node = getNode<Node>();

    GRPH G(Test::numRows,node);
    MAT  A(G);
    // 
    // allocate buffers and fill the graph and matrix
    // must allocate 2D and finalize(true) to get first-touch reallocation
    // 
    {
      ArrayRCP<size_t> numEntries(Test::numRows);
      ArrayRCP<ArrayRCP<Ordinal> >   inds(Test::numRows);
      ArrayRCP<ArrayRCP<Scalar > >   vals(Test::numRows);
      numEntries[0] = 2;
      inds[0] = arcpClone<Ordinal>( tuple<Ordinal>(0,1) );
      vals[0] = arcpClone<Scalar >( tuple<Scalar >(2,-1) );
      for (int i=1; i != Test::numRows-1; ++i) {
        numEntries[i] = 3;
        inds[i] = arcpClone<Ordinal>( tuple<Ordinal>(i-1,i,i+1) );
        vals[i] = arcpClone<Scalar >( tuple<Scalar >( -1,3, -1) );
      }
      numEntries[Test::numRows-1] = 2;
      inds[Test::numRows-1] = arcpClone<Ordinal>( tuple<Ordinal>(Test::numRows-2,Test::numRows-1) );
      vals[Test::numRows-1] = arcpClone<Scalar> ( tuple<Scalar >(-1,2) );
      G.set2DStructure(inds, numEntries);
      inds    = Teuchos::null;
      A.set2DValues(vals);
      vals    = Teuchos::null;
      const bool PleaseOptimizeStorage = true;
      A.finalize(PleaseOptimizeStorage);
    }
    // 
    // fill the matvec
    // 
    typename SparseOps::template rebind<Scalar>::other matvec(node);
    matvec.initializeStructure(G);
    matvec.initializeValues(A);
    // 
    // time the matvec
    // 
    MV X(node), Y(node);
    X.initializeValues( Test::numRows,1, node->template allocBuffer<Scalar>(Test::numRows), Test::numRows);
    Y.initializeValues( Test::numRows,1, node->template allocBuffer<Scalar>(Test::numRows), Test::numRows);
    DefaultArithmetic<MV>::Init( X, ST::one() );
    DefaultArithmetic<MV>::Init( Y, ST::one() );
    Teuchos::RCP<Teuchos::Time> matvectime = Teuchos::TimeMonitor::getNewTimer("LocalTimer");
    {
      Teuchos::TimeMonitor lcltimer(*matvectime);
      for (int i=0; i<Test::numIters; ++i) {
        if ( (i%2) == 0 ) {
          matvec.multiply(Teuchos::NO_TRANS,ST::one(),X,Y);
        }
        else {
          matvec.multiply(Teuchos::NO_TRANS,ST::one(),Y,X);
        }
      }
    }
    out << "Time is " << matvectime->totalElapsedTime() << std::endl;
  }

