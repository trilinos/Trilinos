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

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_DefaultNode.hpp"
#include "Kokkos_DefaultKernels.hpp"
#include "Kokkos_CrsMatrix.hpp"
#include "Kokkos_CrsGraph.hpp"
#include "Kokkos_Version.hpp"

namespace {

  // static graph * matrix-owned graph
  // 1D * 2D
  // OptimizeStorage * none

  // empty

  // CrsGraph exception throwing

  // graph/matrix mismatch throwing

  using Kokkos::DefaultNode;
  using Kokkos::CrsMatrix;
  using Kokkos::CrsGraph;
  using Teuchos::ArrayRCP;
  using Teuchos::RCP;
  using Teuchos::null;
  using Teuchos::arcp;

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
    typedef CrsGraph<int,Node,SparseOps>                       Graph;
    const size_t N = 10;
    RCP<Node> node = Kokkos::DefaultNode::getDefaultNode();
    Graph G(N,node);
    // set1DStructure array sizes
    {
      ArrayRCP<size_t>  b_tooSmall(N), b_justRight(N+1), b_tooBig(N+2);
      ArrayRCP<size_t>  e_tooSmall(N-1), e_justRight(N), e_tooBig(N+1);
      ArrayRCP<int>     inds;
      TEST_THROW( G.set1DStructure(inds,b_tooSmall, e_tooSmall),  std::runtime_error );
      TEST_THROW( G.set1DStructure(inds,b_tooSmall, e_justRight), std::runtime_error );
      TEST_THROW( G.set1DStructure(inds,b_tooSmall, e_tooBig),    std::runtime_error );
      TEST_THROW( G.set1DStructure(inds,b_justRight,e_tooSmall), std::runtime_error );
      TEST_THROW( G.set1DStructure(inds,b_justRight,e_tooBig),   std::runtime_error );
      TEST_THROW( G.set1DStructure(inds,b_tooBig,   e_tooSmall),    std::runtime_error );
      TEST_THROW( G.set1DStructure(inds,b_tooBig,   e_justRight),   std::runtime_error );
      TEST_THROW( G.set1DStructure(inds,b_tooBig,   e_tooBig),      std::runtime_error );
    }
    {
      ArrayRCP<size_t> begs(N+1), ends(N);
      for (size_t i=0; i<N; ++i) {begs[i] = i; ends[i] = i+1;}
      begs[N] = N;
      ArrayRCP<int> tooFewInds(N-1);
      TEST_THROW( G.set1DStructure(tooFewInds, begs,ends), std::runtime_error );
    }
    // set2DStructure array sizes
    {
      ArrayRCP<ArrayRCP<int> > inds(N+1);
      ArrayRCP<size_t> numEntries(N+1);
      TEST_THROW( G.set2DStructure(inds.persistingView(0,N-1), numEntries.persistingView(0,N)), std::runtime_error );
      TEST_THROW( G.set2DStructure(inds.persistingView(0,N+1), numEntries.persistingView(0,N)), std::runtime_error );
      TEST_THROW( G.set2DStructure(inds.persistingView(0,N), numEntries.persistingView(0,N+1)), std::runtime_error );
      TEST_THROW( G.set2DStructure(inds.persistingView(0,N), numEntries.persistingView(0,N-1)), std::runtime_error );
    }
  }

  // intialize using static graph
  TEUCHOS_UNIT_TEST( CrsGraph, StatusTests )
  {
    typedef Kokkos::DefaultKernels<double,int,Node>::SparseOps SparseOps;
    typedef CrsGraph<int,Node,SparseOps>                       Graph;
    const size_t N = 10;
    RCP<Node> node = Kokkos::DefaultNode::getDefaultNode();
    Graph G(N,node);
    G.finalize(false);
    TEST_EQUALITY_CONST( G.isEmpty(), true );
    TEST_EQUALITY_CONST( G.isOptimized(), false );
    TEST_EQUALITY_CONST( G.is1DStructure(), false );
    TEST_EQUALITY_CONST( G.is2DStructure(), false );
    G.finalize(true);
    TEST_EQUALITY_CONST( G.isEmpty(), true );
    TEST_EQUALITY_CONST( G.isOptimized(), false );
    TEST_EQUALITY_CONST( G.is1DStructure(), false );
    TEST_EQUALITY_CONST( G.is2DStructure(), false );
  }

  // intialize using static graph
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( CrsMatrix, StaticGraph, Scalar, Ordinal )
  {
    typedef typename Kokkos::DefaultKernels<Scalar,Ordinal,Node>::SparseOps SparseOps;
    typedef CrsGraph<Ordinal,Node,SparseOps>                                    Graph;
    typedef CrsMatrix<Scalar,Ordinal,Node,SparseOps>                           Matrix;
    const size_t N = 10;
    RCP<Node> node = Kokkos::DefaultNode::getDefaultNode();
    {
      // static graph must be finalized when passed to matrix (empty is okay for this test)
      Graph G(N,node);
      const Graph &staticG = G;
      TEST_THROW( Matrix M1(staticG), std::runtime_error );
      G.finalize(false);  // false, true, doesn't matter here
      TEST_NOTHROW( Matrix M2(staticG) );
    }
    // build a non-empty graph, tridiagonal
    const size_t testNumEntries = 3*N-2;
    ArrayRCP<Ordinal> inds1D(testNumEntries);
    ArrayRCP<Scalar > vals1D(testNumEntries);
    ArrayRCP<ArrayRCP<Ordinal> > inds2D(N);
    ArrayRCP<ArrayRCP<Scalar > > vals2D(N);
    ArrayRCP<size_t> begs(N+1), ends(N), numper(N);
    {
      std::fill( inds1D.begin(), inds1D.end(), 0 );
      std::fill( vals1D.begin(), vals1D.end(), 0 );
      size_t curoffset = 0;
      for (size_t r=0; r < N; ++r) {
        if (r > 0 && r < N-1) numper[r] = 3;
        else numper[r] = 2;
        begs[r] = curoffset;
        ends[r] = begs[r] + numper[r];
        inds2D[r] = inds1D.persistingView(begs[r],numper[r]);
        vals2D[r] = vals1D.persistingView(begs[r],numper[r]);
        curoffset += numper[r];
      }
      begs[N] = curoffset;
    }
    Graph G(N,node);
    TEST_EQUALITY( G.getNumRows(), N );
    TEST_EQUALITY( G.getNode(), node );
    for (int t=0; t < 4; ++t)
    {
      const bool submit1D        = (t && 1);
      const bool OptimizeStorage = (t && 2);
      if (submit1D) {
        G.set1DStructure(inds1D, begs, ends);
        TEST_EQUALITY_CONST( G.is1DStructure(), true );
        TEST_EQUALITY_CONST( G.is2DStructure(), false );
        ArrayRCP<Ordinal> chkInds;
        ArrayRCP<size_t> chkBegs, chkEnds;
        G.get1DStructure(chkInds, chkBegs, chkEnds);
        TEST_EQUALITY( inds1D, chkInds );
        TEST_EQUALITY( begs, chkBegs );
        TEST_EQUALITY( ends, chkEnds );
      }
      else {
        G.set2DStructure(inds2D, numper);
        TEST_EQUALITY_CONST( G.is2DStructure(), true );
        TEST_EQUALITY_CONST( G.is1DStructure(), false );
        ArrayRCP<ArrayRCP<Ordinal> > chkInds;
        ArrayRCP<size_t> chkNumPer;
        G.get2DStructure(chkInds, chkNumPer);
        TEST_EQUALITY( inds2D, chkInds );
        TEST_EQUALITY( numper, chkNumPer );
      }
      TEST_EQUALITY_CONST( G.isFinalized(), false );
      TEST_EQUALITY_CONST( G.isOptimized(), false );
      TEST_EQUALITY_CONST( G.getNumEntries(), testNumEntries );
      G.finalize(OptimizeStorage);
      TEST_EQUALITY_CONST( G.isFinalized(), true );
      TEST_EQUALITY_CONST( G.isEmpty(), false );
      if (OptimizeStorage) {
        TEST_EQUALITY_CONST( G.isOptimized(), true );
        TEST_EQUALITY_CONST( G.is1DStructure(), true );
        TEST_EQUALITY_CONST( G.is2DStructure(), false );
        ArrayRCP<Ordinal> chkInds1D;
        ArrayRCP<ArrayRCP<Ordinal> > chkInds2D;
        ArrayRCP<size_t> chkBegs, chkEnds, chkNumPer;
        G.get1DStructure(chkInds1D, chkBegs, chkEnds);
        G.get2DStructure(chkInds2D, chkNumPer);
        TEST_EQUALITY_CONST( chkInds1D == Teuchos::null, false );
        TEST_EQUALITY_CONST( chkBegs   == Teuchos::null, false );
        TEST_EQUALITY_CONST( chkEnds   == Teuchos::null, false );
        TEST_EQUALITY_CONST( chkInds2D == Teuchos::null, true );
        TEST_EQUALITY_CONST( chkNumPer == Teuchos::null, true );
      }
      else {
        TEST_EQUALITY_CONST( G.isOptimized(), false );
        TEST_EQUALITY( G.is1DStructure(),  submit1D );
        TEST_EQUALITY( G.is2DStructure(), !submit1D );
        ArrayRCP<Ordinal> chkInds1D;
        ArrayRCP<ArrayRCP<Ordinal> > chkInds2D;
        ArrayRCP<size_t> chkBegs, chkEnds, chkNumPer;
        G.get1DStructure(chkInds1D, chkBegs, chkEnds);
        G.get2DStructure(chkInds2D, chkNumPer);
        TEST_EQUALITY( chkInds1D == Teuchos::null, !submit1D );
        TEST_EQUALITY( chkBegs   == Teuchos::null, !submit1D );
        TEST_EQUALITY( chkEnds   == Teuchos::null, !submit1D );
        TEST_EQUALITY( chkInds2D == Teuchos::null,  submit1D );
        TEST_EQUALITY( chkNumPer == Teuchos::null,  submit1D );
      }
      const Graph &staticG = G;
      TEST_EQUALITY( staticG.getNumRows(), N );
      TEST_EQUALITY( staticG.getNode(), node );
      Matrix M(staticG);
      TEST_EQUALITY( M.getNumRows(), G.getNumRows() );
      TEST_EQUALITY( M.getNumEntries(), G.getNumEntries() );
      if (G.is1DStructure()) {
        M.set1DValues(vals1D);
        TEST_EQUALITY_CONST( M.is1DStructure(), true );
        TEST_EQUALITY_CONST( M.is2DStructure(), false );
      }
      else {
        M.set2DValues(vals2D);
        TEST_EQUALITY_CONST( M.is2DStructure(), true );
        TEST_EQUALITY_CONST( M.is1DStructure(), false );
      }
      TEST_EQUALITY_CONST( M.isEmpty(), false );        // can only query this pre-finalize for a static graph
      TEST_EQUALITY_CONST( M.isFinalized(), false );
      if (G.isOptimized() == false) {
        TEST_THROW( M.finalize(true), std::runtime_error );
        M.finalize(false);
      }
      else {
        // true or false doesn't matter; switch it up (submit1D alternates w.r.t. OptimizeStorge)
        M.finalize( submit1D );
      }
      TEST_EQUALITY_CONST( M.isEmpty(), false );
      TEST_EQUALITY_CONST( M.isFinalized(), true );
      TEST_EQUALITY_CONST( M.isOptimized(), G.isOptimized() );
      {
        // test graph clear
        G.clear();
        TEST_EQUALITY_CONST( G.isFinalized(), false );
        TEST_EQUALITY( G.getNumRows(), N );
        TEST_EQUALITY( G.getNumEntries(), 0 );
        TEST_EQUALITY_CONST( G.is1DStructure(), false );
        TEST_EQUALITY_CONST( G.is2DStructure(), false );
        ArrayRCP<Ordinal> chkInds1D;
        ArrayRCP<ArrayRCP<Ordinal> > chkInds2D;
        ArrayRCP<size_t> chkBegs, chkEnds, chkNumPer;
        G.get1DStructure(chkInds1D, chkBegs, chkEnds);
        G.get2DStructure(chkInds2D, chkNumPer);
        TEST_EQUALITY_CONST( chkInds1D == Teuchos::null, true );
        TEST_EQUALITY_CONST( chkBegs   == Teuchos::null, true );
        TEST_EQUALITY_CONST( chkEnds   == Teuchos::null, true );
        TEST_EQUALITY_CONST( chkInds2D == Teuchos::null, true );
        TEST_EQUALITY_CONST( chkNumPer == Teuchos::null, true );
      }
    }
  }

  // intialize using static graph
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( CrsMatrix, DynamicGraph, Scalar, Ordinal )
  {
    typedef typename Kokkos::DefaultKernels<Scalar,Ordinal,Node>::SparseOps SparseOps;
    typedef CrsGraph<Ordinal,Node,SparseOps>                                   Graph;
    typedef CrsMatrix<Scalar,Ordinal,Node,SparseOps>                           Matrix;
    const size_t N = 10;
    RCP<Node> node = Kokkos::DefaultNode::getDefaultNode();

    // build a non-empty graph, tridiagonal
    const size_t testNumEntries = 3*N-2;
    ArrayRCP<Ordinal> inds1D(testNumEntries);
    ArrayRCP<Scalar > vals1D(testNumEntries);
    ArrayRCP<ArrayRCP<Ordinal> > inds2D(N);
    ArrayRCP<ArrayRCP<Scalar > > vals2D(N);
    ArrayRCP<size_t> begs(N+1), ends(N), numper(N);
    {
      std::fill( inds1D.begin(), inds1D.end(), 0 );
      std::fill( vals1D.begin(), vals1D.end(), 0 );
      size_t curoffset = 0;
      for (size_t r=0; r < N; ++r) {
        if (r > 0 && r < N-1) numper[r] = 3;
        else numper[r] = 2;
        begs[r] = curoffset;
        ends[r] = begs[r] + numper[r];
        inds2D[r] = inds1D.persistingView(begs[r],numper[r]);
        vals2D[r] = vals1D.persistingView(begs[r],numper[r]);
        curoffset += numper[r];
      }
      begs[N] = curoffset;
    }
    Graph G(N,node);
    Matrix M(G);
    for (int t=0; t < 4; ++t) {
      const bool submit1D        = (t && 1);
      const bool OptimizeStorage = (t && 2);
      if (submit1D) {
        G.set1DStructure(inds1D, begs, ends);
        M.set1DValues(vals1D);
        TEST_EQUALITY_CONST( G.is1DStructure(), true );
        TEST_EQUALITY_CONST( G.is2DStructure(), false );
        TEST_EQUALITY_CONST( M.is1DStructure(), true );
        TEST_EQUALITY_CONST( M.is2DStructure(), false );
        ArrayRCP<Ordinal> chkInds;
        ArrayRCP<size_t> chkBegs, chkEnds;
        ArrayRCP<Scalar> chkVals;
        G.get1DStructure(chkInds, chkBegs, chkEnds);
        M.get1DValues(chkVals);
        TEST_EQUALITY( inds1D, chkInds );
        TEST_EQUALITY( begs, chkBegs );
        TEST_EQUALITY( ends, chkEnds );
        TEST_EQUALITY( vals1D, chkVals );
      }
      else {
        G.set2DStructure(inds2D, numper);
        M.set2DValues(vals2D);
        TEST_EQUALITY_CONST( G.is2DStructure(), true );
        TEST_EQUALITY_CONST( G.is1DStructure(), false );
        TEST_EQUALITY_CONST( M.is2DStructure(), true );
        TEST_EQUALITY_CONST( M.is1DStructure(), false );
        ArrayRCP<ArrayRCP<Ordinal> > chkInds;
        ArrayRCP<ArrayRCP<Scalar> > chkVals;
        ArrayRCP<size_t> chkNumPer;
        G.get2DStructure(chkInds, chkNumPer);
        M.get2DValues(chkVals);
        TEST_EQUALITY( inds2D, chkInds );
        TEST_EQUALITY( numper, chkNumPer );
        TEST_EQUALITY( vals2D, chkVals );
      }
      TEST_EQUALITY_CONST( G.isFinalized(), false );
      TEST_EQUALITY_CONST( G.isOptimized(), false );
      TEST_EQUALITY_CONST( G.getNumEntries(), testNumEntries );
      TEST_EQUALITY_CONST( M.isFinalized(), false );
      TEST_EQUALITY_CONST( M.isOptimized(), false );
      M.finalize(OptimizeStorage);
      TEST_EQUALITY_CONST( G.isFinalized(), true );
      TEST_EQUALITY_CONST( G.isEmpty(), false );
      TEST_EQUALITY_CONST( M.isFinalized(), true );
      TEST_EQUALITY_CONST( M.isEmpty(), false );
      if (OptimizeStorage) {
        TEST_EQUALITY_CONST( G.isOptimized(), true );
        TEST_EQUALITY_CONST( G.is1DStructure(), true );
        TEST_EQUALITY_CONST( G.is2DStructure(), false );
        TEST_EQUALITY_CONST( M.isOptimized(), true );
        TEST_EQUALITY_CONST( M.is1DStructure(), true );
        TEST_EQUALITY_CONST( M.is2DStructure(), false );
        ArrayRCP<Ordinal> chkInds1D;
        ArrayRCP<ArrayRCP<Ordinal> > chkInds2D;
        ArrayRCP<ArrayRCP<Scalar> >  chkVals2D;
        ArrayRCP<size_t> chkBegs, chkEnds, chkNumPer;
        ArrayRCP<Scalar> chkVals1D;
        G.get1DStructure(chkInds1D, chkBegs, chkEnds);
        G.get2DStructure(chkInds2D, chkNumPer);
        M.get1DValues(chkVals1D);
        M.get2DValues(chkVals2D);
        TEST_EQUALITY_CONST( chkInds1D == Teuchos::null, false );
        TEST_EQUALITY_CONST( chkVals1D == Teuchos::null, false );
        TEST_EQUALITY_CONST( chkBegs   == Teuchos::null, false );
        TEST_EQUALITY_CONST( chkEnds   == Teuchos::null, false );
        TEST_EQUALITY_CONST( chkInds2D == Teuchos::null, true );
        TEST_EQUALITY_CONST( chkVals2D == Teuchos::null, true );
        TEST_EQUALITY_CONST( chkNumPer == Teuchos::null, true );
      }
      else {
        TEST_EQUALITY_CONST( G.isOptimized(), false );
        TEST_EQUALITY( G.is1DStructure(),  submit1D );
        TEST_EQUALITY( G.is2DStructure(), !submit1D );
        TEST_EQUALITY_CONST( M.isOptimized(), false );
        TEST_EQUALITY( M.is1DStructure(),  submit1D );
        TEST_EQUALITY( M.is2DStructure(), !submit1D );
        ArrayRCP<Ordinal> chkInds1D;
        ArrayRCP<Scalar>  chkVals1D;
        ArrayRCP<ArrayRCP<Ordinal> > chkInds2D;
        ArrayRCP<ArrayRCP<Scalar> >  chkVals2D;
        ArrayRCP<size_t> chkBegs, chkEnds, chkNumPer;
        G.get1DStructure(chkInds1D, chkBegs, chkEnds);
        G.get2DStructure(chkInds2D, chkNumPer);
        M.get1DValues(chkVals1D);
        M.get2DValues(chkVals2D);
        TEST_EQUALITY( chkInds1D == Teuchos::null, !submit1D );
        TEST_EQUALITY( chkVals1D == Teuchos::null, !submit1D );
        TEST_EQUALITY( chkBegs   == Teuchos::null, !submit1D );
        TEST_EQUALITY( chkEnds   == Teuchos::null, !submit1D );
        TEST_EQUALITY( chkInds2D == Teuchos::null,  submit1D );
        TEST_EQUALITY( chkVals2D == Teuchos::null,  submit1D );
        TEST_EQUALITY( chkNumPer == Teuchos::null,  submit1D );
      }
      {
        // test matrix and graph clear
        M.clear();
        TEST_EQUALITY( M.getNumRows(), N );
        TEST_EQUALITY( G.getNumRows(), N );
        TEST_EQUALITY_CONST( G.getNumEntries() != 0, true );
        TEST_EQUALITY_CONST( G.is1DStructure() || G.is2DStructure(), true );
        TEST_EQUALITY_CONST( M.isFinalized(), false );
        TEST_EQUALITY_CONST( G.isFinalized(), true );
        G.clear();
        TEST_EQUALITY( G.getNumRows(), N );
        TEST_EQUALITY_CONST( G.getNumEntries() == 0, true );
        TEST_EQUALITY_CONST( G.is1DStructure(), false );
        TEST_EQUALITY_CONST( G.is2DStructure(), false );
        TEST_EQUALITY_CONST( G.isFinalized(), false );
        ArrayRCP<Ordinal> chkInds1D;
        ArrayRCP<Scalar>  chkVals1D;
        ArrayRCP<ArrayRCP<Ordinal> > chkInds2D;
        ArrayRCP<ArrayRCP<Scalar> >  chkVals2D;
        ArrayRCP<size_t> chkBegs, chkEnds, chkNumPer;
        G.get1DStructure(chkInds1D, chkBegs, chkEnds);
        G.get2DStructure(chkInds2D, chkNumPer);
        M.get1DValues(chkVals1D);
        M.get2DValues(chkVals2D);
        TEST_EQUALITY_CONST( chkInds1D == Teuchos::null, true );
        TEST_EQUALITY_CONST( chkVals1D == Teuchos::null, true );
        TEST_EQUALITY_CONST( chkBegs   == Teuchos::null, true );
        TEST_EQUALITY_CONST( chkEnds   == Teuchos::null, true );
        TEST_EQUALITY_CONST( chkInds2D == Teuchos::null, true );
        TEST_EQUALITY_CONST( chkVals2D == Teuchos::null, true );
        TEST_EQUALITY_CONST( chkNumPer == Teuchos::null, true );
      }
    }
  }

  // bad ends
  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( CrsGraph, BadEnds, Ordinal )
  {
    typedef typename Kokkos::DefaultKernels<float,Ordinal,Node>::SparseOps SparseOps;
    typedef CrsGraph<Ordinal,Node,SparseOps>                                   Graph;
    const size_t N = 1;
    RCP<Node> node = Kokkos::DefaultNode::getDefaultNode();

    // build a null-size graph
    Graph G(N,node);
    {
      ArrayRCP<size_t> begs, ends;
      ArrayRCP<Ordinal> inds1D;
      begs = arcp<size_t>(N+1);
      ends = arcp<size_t>(N);
      //
#ifdef HAVE_KOKKOS_DEBUG
      begs[0] = 0; begs[1] = 0;
      ends[0] = 1;  // ends causes rows to overlap; not consistent
      TEST_THROW( G.set1DStructure(inds1D, begs, ends), std::runtime_error );
#endif
      //
      begs[0] = 0; begs[1] = 2; // begs size is larger than allocation in inds1D
      ends[0] = 0;
      inds1D = arcp<Ordinal>(1);
      TEST_THROW( G.set1DStructure(inds1D, begs, ends), std::runtime_error );
      // 
      begs[0] = 0; begs[1] = 1;
      ends[0] = 0;
      // this allocation is allowed to be too large, w.r.t. begs
      inds1D = arcp<Ordinal>(2);
      G.set1DStructure(inds1D, begs, ends);
    }
  }

  // no rows
  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( CrsGraph, NoRows, Ordinal )
  {
    typedef typename Kokkos::DefaultKernels<float,Ordinal,Node>::SparseOps SparseOps;
    typedef CrsGraph<Ordinal,Node,SparseOps>                                   Graph;
    const size_t N = 0;
    RCP<Node> node = Kokkos::DefaultNode::getDefaultNode();

    // build a null-size graph
    Graph G(N,node);
    {
      ArrayRCP<size_t> begs, ends;
      ArrayRCP<Ordinal> inds1D;
      begs = null;
      ends = null;
      inds1D = null;
      // begs is null; not allowed
      TEST_THROW( G.set1DStructure(inds1D, begs, ends), std::runtime_error );
      begs = arcp<size_t>(1);
      begs[0] = 0;
      G.set1DStructure(inds1D, begs, ends);
      begs = null;
      ends = null;
      inds1D = null;
      G.get1DStructure(inds1D, begs, ends);
      TEST_EQUALITY_CONST( begs == null, false );
      TEST_EQUALITY_CONST( ends == null, true );
      TEST_EQUALITY_CONST( inds1D == null, true );
    }
  }

 #define UNIT_TEST_GROUP_ORDINAL_SCALAR( ORDINAL, SCALAR ) \
       TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsMatrix, StaticGraph,  SCALAR, ORDINAL ) \
       TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsMatrix, DynamicGraph, SCALAR, ORDINAL )

 #define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
          TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( CrsGraph, NoRows,  ORDINAL ) \
          TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( CrsGraph, BadEnds,  ORDINAL ) \
          UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, int) \
          UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, float)
      UNIT_TEST_GROUP_ORDINAL(int)
      typedef short int ShortInt; UNIT_TEST_GROUP_ORDINAL(ShortInt)

}
