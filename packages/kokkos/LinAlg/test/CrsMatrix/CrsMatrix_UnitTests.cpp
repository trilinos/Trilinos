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
#include "Kokkos_CrsMatrix.hpp"
#include "Kokkos_CrsGraph.hpp"
#include "Kokkos_Version.hpp"

namespace {

  using Kokkos::DefaultNode;
  using Kokkos::CrsMatrix;
  using Kokkos::CrsGraph;
  using Teuchos::ArrayRCP;
  using Teuchos::RCP;

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

  // intialize using packed data
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( CrsMatrix, PackedData, Scalar, Ordinal )
  {
    const int N = 10;
    typedef typename Node::size_t size;
    RCP<Node> node = Kokkos::DefaultNode::getDefaultNode();
    // test non-empty
    {
      CrsGraph<Ordinal,Node> G(N,node);
      CrsMatrix<Scalar,Node> A(N,node);
      TEST_EQUALITY( G.getNumRows(), N );
      TEST_EQUALITY( A.getNumRows(), N );
      TEST_EQUALITY_CONST( G.isPacked(), false );
      TEST_EQUALITY_CONST( A.isPacked(), false );
      TEST_EQUALITY_CONST( G.isEmpty(),  true );
      TEST_EQUALITY_CONST( A.isEmpty(),  true );
      //
      ArrayRCP<size_t> offsets = node->template allocBuffer<size_t>(N+1);
      ArrayRCP<Ordinal> inds   = node->template allocBuffer<Ordinal>(2*N);
      ArrayRCP<Scalar>  vals   = node->template allocBuffer<Scalar>(2*N);
      G.setPackedStructure(offsets,inds);
      A.setPackedValues(vals);
      TEST_EQUALITY_CONST( G.isPacked(), true );
      TEST_EQUALITY_CONST( A.isPacked(), true );
      TEST_EQUALITY_CONST( G.isEmpty(), false );
      TEST_EQUALITY_CONST( A.isEmpty(), false );
      TEST_EQUALITY_CONST( G.getPackedOffsets(), offsets );
      TEST_EQUALITY_CONST( G.getPackedIndices(), inds );
      TEST_EQUALITY_CONST( A.getPackedValues(),  vals );
      //
      A.clear();
      G.clear();
      TEST_EQUALITY( G.getNumRows(), N );
      TEST_EQUALITY( A.getNumRows(), N );
      TEST_EQUALITY_CONST( G.isPacked(), false );
      TEST_EQUALITY_CONST( A.isPacked(), false );
      TEST_EQUALITY_CONST( G.isEmpty(),  true );
      TEST_EQUALITY_CONST( A.isEmpty(),  true );
    }
    // test empty
    {
      CrsGraph<Ordinal,Node> G(N,node);
      CrsMatrix<Scalar,Node> A(N,node);
      TEST_EQUALITY( G.getNumRows(), N );
      TEST_EQUALITY( A.getNumRows(), N );
      TEST_EQUALITY_CONST( G.isPacked(), false );
      TEST_EQUALITY_CONST( A.isPacked(), false );
      TEST_EQUALITY_CONST( G.isEmpty(),  true );
      TEST_EQUALITY_CONST( A.isEmpty(),  true );
      //
      ArrayRCP<size_t> offsets;
      ArrayRCP<Ordinal> inds;
      ArrayRCP<Scalar>  vals;
      G.setPackedStructure(offsets,inds);
      A.setPackedValues(vals);
      TEST_EQUALITY_CONST( G.isPacked(), false );
      TEST_EQUALITY_CONST( A.isPacked(), false );
      TEST_EQUALITY_CONST( G.isEmpty(),  true );
      TEST_EQUALITY_CONST( A.isEmpty(),  true );
      TEST_EQUALITY_CONST( G.getPackedOffsets(), offsets );
      TEST_EQUALITY_CONST( G.getPackedIndices(), inds );
      TEST_EQUALITY_CONST( A.getPackedValues(),  vals );
      //
      A.clear();
      G.clear();
      TEST_EQUALITY( G.getNumRows(), N );
      TEST_EQUALITY( A.getNumRows(), N );
      TEST_EQUALITY_CONST( G.isPacked(), false );
      TEST_EQUALITY_CONST( A.isPacked(), false );
      TEST_EQUALITY_CONST( G.isEmpty(),  true );
      TEST_EQUALITY_CONST( A.isEmpty(),  true );
    }
  }

  // intialize using nonpacked data
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( CrsMatrix, NonPackedData, Scalar, Ordinal )
  {
    const int N = 3;
    typedef typename Node::size_t size;
    RCP<Node> node = Kokkos::DefaultNode::getDefaultNode();
    // test sort of empty
    {
      CrsGraph<Ordinal,Node> G(N,node);
      CrsMatrix<Scalar,Node> A(N,node);
      TEST_EQUALITY( G.getNumRows(), N );
      TEST_EQUALITY( A.getNumRows(), N );
      TEST_EQUALITY_CONST( G.isPacked(), false );
      TEST_EQUALITY_CONST( A.isPacked(), false );
      TEST_EQUALITY_CONST( G.isEmpty(),  true );
      TEST_EQUALITY_CONST( A.isEmpty(),  true );
      //
      ArrayRCP<Ordinal> inds   = node->template allocBuffer<Ordinal>(1);
      ArrayRCP<Scalar>  vals   = node->template allocBuffer<Scalar>(1);
      G.set2DIndices(1,inds);
      A.set2DValues(1,vals);
      TEST_EQUALITY_CONST( G.isPacked(), false );
      TEST_EQUALITY_CONST( A.isPacked(), false );
      TEST_EQUALITY_CONST( G.isEmpty(),  false );
      TEST_EQUALITY_CONST( A.isEmpty(),  false );
      TEST_EQUALITY_CONST( G.get2DIndices(0), Teuchos::null );
      TEST_EQUALITY_CONST( A.get2DValues( 0), Teuchos::null );
      TEST_EQUALITY_CONST( G.get2DIndices(1), inds );
      TEST_EQUALITY_CONST( A.get2DValues( 1), vals );
      TEST_EQUALITY_CONST( G.get2DIndices(2), Teuchos::null );
      TEST_EQUALITY_CONST( A.get2DValues( 2), Teuchos::null );
      //
      A.clear();
      G.clear();
      TEST_EQUALITY( G.getNumRows(), N );
      TEST_EQUALITY( A.getNumRows(), N );
      TEST_EQUALITY_CONST( G.isPacked(), false );
      TEST_EQUALITY_CONST( A.isPacked(), false );
      TEST_EQUALITY_CONST( G.isEmpty(),  true );
      TEST_EQUALITY_CONST( A.isEmpty(),  true );
    }
    // test empty
    {
      CrsGraph<Ordinal,Node> G(N,node);
      CrsMatrix<Scalar,Node> A(N,node);
      TEST_EQUALITY( G.getNumRows(), N );
      TEST_EQUALITY( A.getNumRows(), N );
      TEST_EQUALITY_CONST( G.isPacked(), false );
      TEST_EQUALITY_CONST( A.isPacked(), false );
      TEST_EQUALITY_CONST( G.isEmpty(),  true );
      TEST_EQUALITY_CONST( A.isEmpty(),  true );
      // 
      ArrayRCP<Ordinal> inds;
      ArrayRCP<Scalar>  vals;
      G.set2DIndices(0,inds);
      G.set2DIndices(1,inds);
      G.set2DIndices(2,inds);
      A.set2DValues(0,vals);
      A.set2DValues(1,vals);
      A.set2DValues(2,vals);
      TEST_EQUALITY_CONST( G.isPacked(), false );
      TEST_EQUALITY_CONST( A.isPacked(), false );
      TEST_EQUALITY_CONST( G.isEmpty(),  true );
      TEST_EQUALITY_CONST( A.isEmpty(),  true );
      for (int i=0; i<3; ++i) {
        TEST_EQUALITY_CONST( G.get2DIndices(i), Teuchos::null );
        TEST_EQUALITY_CONST( A.get2DValues( i), Teuchos::null );
      }
      //
      A.clear();
      G.clear();
      TEST_EQUALITY( G.getNumRows(), N );
      TEST_EQUALITY( A.getNumRows(), N );
      TEST_EQUALITY_CONST( G.isPacked(), false );
      TEST_EQUALITY_CONST( A.isPacked(), false );
      TEST_EQUALITY_CONST( G.isEmpty(),  true );
      TEST_EQUALITY_CONST( A.isEmpty(),  true );
    }
  }

 #define UNIT_TEST_GROUP_ORDINAL_SCALAR( ORDINAL, SCALAR ) \
       TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsMatrix,    PackedData,  SCALAR, ORDINAL ) \
       TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsMatrix, NonPackedData,  SCALAR, ORDINAL )
 
 #define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
          UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, int) \
          UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, float)
      UNIT_TEST_GROUP_ORDINAL(int)
      typedef short int ShortInt; UNIT_TEST_GROUP_ORDINAL(ShortInt)

}
