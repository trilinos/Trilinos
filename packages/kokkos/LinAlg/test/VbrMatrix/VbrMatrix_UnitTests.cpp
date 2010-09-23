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
#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayRCP.hpp>

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_DefaultNode.hpp"
#include "Kokkos_VbrMatrix.hpp"
#include "Kokkos_Version.hpp"

namespace {

  using Kokkos::DefaultNode;
  using Kokkos::VbrMatrix;
  using Teuchos::ArrayRCP;
  using Teuchos::Array;
  using Teuchos::RCP;
  using Teuchos::arcpFromArray;

  typedef Kokkos::DefaultNode::DefaultNodeType Node;

  int N;
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
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( VbrMatrix, PackedData, Scalar, Ordinal )
  {
    const size_t N = 10;
    RCP<Node> node = Kokkos::DefaultNode::getDefaultNode();
    // test non-empty
    {
      VbrMatrix<Scalar,Ordinal,Node> A(N,node);
      TEST_EQUALITY( A.getNumBlockRows(), N );
      TEST_EQUALITY_CONST( A.isPacked(), false );
      TEST_EQUALITY_CONST( A.isEmpty(),  true );
      //
      ArrayRCP<Scalar>  vals   = node->template allocBuffer<Scalar>(2*N);
      ArrayRCP<Ordinal> rptr   = node->template allocBuffer<Ordinal>(2*N);
      ArrayRCP<Ordinal> cptr   = node->template allocBuffer<Ordinal>(2*N);
      ArrayRCP<size_t> bptr   = node->template allocBuffer<size_t>(2*N);
      ArrayRCP<Ordinal> bindx   = node->template allocBuffer<Ordinal>(2*N);
      ArrayRCP<Ordinal> indx   = node->template allocBuffer<Ordinal>(2*N);
      A.setPackedValues(vals,rptr,cptr,bptr,bindx,indx);
      TEST_EQUALITY_CONST( A.isPacked(), true );
      TEST_EQUALITY_CONST( A.isEmpty(), false );
      //
      A.clear();
      TEST_EQUALITY( A.getNumBlockRows(), N );
      TEST_EQUALITY_CONST( A.isPacked(), false );
      TEST_EQUALITY_CONST( A.isEmpty(),  true );
    }
    // test empty
    {
      VbrMatrix<Scalar,Ordinal,Node> A(N,node);
      TEST_EQUALITY( A.getNumBlockRows(), N );
      TEST_EQUALITY_CONST( A.isPacked(), false );
      TEST_EQUALITY_CONST( A.isEmpty(),  true );
      //
      ArrayRCP<Scalar>  vals;
      ArrayRCP<Ordinal> rptr;
      ArrayRCP<Ordinal> cptr;
      ArrayRCP<size_t> bptr;
      ArrayRCP<Ordinal> bindx;
      ArrayRCP<Ordinal> indx;
      A.setPackedValues(vals,rptr,cptr,bptr,bindx,indx);
      TEST_EQUALITY_CONST( A.isPacked(), false );
      TEST_EQUALITY_CONST( A.isEmpty(),  true );
      //
      A.clear();
      TEST_EQUALITY( A.getNumBlockRows(), N );
      TEST_EQUALITY_CONST( A.isPacked(), false );
      TEST_EQUALITY_CONST( A.isEmpty(),  true );
    }
  }

 #define UNIT_TEST_GROUP_ORDINAL_SCALAR( ORDINAL, SCALAR ) \
       TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( VbrMatrix,    PackedData,  SCALAR, ORDINAL )
 
 #define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
          UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, int) \
          UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, float)
      UNIT_TEST_GROUP_ORDINAL(int)
      typedef short int ShortInt; UNIT_TEST_GROUP_ORDINAL(ShortInt)

}
