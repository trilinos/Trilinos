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
