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

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_DefaultNode.hpp"
#include "Kokkos_MultiVector.hpp"
#include "Kokkos_Version.hpp"

namespace {

  using Kokkos::DefaultNode;
  using Kokkos::MultiVector;
  using Teuchos::ArrayRCP;
  using Teuchos::RCP;

  size_t N = 1000;

  typedef Kokkos::DefaultNode::DefaultNodeType Node;

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
    int n = N;
    clp.setOption("test-size",&n,"Vector length for tests.");
    N = n;
  }

  //
  // UNIT TESTS
  // 

  // check that default constructor zeros out, for both V and MV
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, DefaultConstructor, Scalar, Ordinal )
  {
    RCP<Node> node = Kokkos::DefaultNode::getDefaultNode();
    MultiVector<Scalar,Node> A(node);
    TEST_EQUALITY_CONST(A.getNumRows(), 0);
    TEST_EQUALITY_CONST(A.getNumCols(), 0);
    TEST_EQUALITY_CONST(A.getStride(), 0);
  }

  // check copy constructor
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, CopyConstructor, Scalar, Ordinal )
  {
    RCP<Node> node = Kokkos::DefaultNode::getDefaultNode();
    MultiVector<Scalar,Node> A(node);
    ArrayRCP<Scalar> buf = A.getNode()->template allocBuffer<Scalar>(2*N);
    A.initializeValues(N,2,buf,N);
    {
      MultiVector<Scalar,Node> Acopy(A);
      TEST_EQUALITY_CONST(Acopy.getNumRows(), N);
      TEST_EQUALITY_CONST(Acopy.getNumCols(), 2);
      TEST_EQUALITY_CONST(Acopy.getStride(), N);
      TEST_EQUALITY(Acopy.getValues(0), buf);
      TEST_INEQUALITY(Acopy.getValues(1), buf);
    }
  }

  // check that non-default constructor honors given parameters
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, InitializeAndAccess, Scalar, Ordinal )
  {
    RCP<Node> node = Kokkos::DefaultNode::getDefaultNode();
    MultiVector<Scalar,Node> A(node);
    ArrayRCP<Scalar> buf = A.getNode()->template allocBuffer<Scalar>(2*N);
    A.initializeValues(N,2,buf,N);
    TEST_EQUALITY_CONST(A.getNumRows(), N);
    TEST_EQUALITY_CONST(A.getNumCols(), 2);
    TEST_EQUALITY_CONST(A.getStride(), N);
    TEST_EQUALITY(A.getValues(0), buf);
    TEST_INEQUALITY(A.getValues(1), buf);
    buf = Teuchos::null;
  }


#define UNIT_TEST_GROUP_ORDINAL_SCALAR( ORDINAL, SCALAR ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, DefaultConstructor, SCALAR, ORDINAL ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, CopyConstructor   , SCALAR, ORDINAL ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, InitializeAndAccess, SCALAR, ORDINAL )

#define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
  UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, int) \
  UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, float)
  UNIT_TEST_GROUP_ORDINAL(int)
    typedef short int ShortInt; UNIT_TEST_GROUP_ORDINAL(ShortInt)

}
