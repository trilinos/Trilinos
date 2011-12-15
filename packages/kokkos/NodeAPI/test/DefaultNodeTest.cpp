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

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_Tuple.hpp>

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_DefaultNode.hpp"
#include "Kokkos_NodeHelpers.hpp"
#include "TestOps.hpp"

namespace {

  using Teuchos::Time;
  using Teuchos::TimeMonitor;
  using Kokkos::DefaultNode;
  using Kokkos::ReadyBufferHelper;
  using Teuchos::tuple;
  using Teuchos::ArrayRCP;
  using Teuchos::Tuple;
  using Teuchos::RCP;
  using Teuchos::rcp;

  int N = 1000;

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
    clp.setOption("test-size",&N,"Vector length for tests.");
  }

  //
  // UNIT TESTS
  // 

  ////
  TEUCHOS_UNIT_TEST( NodeAPI, DefaultNodeTest )
  {
    Time tAlloc("Alloc Time"), tInit("Init Op"), tSum("Sum Op"), tFree("Free Time");
    typedef DefaultNode::DefaultNodeType NODE;
    typedef ArrayRCP<const char>  cbuf;
    typedef ArrayRCP<      char> ncbuf;
    Teuchos::ArrayRCP<int> x;
    RCP<NODE> node = DefaultNode::getDefaultNode();
    ReadyBufferHelper<NODE> rbh(node);
    out << "Default Node Type: " << Teuchos::typeName(*node) << std::endl;
    int result;
    {
      TimeMonitor localTimer(tAlloc);
      x = node->allocBuffer<int>(N);
    }
    // set x[i] = 1, i=0:N-1
    {
      TimeMonitor localTimer(tInit);
      InitOp<int> wdp;
      rbh.begin();
      wdp.x = rbh.addNonConstBuffer(x);
      rbh.end();
      node->parallel_for(0,N,wdp);
    }
    // compute sum x[i], i=0:N-1
    {
      TimeMonitor localTimer(tSum);
      SumOp<int> wdp;
      rbh.begin();
      wdp.x = rbh.addConstBuffer<int>(x);
      rbh.end();
      result = node->parallel_reduce(0,N,wdp);
    }
    int expectedResult = (int)(N);
    TEST_EQUALITY(result, expectedResult);
    {
      TimeMonitor localTimer(tFree);
      x = Teuchos::null;
    }
    out << "allocBuffer Time: " << tAlloc.totalElapsedTime() << std::endl;
    out << "InitOp Time: " << tInit.totalElapsedTime() << std::endl;
    out << "SumOp Time: " << tSum.totalElapsedTime() << std::endl;
    out << "freeBuffer Time: " << tFree.totalElapsedTime() << std::endl;
  }

}
