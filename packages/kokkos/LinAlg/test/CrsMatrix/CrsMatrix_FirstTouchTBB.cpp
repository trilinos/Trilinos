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

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_DefaultArithmetic.hpp"
#include "Kokkos_DefaultKernels.hpp"
#include "Kokkos_CrsGraph.hpp"
#include "Kokkos_CrsMatrix.hpp"
#include "Kokkos_Version.hpp"

#include "Kokkos_TBBNode.hpp"

namespace {

  using Kokkos::MultiVector;
  using Kokkos::CrsMatrixHostCompute;
  using Kokkos::CrsGraphHostCompute;
  using Kokkos::FirstTouchHostCrsMatrix;
  using Kokkos::FirstTouchHostCrsGraph;
  using Kokkos::DefaultArithmetic;
  using Kokkos::SerialNode;
  using Kokkos::DefaultKernels;
  using Teuchos::ArrayRCP;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::null;
  using std::endl;

#include "CrsMatrix_CrsTimingTests.hpp"

  using Kokkos::TBBNode;
  RCP<TBBNode> tbbnode;

  template <>
  RCP<TBBNode> getNode<TBBNode>() {
    if (tbbnode == null) {
      Teuchos::ParameterList pl;
      pl.set<int>("Num Threads",Test::numThreads);
      tbbnode = rcp(new TBBNode(pl));
    }
    return tbbnode;
  }

  //
  // Test
  // 

  typedef TBBNode                                             Node;
  // mfh 23 Feb 2011:
  // You don't need "typename" here, since none of the three 
  // parameters of DefaultKernels are template parameters.
  //typedef typename DefaultKernels<void,int,Node>::SparseOps   DSM;
  typedef DefaultKernels<void,int,Node>::SparseOps            DSM;
  typedef CrsMatrixHostCompute<double,int,Node,DSM>           StandardMat;
  typedef CrsGraphHostCompute<int,Node,DSM>                   StandardGraph;
  typedef FirstTouchHostCrsMatrix<double,int,Node,DSM>        FirstTouchMat;
  typedef FirstTouchHostCrsGraph<int,Node,DSM>                FirstTouchGraph;

  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsTiming, PowerTriDiag, StandardGraph,   StandardMat   )
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsTiming, PowerTriDiag, FirstTouchGraph, FirstTouchMat )

}

