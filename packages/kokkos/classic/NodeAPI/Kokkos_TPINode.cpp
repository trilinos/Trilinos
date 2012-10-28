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

#include "Kokkos_TPINode.hpp"
#include <iostream>
#include <Teuchos_Assert.hpp>

namespace Kokkos {

  TPINode::TPINode(ParameterList &plist) 
  {
    using std::cout;
    using std::cerr;
    using std::endl;

    ParameterList params = getDefaultParameters();
    params.setParameters(plist);
    curNumThreads_ = params.get<int>("Num Threads");
    int verbose = params.get<int>("Verbose");
    TEUCHOS_TEST_FOR_EXCEPTION(curNumThreads_ < 0, std::runtime_error, 
        "TPINode::TPINode(): invalid ""Num Threads"" specification.");
    if (verbose) {
      cout << "TPINode initializing with numThreads == " << curNumThreads_ << std::endl;
    }
    init(curNumThreads_);
  }

  ParameterList TPINode::getDefaultParameters() 
  {
    ParameterList params;
    params.set("Num Threads", 0);
    params.set("Verbose",     0);
    return params;
  }

  void TPINode::init(int numThreads) {
    if (curNumThreads_ >= 1) {
      TPI_Finalize();
    }
    curNumThreads_ = numThreads;
    if (curNumThreads_ >= 1) {
      TPI_Init(curNumThreads_);
    }
  }

  TPINode::~TPINode()
  {
    if (curNumThreads_ >= 1) {
      TPI_Finalize();
    }
  }

}
