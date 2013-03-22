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

#include <Teuchos_ParameterList.hpp>
#include "Kokkos_OpenMPNode.hpp"

namespace Kokkos {

  OpenMPNode::OpenMPNode () :
    curNumThreads_ (-1), // Default: Let OpenMP pick the number of threads
    verbose_ (false)     // Default: No verbose status output
  {
    init (curNumThreads_);
  }


  OpenMPNode::OpenMPNode (Teuchos::ParameterList &pl) :
    curNumThreads_ (-1), // Default: Let OpenMP pick the number of threads
    verbose_ (false)     // Default: No verbose status output
  {
    ParameterList params = getDefaultParameters();
    params.setParameters(pl);
    const int curNumThreads = params.get<int>("Num Threads");
    int verboseInt = params.get<int>("Verbose");
    bool verbose = (verboseInt != 0);
    if (verbose) {
      std::cout << "OpenMPNode initializing with \"Num Threads\" = "
                << curNumThreads << std::endl;
    }
    init (curNumThreads);
    curNumThreads_ = curNumThreads; // Now it's safe to set state.
    verbose_ = verbose;
  }

  ParameterList OpenMPNode::getDefaultParameters() 
  {
    ParameterList params;
    params.set("Verbose",      0);
    params.set("Num Threads", -1);
    return params;
  }

  OpenMPNode::~OpenMPNode() {}

  void OpenMPNode::init (int numThreads) {
    // mfh 10 Jul 2012: Don't set the number of threads if it's
    // negative (the default value of the "Num Threads" parameter is
    // -1) or 0.
    if (numThreads > 0) {
      omp_set_num_threads (numThreads);
    }
  }

}
