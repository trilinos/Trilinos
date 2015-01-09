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

#include "Kokkos_TBBNode.hpp"

// mfh 08 Jan 2015: Don't enable the contents of this file unless the
// appropriate CMake option is enabled.  This avoids deprecation
// warnings once we deprecate this Node type.
#ifdef HAVE_TPETRACLASSIC_TBB

#include <Teuchos_ParameterList.hpp>
#include <iostream>

// tbb::task_scheduler_init KokkosClassic::TBBNode::tsi_(tbb::task_scheduler_init::deferred);

namespace KokkosClassic {

  TBBNode::TBBNode () :
    tsi_ (tbb::task_scheduler_init::deferred),
    alreadyInit_ (false)
  {
    Teuchos::ParameterList params = getDefaultParameters ();
    const int numThreads = params.get<int> ("Num Threads");
    //
    // Allow either int or bool value of "Verbose" parameter.
    //
    int verbose = 0; // default value
    bool gotVerbose = false;
    try {
      verbose = params.get<int> ("Verbose");
      gotVerbose = true;
    }
    catch (...) {}
    if (! gotVerbose) {
      try {
        const bool boolVerbose = params.get<bool> ("Verbose");
        verbose = boolVerbose ? 1 : 0;
        // gotVerbose = true; // no longer needed
      }
      catch (...) {}
    }

    if (numThreads >= 0) {
      if (verbose) {
        std::cout << "TBBNode initializing with numThreads == " << numThreads
                  << std::endl;
      }
      init (numThreads);
    }
  }


  TBBNode::TBBNode (Teuchos::ParameterList &pl) :
    tsi_ (tbb::task_scheduler_init::deferred),
    alreadyInit_ (false)
  {
    Teuchos::ParameterList params = getDefaultParameters ();
    params.setParameters (pl);
    const int numThreads = params.get<int> ("Num Threads");
    //
    // Allow either int or bool value of "Verbose" parameter.
    //
    int verbose = 0; // default value
    bool gotVerbose = false;
    try {
      verbose = params.get<int> ("Verbose");
      gotVerbose = true;
    }
    catch (...) {}
    if (! gotVerbose) {
      try {
        const bool boolVerbose = params.get<bool> ("Verbose");
        verbose = boolVerbose ? 1 : 0;
        // gotVerbose = true; // no longer needed
      }
      catch (...) {}
    }

    if (numThreads >= 0) {
      if (verbose) {
        std::cout << "TBBNode initializing with numThreads == " << numThreads
                  << std::endl;
      }
      init (numThreads);
    }
  }

  Teuchos::ParameterList TBBNode::getDefaultParameters () {
    Teuchos::ParameterList params;
    params.set ("Verbose"    , 0);
    params.set ("Num Threads",-1);
    return params;
  }

  void TBBNode::init (const int numThreads) {
    if (alreadyInit_) {
      tsi_.terminate ();
    }
    //
    if (numThreads >= 1) {
      tsi_.initialize (numThreads);
    }
    else {
      tsi_.initialize (tbb::task_scheduler_init::automatic);
    }
  }

  TBBNode::~TBBNode () {}

  std::string TBBNode::name () {
    return "TBB";
  }

} // namespace KokkosClassic

#endif // HAVE_TPETRACLASSIC_TBB
