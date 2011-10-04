// @HEADER
// ***********************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//                 Copyright (2008) Sandia Corporation
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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "Kokkos_TBBNode.hpp"
#include <Teuchos_ParameterList.hpp>
#include <iostream>

// tbb::task_scheduler_init Kokkos::TBBNode::tsi_(tbb::task_scheduler_init::deferred);

namespace Kokkos {

  TBBNode::TBBNode(Teuchos::ParameterList &pl) : alreadyInit_(false), tsi_(tbb::task_scheduler_init::deferred) {
    int numThreads = pl.get<int>("Num Threads",-1);
    int verbose = pl.get<int>("Verbose",0);
    if (numThreads >= 0) {
      if (verbose) {
        std::cout << "TBBNode initializing with numThreads == " << numThreads << std::endl;
      }
      init(numThreads);
    }
  }

  void TBBNode::init(int numThreads) {
    if (alreadyInit_) {
      tsi_.terminate();
    }
    // 
    if (numThreads >= 1) {
      tsi_.initialize(numThreads);
    }
    else {
      tsi_.initialize(tbb::task_scheduler_init::automatic);
    }
  }

  TBBNode::~TBBNode() {}

}
