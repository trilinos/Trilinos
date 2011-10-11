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

#include "Kokkos_TPINode.hpp"
#include <iostream>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_TestForException.hpp>

namespace Kokkos {

  TPINode::TPINode(Teuchos::ParameterList &plist) {
    using std::cout;
    using std::cerr;
    using std::endl;

    curNumThreads_ = plist.get<int>("Num Threads", 0);
    int verbose = plist.get<int>("Verbose",0);
    TEST_FOR_EXCEPTION(curNumThreads_ < 0, std::runtime_error, 
        "TPINode::TPINode(): invalid ""Num Threads"" specification.");
    if (verbose) {
      cout << "TPINode initializing with numThreads == " << curNumThreads_ << std::endl;
    }
    init(curNumThreads_);
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
