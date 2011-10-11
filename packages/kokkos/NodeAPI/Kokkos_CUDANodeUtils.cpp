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

#include "Kokkos_CUDANodeUtils.hpp"
#include <iostream>
#include <cuda_runtime.h>

namespace Kokkos {

  CUDANodeDeallocator::CUDANodeDeallocator(size_t sizeInBytes, const Teuchos::RCP<CUDANodeMemoryModel> &node)
#ifdef HAVE_KOKKOS_CUDA_NODE_MEMORY_PROFILING
  : node_(node)
  , allocSize_(sizeInBytes)
#endif
  {
    (void)sizeInBytes;
    (void)node;
  }

  void CUDANodeDeallocator::free(void *ptr) {
    cudaError_t err = cudaFree(ptr);
    TEST_FOR_EXCEPTION( cudaSuccess != err, std::runtime_error, 
        "Kokkos::CUDANodeDeallocator::free(): cudaFree() returned error:\n"
        << cudaGetErrorString(err) 
      );
#ifdef HAVE_KOKKOS_CUDA_NODE_MEMORY_PROFILING
    node_->allocSize_ -= allocSize_;
#endif
  }

}
