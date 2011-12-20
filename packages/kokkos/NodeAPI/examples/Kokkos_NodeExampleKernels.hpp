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

#ifndef KOKKOS_NODE_EXAMPLE_KERNELS
#define KOKKOS_NODE_EXAMPLE_KERNELS

#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Kokkos_NodeHelpers.hpp>

/** \file Kokkos_NodeExampleKernels.hpp
    \brief A file containing example kernels for parallel_for and parallel_reduce.
 */

#ifndef KERNEL_PREFIX
#define KERNEL_PREFIX
#endif

namespace KokkosExamples {

  /** \brief A kernel for parallel_for, with a simple vector initialization.
    */
  struct VecInit {
    //! 
    int * x;
    //! 
    KERNEL_PREFIX inline void execute(int i) {x[i] = i;}
  };

  /** \brief A wrapper for kernel VecInit
   */
  template <class NODE>
  void initVec(Teuchos::RCP<NODE> node, Teuchos::ArrayRCP<int> data) {
    Kokkos::ReadyBufferHelper<NODE> rbh(node);
    VecInit init;
    // ready the buffer and encapsulate the kernel arguments
    rbh.begin();
    init.x = rbh.addNonConstBuffer(data);
    rbh.end();  // this call node->readyBuffers()
    // execute the kernel
    node->parallel_for(0,data.size(),init);
  }

  /** \brief A kernel for parallel_reduce, with a simple sum reduction.
   */
  struct VecReduce {
    //!
    const int * x;
    //!
    typedef int ReductionType;
    //!
    KERNEL_PREFIX static inline int identity()            { return 0;    }
    //!
    KERNEL_PREFIX        inline int generate(int i)       { return x[i]; }
    //!
    KERNEL_PREFIX        inline int reduce  (int a, int b){ return a+b;  }
  };

  /** \brief A wrapper for for VecReduce 
    */
  template <class NODE>
  int reduceVec(Teuchos::RCP<NODE> node, Teuchos::ArrayRCP<const int> data) {
    Kokkos::ReadyBufferHelper<NODE> rbh(node);
    VecReduce reduce;
    // ready the buffer and encapsulate the kernel arguments
    rbh.begin();
    reduce.x = rbh.addConstBuffer(data);
    rbh.end();  // this call node->readyBuffers()
    int ret = node->parallel_reduce(0,data.size(),reduce);
    return ret;
  }

  /** \example SimpleNodeExample.cpp
    * This is an example of the simple kernels KokkosExamples::VecInit and KokkosExamples::VecReduce, as well as their wrappers KokkosExamples::reduceVec() and KokkosExamples::initVec().
    */

} // end of namespace KokkosExamples

#endif
