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

#ifndef KOKKOS_THRUSTGPUNODE_HPP_
#define KOKKOS_THRUSTGPUNODE_HPP_

#include "Kokkos_NodeAPIConfigDefs.hpp"
#include "Kokkos_CUDANodeMemoryModel.hpp"
#include "Kokkos_ThrustGPUWrappers.hpp"

namespace Kokkos {

/** \brief %Kokkos node interface to the Thrust library for NVIDIA CUDA-capable GPUs.
    \ingroup kokkos_node_api
    \ingroup kokkos_cuda_support
 */
class ThrustGPUNode : public CUDANodeMemoryModel {
public:
  //! Constructor that sets default parameters.
  ThrustGPUNode ();

  /*! \brief Constructor that takes a list of parameters.
          
    We accept the following parameters:
    - "Device Number" [int] The CUDA device to which the Node will attach.
    - "Verbose" [int] Non-zero parameter specifies that the
      constructor is verbose, printing information about the the
      attached device. Default: 0.

    The constructor will throw std::runtime_error if "Device Number"
    is outside the range \f$[0,numDevices)\f$, where \c numDevices is
    the number of CUDA devices reported by cudaGetDeviceCount().
  */
  ThrustGPUNode(Teuchos::ParameterList &pl);

  //! \brief Destructor has no effect.
  ~ThrustGPUNode();

  /*! \brief Get default parameters for this node */
  static ParameterList getDefaultParameters();

  //@{ Computational methods

  //! \begin parallel for skeleton, a wrapper around thrust::for_each. See \ref kokkos_node_api "Kokkos Node API"
  template <class WDP>
  static void parallel_for(int begin, int end, WDP wdp);

  //! \begin parallel reduction skeleton, a wrapper around thrust::transform_reduce. See \ref kokkos_node_api "Kokkos Node API"
  template <class WDP>
  static typename WDP::ReductionType
  parallel_reduce(int begin, int end, WDP wd);

  //! \begin Block until all node work is complete. Aids in accurate timing of multiple kernels.
  void sync() const;

  //@} 

private:
  int totalMem_;
};

template <class WDP>
void ThrustGPUNode::parallel_for(int begin, int end, WDP wd) {
  ThrustGPUNodeDetails::parallel_for(begin,end,wd);
}

template <class WDP>
typename WDP::ReductionType
ThrustGPUNode::parallel_reduce(int begin, int end, WDP wd) 
{
  return ThrustGPUNodeDetails::parallel_reduce(begin,end,wd);
}

} // namespace Kokkos

#endif
