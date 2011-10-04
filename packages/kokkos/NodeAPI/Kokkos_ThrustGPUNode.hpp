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

#ifndef KOKKOS_THRUSTGPUNODE_HPP_
#define KOKKOS_THRUSTGPUNODE_HPP_

#include "Kokkos_NodeAPIConfigDefs.hpp"
#include "Kokkos_CUDANodeMemoryModel.hpp"

namespace Kokkos {

/** \brief %Kokkos node interface to the Thrust library for NVIDIA CUDA-capable GPUs.
    \ingroup kokkos_node_api
    \ingroup kokkos_cuda_support
 */
class ThrustGPUNode : public CUDANodeMemoryModel {
  public:

      /*! \brief Constructor acceptings a list of parameters
          
          This constructor accepts the parameters:
          \param "Device Number" [int] Specifies the CUDA device to which the node will attach.
          \param "Verbose"       [int] Non-zero parameter specifies that the constructor is verbose, printing information about the the attached device. Default: 0.

          The constructor throw std::runtime_error if "Device Number" is outside the range \f$[0,numDevices)\f$, where \c numDevices is the number of CUDA devices reported by cudaGetDeviceCount().
          
       */
    ThrustGPUNode(Teuchos::ParameterList &pl);

    //! \brief Destructor has no effect.
    ~ThrustGPUNode();

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

} // namespace Kokkos

#endif
