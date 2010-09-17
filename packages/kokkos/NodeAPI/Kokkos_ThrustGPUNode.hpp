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

    //@} 

  private:
    int totalMem_;
};

} // namespace Kokkos

#endif
