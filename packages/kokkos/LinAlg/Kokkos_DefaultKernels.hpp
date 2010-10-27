#ifndef KOKKOS_DEFAULT_KERNELS_
#define KOKKOS_DEFAULT_KERNELS_

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_DefaultSparseOps.hpp"
#include "Kokkos_DefaultBlockSparseOps.hpp"
#include "Kokkos_DefaultRelaxation.hpp"
#if HAVE_KOKKOS_CUSP
#include "Kokkos_CUSPSparseOps.hpp"
#endif

namespace Kokkos {

  /** \brief Traits class providing default kernel types for CRS, block CRS and relaxation kernels.
      \ingroup kokkos_crs_ops
   */
  template <class Scalar, class Ordinal, class Node>
  struct DefaultKernels {
    typedef DefaultHostSparseOps <void  ,Ordinal,Node>  SparseOps;
    typedef DefaultBlockSparseOps<Scalar,Ordinal,Node>  BlockSparseOps;
    typedef DefaultRelaxation    <Scalar,Ordinal,Node>  Relaxations;
  };

  /** \brief Traits class providing default kernel types for CRS, block CRS and relaxation kernels.
      \ingroup kokkos_crs_ops
    
      For ThrustGPUNode, defaults are the same as in general, except that the default sparse ops should be provided by 
      DefaultDeviceSparseOps.
   */
#ifdef HAVE_KOKKOS_CUSP
  class ThrustGPUNode;
  template <class Scalar, class Ordinal>
  struct DefaultKernels<Scalar,Ordinal,ThrustGPUNode> {
    typedef CUSPSparseOps<void  ,Ordinal,ThrustGPUNode>           SparseOps;
    typedef DefaultBlockSparseOps <Scalar,Ordinal,ThrustGPUNode>  BlockSparseOps;
    typedef DefaultRelaxation     <Scalar,Ordinal,ThrustGPUNode>  Relaxations;
  };
#else
  class ThrustGPUNode;
  template <class Scalar, class Ordinal>
  struct DefaultKernels<Scalar,Ordinal,ThrustGPUNode> {
    typedef DefaultDeviceSparseOps<void  ,Ordinal,ThrustGPUNode>  SparseOps;
    typedef DefaultBlockSparseOps <Scalar,Ordinal,ThrustGPUNode>  BlockSparseOps;
    typedef DefaultRelaxation     <Scalar,Ordinal,ThrustGPUNode>  Relaxations;
  };
#endif

#ifdef HAVE_KOKKOS_TREAT_SERIALNODE_AS_DEVICE
  /** \brief Traits class providing default kernel types for CRS, block CRS and relaxation kernels.
      \ingroup kokkos_crs_ops
      
      If compiled with HAVE_KOKKOS_TREAT_SERIALNODE_AS_DEVICE, then the default sparse ops for SerialNode should be provided 
      by DefaultDeviceSparseOps. 
   */
  class SerialNode;
  template <class Scalar, class Ordinal>
  struct DefaultKernels<Scalar,Ordinal,SerialNode> {
    typedef DefaultDeviceSparseOps<void  ,Ordinal,SerialNode>     SparseOps;
    typedef DefaultBlockSparseOps <Scalar,Ordinal,ThrustGPUNode>  BlockSparseOps;
    typedef DefaultRelaxation     <Scalar,Ordinal,ThrustGPUNode>  Relaxations;
  };
#endif

}

#endif
