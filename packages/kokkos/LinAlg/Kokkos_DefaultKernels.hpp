#ifndef KOKKOS_DEFAULT_KERNELS_
#define KOKKOS_DEFAULT_KERNELS_

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_DefaultSparseOps.hpp"
#include "Kokkos_DefaultBlockSparseOps.hpp"

namespace Kokkos {

  template <class Scalar, class Ordinal, class Node>
  struct DefaultKernels {
    typedef Kokkos::DefaultSparseOps     <Scalar,Ordinal,Node>      SparseOps;
    typedef Kokkos::DefaultBlockSparseOps<Scalar,Ordinal,Node> BlockSparseOps;
  };

}


#endif
