#ifndef FAKEKOKKOS_DEFAULT_KERNELS_
#define FAKEKOKKOS_DEFAULT_KERNELS_

namespace Kokkos {

  template <class Scalar, class Ordinal, class Node>
  struct DefaultKernels {
    typedef void SparseOps;
  };

}

#endif
