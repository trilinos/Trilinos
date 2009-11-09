#ifndef KOKKOS_CUDANODEUTILS_HPP_
#define KOKKOS_CUDANODEUTILS_HPP_

namespace Kokkos {

  class CUDANodeDeallocator {
    public:
      static void free(void *ptr);
  };

}

#endif // KOKKOS_CUDANODEUTILS_HPP_
