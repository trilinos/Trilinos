#ifndef PACKAGES_STK_STK_LEARNING_KOKKOS_NGPATOMICS_H_
#define PACKAGES_STK_STK_LEARNING_KOKKOS_NGPATOMICS_H_

#include <stk_util/stk_config.h>
#include <Kokkos_Core.hpp>

namespace ngp {

template <typename T> STK_FUNCTION
void atomic_add(T *dest, const T src)
{
#if defined(KOKKOS_HAVE_CUDA) || defined(KOKKOS_HAVE_OPENMP)
    Kokkos::atomic_add(dest, src);
#else
    *dest += src;
#endif
}

}

#endif /* PACKAGES_STK_STK_LEARNING_KOKKOS_NGPATOMICS_H_ */
