#ifndef KOKKOSBLAS_CUDA_TPL_HPP_
#define KOKKOSBLAS_CUDA_TPL_HPP_

#if defined (KOKKOSKERNELS_ENABLE_TPL_CUBLAS)
#include<KokkosBlas_tpl_spec.hpp>

namespace KokkosBlas {
namespace Impl {

CudaBlasSingleton::CudaBlasSingleton()
{
    cublasStatus_t stat = cublasCreate(&handle);
    if (stat != CUBLAS_STATUS_SUCCESS)
        Kokkos::abort("CUBLAS initialization failed\n");

    Kokkos::push_finalize_hook ([&] () { 
        cublasDestroy(handle);
    });
}

CudaBlasSingleton & CudaBlasSingleton::singleton()
{ static CudaBlasSingleton s ; return s ; }

}
}
#endif

#if defined (KOKKOSKERNELS_ENABLE_TPL_MAGMA)
#include<KokkosBlas_tpl_spec.hpp>

namespace KokkosBlas {
namespace Impl {

MagmaSingleton::MagmaSingleton()
{
    magma_int_t stat = magma_init();
    if (stat != MAGMA_SUCCESS)
        Kokkos::abort("MAGMA initialization failed\n");

    Kokkos::push_finalize_hook ([&] () { 
        magma_finalize();
    });
}

MagmaSingleton & MagmaSingleton::singleton()
{ static MagmaSingleton s ; return s ; }

}
}
#endif

#endif