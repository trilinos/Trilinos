#ifndef __TACHO_EXAMPLE_DEVICE_DENSE_CHOLESKY_HPP__
#define __TACHO_EXAMPLE_DEVICE_DENSE_CHOLESKY_HPP__

#include "Kokkos_Random.hpp"

#include "Tacho_Util.hpp"
#include "Tacho_Blas_External.hpp"
#include "Tacho_Lapack_External.hpp"

template<typename value_type>
int driver_chol (const int m, const bool verbose) {
  int max_iter = 1;

  Kokkos::Impl::Timer timer;
  const bool detail = false;

  typedef typename Tacho::UseThisDevice<Kokkos::DefaultExecutionSpace>::device_type device_type;
  typedef typename Tacho::UseThisDevice<Kokkos::DefaultHostExecutionSpace>::device_type host_device_type;

  Tacho::printExecSpaceConfiguration<typename device_type::execution_space>("DeviceSpace", detail);
  Tacho::printExecSpaceConfiguration<typename host_device_type::execution_space>("HostSpace",   detail);
  printf("\n\n");

#if defined (KOKKOS_ENABLE_CUDA) 
  printf("CUDA testing\n");
#else
  printf("Host testing\n");
#endif
  int r_val = 0;  
  {
    const value_type one(1), zero(0);
    
#if defined (KOKKOS_ENABLE_CUDA) 
    cublasHandle_t handle_blas;
    cusolverDnHandle_t handle_lapack;
    {
      {
        const int status = cublasCreate(&handle_blas); 
        if (status) printf("Nonzero error from cublasCreate %d\n", status);
      }
      {
        const int status = cusolverDnCreate(&handle_lapack); 
        if (status) printf("Nonzero error from cusolverDnCreate %d\n", status);
      }
    }
#endif

    Kokkos::View<value_type**,Kokkos::LayoutLeft,device_type> Arand("Arand", m, m);
    {
      Kokkos::Random_XorShift64_Pool<device_type> random(13718);
      Kokkos::fill_random(Arand, random, one);
    }

    Kokkos::View<value_type**,Kokkos::LayoutLeft,device_type> A("A", m, m);
    {
      const value_type * Arand_ptr = Arand.data();
      value_type * A_ptr = A.data();
      
      timer.reset();
#if defined (KOKKOS_ENABLE_CUDA) 
      Tacho::Blas<value_type>::gemm(handle_blas, 
                                    CUBLAS_OP_N, CUBLAS_OP_T,
                                    m, m, m,
                                    one,
                                    Arand_ptr, m,
                                    Arand_ptr, m,
                                    zero,
                                    A_ptr, m);
#else
      Tacho::Blas<value_type>::gemm('N', 'T', 
                                    m, m, m,
                                    one,
                                    Arand_ptr, m, 
                                    Arand_ptr, m, 
                                    zero,
                                    A_ptr, m);
#endif
      Kokkos::fence();
      const double t = timer.seconds();
      printf("hermitianize time %e\n", t);
    } 

    Kokkos::View<value_type**,Kokkos::LayoutLeft,device_type> Aback("Aback", m, m);    
    {
      Kokkos::deep_copy(Aback, A);
      auto AA = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), Aback);
      if (0) {
        printf("AA  = \n");
        for (int i=0;i<m;++i) {
          for (int j=0;j<m;++j)
            std::cout << AA(i,j) << " ";
          printf("\n");
        }
      }
    }

    Kokkos::View<value_type*,Kokkos::LayoutLeft,device_type> x("x", m);
    {
      Kokkos::parallel_for
        (Kokkos::RangePolicy<typename device_type::execution_space>(0,m),
         KOKKOS_LAMBDA(const int i) {
          x(i) = i;
        });
    }

    Kokkos::View<value_type*,Kokkos::LayoutLeft,device_type> b("b", m);
    {
#if defined (KOKKOS_ENABLE_CUDA) 
      Tacho::Blas<value_type>::gemv(handle_blas,
                                    CUBLAS_OP_N,
                                    m, m,
                                    one,
                                    A.data(), m,
                                    x.data(), 1,
                                    zero,
                                    b.data(), 1);
#else
      Tacho::Blas<value_type>::gemv('N',
                                    m, m,
                                    one,
                                    A.data(), m,
                                    x.data(), 1,
                                    zero,
                                    b.data(), 1);
#endif
      Kokkos::fence();
    }
    
    /// factorizeCholesky
    for (int iter=0;iter<max_iter;++iter) {
      Kokkos::deep_copy(A, Aback);
      Kokkos::fence();

      value_type * A_ptr = A.data();
      value_type * x_ptr = x.data();
      Kokkos::View<int,device_type> dev("dev");
      
      timer.reset();
#if defined (KOKKOS_ENABLE_CUDA) 
      int lwork(0);
      Tacho::Lapack<value_type>::potrf_buffersize(handle_lapack, 
                                                  CUBLAS_FILL_MODE_UPPER,
                                                  m, 
                                                  A_ptr, m, 
                                                  &lwork);
      printf("Cholesky lwork %d\n", lwork);
      Kokkos::View<value_type*,Kokkos::LayoutLeft,device_type> W("W", lwork);
      value_type * W_ptr = W.data();
      
      Tacho::Lapack<value_type>::potrf(handle_lapack,
                                       CUBLAS_FILL_MODE_UPPER,
                                       m, 
                                       A_ptr, m,
                                       W_ptr, lwork,
                                       dev.data());
#else
      Tacho::Lapack<value_type>::potrf('U',
                                       m,
                                       A_ptr, m,
                                       dev.data());
#endif
      Kokkos::fence();
      const double t = timer.seconds();
      printf("Cholesky time %e\n", t);
      {
        const auto dev_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), dev);
        if (dev_host()) printf("Cholesky returns non-zero dev info %d\n", dev_host());
      }

      Kokkos::deep_copy(x, b);
      Kokkos::fence();

#if defined (KOKKOS_ENABLE_CUDA) 
      Tacho::Blas<value_type>::trsv(handle_blas,
                                    CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_C, CUBLAS_DIAG_NON_UNIT,
                                    m,
                                    A_ptr, m,
                                    x, 1);
      Kokkos::fence();
      Tacho::Blas<value_type>::trsv(handle_blas,
                                    CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT,
                                    m,
                                    A_ptr, m,
                                    x, 1);      
#else
      Tacho::Blas<value_type>::trsv('U', 'C', 'N',
                                    m,
                                    A_ptr, m,
                                    x_ptr, 1);
      Kokkos::fence();
      Tacho::Blas<value_type>::trsv('U', 'N', 'N',
                                    m,
                                    A_ptr, m,
                                    x_ptr, 1);            
#endif
      {
        const auto x_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), x);
        if (0) {
          printf("x = \n");
          for (int i=0;i<m;++i)
            std::cout << x(i) << std::endl;
        }
      }
    }

#if defined (KOKKOS_ENABLE_CUDA) 
    {
      {
        const int status = cublasDestroy(handle_blas);
        if (status) printf("Nonzero error from cublasDestroy %d\n", status);
      }
      {
        const int status = cusolverDnDestroy(handle_lapack);
        if (status) printf("Nonzero error from cusolverDnDestroy %d\n", status);
      }
    }
#endif
  }


  return r_val;
}

#endif
