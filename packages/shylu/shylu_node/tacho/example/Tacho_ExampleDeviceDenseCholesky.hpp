#ifndef __TACHO_EXAMPLE_DEVICE_DENSE_CHOLESKY_HPP__
#define __TACHO_EXAMPLE_DEVICE_DENSE_CHOLESKY_HPP__

#include "Kokkos_Random.hpp"

#include "Tacho_Util.hpp"
#include "Tacho_Blas_External.hpp"
#include "Tacho_Lapack_External.hpp"

#include "Tacho_Chol.hpp"
#include "Tacho_Chol_OnDevice.hpp"

#include "Tacho_Trsv.hpp"
#include "Tacho_Trsv_OnDevice.hpp"

#include "Tacho_Gemv.hpp"
#include "Tacho_Gemv_OnDevice.hpp"

template<typename value_type>
int driver_chol (const int m, const bool verbose) {
  int max_iter = 1;

  Kokkos::Timer timer;
  const bool detail = false;

  typedef typename Tacho::UseThisDevice<Kokkos::DefaultExecutionSpace>::type device_type;
  typedef typename Tacho::UseThisDevice<Kokkos::DefaultHostExecutionSpace>::type host_device_type;

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
#else
    int handle_blas, handle_lapack; // dummy
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
      if (m < 20) {
        printf("AA  = \n");
        for (int i=0;i<m;++i) {
          for (int j=0;j<m;++j)
            std::cout << AA(i,j) << " ";
          printf("\n");
        }
      }
    }

    Kokkos::View<value_type**,Kokkos::LayoutLeft,device_type> x("x", m, 1);
    {
      Kokkos::parallel_for
        (Kokkos::RangePolicy<typename device_type::execution_space>(0,m),
         KOKKOS_LAMBDA(const int i) {
          x(i,0) = i+1;
        });
      Kokkos::fence();
    }

    Kokkos::View<value_type**,Kokkos::LayoutLeft,device_type> b("b", m, 1);
    {
      Tacho::Gemv<Tacho::Trans::NoTranspose,Tacho::Algo::OnDevice>
        ::invoke(handle_blas,
                 one,
                 A, x, 
                 zero,
                 b);
      Kokkos::fence();
    }
    
    /// factorizeCholesky
    for (int iter=0;iter<max_iter;++iter) {
      Kokkos::deep_copy(A, Aback);
      Kokkos::fence();

      Kokkos::View<int,device_type> dev("dev");
      
      timer.reset();
      int lwork(0);
#if defined (KOKKOS_ENABLE_CUDA) 
      Tacho::Lapack<value_type>::potrf_buffersize(handle_lapack, 
                                                  CUBLAS_FILL_MODE_UPPER,
                                                  m, 
                                                  A.data(), m, 
                                                  &lwork);
      printf("Cholesky lwork %d\n", lwork);
#endif
      Kokkos::View<value_type*,Kokkos::LayoutLeft,device_type> W("W", lwork);
      Tacho::Chol<Tacho::Uplo::Upper,Tacho::Algo::OnDevice>::invoke(handle_lapack, A, W);
      Kokkos::fence();
      {
        const double t = timer.seconds();
        printf("Cholesky factorization time %e\n", t);
      }
      {
        const auto dev_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), dev);
        if (dev_host()) printf("Cholesky returns non-zero dev info %d\n", dev_host());
      }
      
      Kokkos::deep_copy(x, b);
      Kokkos::fence();

      timer.reset();
      {
        Tacho::Trsv<Tacho::Uplo::Upper,Tacho::Trans::ConjTranspose,Tacho::Algo::OnDevice>
          ::invoke(handle_blas, Tacho::Diag::NonUnit(), A, x);
        Kokkos::fence();
        Tacho::Trsv<Tacho::Uplo::Upper,Tacho::Trans::NoTranspose,Tacho::Algo::OnDevice>
          ::invoke(handle_blas, Tacho::Diag::NonUnit(), A, x);
        Kokkos::fence();
      }
      {
        const double t = timer.seconds();
        printf("Cholesky solve time %e\n", t);
      }
      
      {
        const auto x_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), x);
        if (m < 20) {
          printf("x = \n");
          for (int i=0;i<m;++i)
            std::cout << x(i,0) << std::endl;
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
