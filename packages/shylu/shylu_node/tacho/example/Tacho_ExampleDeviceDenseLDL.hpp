#ifndef __TACHO_EXAMPLE_DEVICE_DENSE_LDL_HPP__
#define __TACHO_EXAMPLE_DEVICE_DENSE_LDL_HPP__

#include "Kokkos_Random.hpp"

#include "Tacho_Util.hpp"
#include "Tacho_Blas_External.hpp"
#include "Tacho_Lapack_External.hpp"

#include "Tacho_LDL.hpp"
#include "Tacho_LDL_External.hpp"
#include "Tacho_LDL_OnDevice.hpp"

#include "Tacho_Gemv.hpp"
#include "Tacho_Gemv_OnDevice.hpp"

template<typename value_type>
int driver_ldl (const int m, const bool verbose) {
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
    Kokkos::Cuda exec_instance;
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
    int handle_blas(0), handle_lapack(0), exec_instance(0); // dummy
#endif

    Kokkos::View<value_type**,Kokkos::LayoutLeft,device_type> A("A", m, m);
    {
      Kokkos::Random_XorShift64_Pool<device_type> random(13718);
      Kokkos::fill_random(A, random, one);
      Kokkos::fence();
     
      Kokkos::parallel_for
        (Kokkos::RangePolicy<typename device_type::execution_space>(0,m*m),
         KOKKOS_LAMBDA(const int ij) {
          const int i = ij%m, j = ij/m;
          A(i,j) = (i < j ? A(j,i) : A(i,j));
        });
    } 

    Kokkos::View<value_type**,Kokkos::LayoutLeft,device_type> Aback("Aback", m, m);    
    {
      Kokkos::deep_copy(Aback, A);
      if (m < 20) {
        auto AA = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), Aback);
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
          x(i, 0) = i+1;
        });
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

    /// factorize LDLt
    for (int iter=0;iter<max_iter;++iter) {
      Kokkos::deep_copy(A, Aback);
      Kokkos::fence();

      /// pivot needs 3m array size: lapack, flame, permutation
      Kokkos::View<value_type**,Kokkos::LayoutRight,device_type> D("D", m, 2);
      Kokkos::View<int*,Kokkos::LayoutLeft,device_type> p("pivots", 4*m);

      value_type * A_ptr = A.data();
      value_type * x_ptr = x.data();
      int * p_ptr = p.data();
      Kokkos::View<int,device_type> dev("dev");
        
      timer.reset();
      Kokkos::View<value_type*,Kokkos::LayoutLeft,device_type> W;
      const int lwork = Tacho::LDL<Tacho::Uplo::Lower,Tacho::Algo::OnDevice>
        ::invoke(handle_lapack, A, p, W);
      W = Kokkos::View<value_type*,Kokkos::LayoutLeft,device_type>("W", lwork);

      printf("LDLt lwork %d\n", lwork);
      Tacho::LDL<Tacho::Uplo::Lower,Tacho::Algo::OnDevice>
        ::invoke(handle_lapack, A, p, W);
      Kokkos::fence();
      
      Tacho::LDL<Tacho::Uplo::Lower,Tacho::Algo::OnDevice>
        ::modify(exec_instance, A, p, D);
      Kokkos::fence();

      const double t = timer.seconds();
      {
        const auto dev_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), dev);
        if (dev_host()) printf("LDL returns non-zero dev info %d\n", dev_host());
      }
        
      auto AA = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A);
      auto pp = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), p);
      {
        printf("LDL  = \n");
        for (int i=0;i<m;++i) {
          for (int j=0;j<m;++j)
            std::cout << AA(i,j) << " ";
          printf("\n");
        }
        printf("lapack piv = \n");
        for (int i=0;i<m;++i) {
          printf("%d\n", pp(i));
        }
        printf("flame piv = \n");
        for (int i=0;i<m;++i) {
          printf("%d\n", pp(i+m));
        }
        printf("permutation vector = \n");
        for (int i=0;i<m;++i) {
          printf("%d\n", pp(i+2*m));
        }
      }
      printf("LDLt time %e\n", t);

      Kokkos::fence();
        
      /// compute diagonals
      auto DD = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), D);
      {
        printf("D = \n");
        for (int i=0;i<m;++i) {
          printf("lapack %d, flame %d, perm %d, D %e %e\n", pp(i), pp(i+m), pp(i+2*m), DD(i,0), DD(i,1));
        }
        printf("L  = \n");
        for (int i=0;i<m;++i) {
          for (int j=0;j<m;++j)
            printf("%e ", AA(i,j));
          printf("\n");
        }
      }

      /// solve 
      {
        Kokkos::fence();
        auto perm = Kokkos::subview(p, Kokkos::pair<int,int>(2*m, 3*m));
        auto peri = Kokkos::subview(p, Kokkos::pair<int,int>(3*m, 4*m));
        /// copy and transpose
        Kokkos::parallel_for
          (Kokkos::RangePolicy<typename device_type::execution_space>(0,m),
           KOKKOS_LAMBDA(const int i) {           
            x(i,0) = b(perm(i),0);
          });
        {
          auto x_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), x);
          printf("x (permuted b) = \n");
          for (int i=0;i<m;++i) {
            printf("%e\n", x_host(i,0));
          }
        }
        
        /// trsv
#if defined (KOKKOS_ENABLE_CUDA) 
        Tacho::Blas<value_type>::trsv(handle_blas,
                                      CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_N, CUBLAS_DIAG_UNIT,
                                      m,
                                      A_ptr, m,
                                      x_ptr, 1);
#else
        Tacho::Blas<value_type>::trsv('L', 'N', 'U',
                                      m,
                                      A_ptr, m,
                                      x_ptr, 1);
#endif      
        Kokkos::fence();
        {
          auto x_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), x);
          printf("x (solve first) = \n");
          for (int i=0;i<m;++i) {
            printf("%e\n", x_host(i,0));
          }
        }

        Kokkos::parallel_for
          (Kokkos::RangePolicy<typename device_type::execution_space>(0,m),
           KOKKOS_LAMBDA(const int i) {
            if (p(i) == 0) { 
              /// do nothing
            } else if (p(i) < 0) {
              /// take 2x2 block to D
              const value_type 
                a00 = D(i-1, 0), a01 = D(i-1, 1),
                a10 = D(i  , 0), a11 = D(i  , 1);
              const value_type 
                det = a00*a11-a10*a01;
              const value_type 
                x0 = x(i-1,0),
                x1 = x(i,0);
              
              printf("i %d; %e %e; %e %e\n", a11/det, -a10/det, -a10/det, a00/det);
              x(i-1,0) = ( a11*x0 - a10*x1)/det;
              x(i  ,0) = (-a10*x0 + a00*x1)/det;
            } else {
              const value_type
                a00 = D(i,0);
              x(i,0) /= a00;
            }
          });
        Kokkos::fence();
        {
          auto x_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), x);
          printf("x (inverse scale diagoanl) = \n");
          for (int i=0;i<m;++i) {
            printf("%e\n", x_host(i,0));
          }
        }

#if defined (KOKKOS_ENABLE_CUDA) 
        Tacho::Blas<value_type>::trsv(handle_blas,
          CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_T, CUBLAS_DIAG_UNIT,
          m,
          A_ptr, m,
          x, 1);      
#else
        Tacho::Blas<value_type>::trsv('L', 'T', 'U',
          m,
          A_ptr, m,
          x_ptr, 1);            
#endif
        Kokkos::fence();
        {
          auto x_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), x);
          printf("x (solve b) = \n");
          for (int i=0;i<m;++i) {
            printf("%e\n", x_host(i,0));
          }
        }
        

        /// permute back to z
        Kokkos::View<value_type**,device_type> z("z",m,1);
        Kokkos::parallel_for
          (Kokkos::RangePolicy<typename device_type::execution_space>(0,m),
           KOKKOS_LAMBDA(const int i) {           
            z(i,0) = x(peri(i),0);
          });
        
        {
          auto x_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), z);
          printf("x (permute) = \n");
          for (int i=0;i<m;++i) {
            printf("%e\n", x_host(i,0));
          }
        }

      }
    }



#if defined (KOKKOS_ENABLE_CUDA) 
    {
      if (verbose) printf("cublas/cusolver handle destroy begin\n");
      {
        const int status = cublasDestroy(handle_blas);
        if (status) printf("Nonzero error from cublasDestroy %d\n", status);
      }
      {
        const int status = cusolverDnDestroy(handle_lapack);
        if (status) printf("Nonzero error from cusolverDnDestroy %d\n", status);
      }
      if (verbose) printf("cublas/cusolver handle destroy end\n");
    }
#endif
  }

  return r_val;
}

#endif
