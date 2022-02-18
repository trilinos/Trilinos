#ifndef __TACHO_EXAMPLE_DEVICE_DENSE_LDL_HPP__
#define __TACHO_EXAMPLE_DEVICE_DENSE_LDL_HPP__

#include "Kokkos_Random.hpp"

#include "Tacho_Util.hpp"
#include "Tacho_Blas_External.hpp"
#include "Tacho_Lapack_External.hpp"

#include "Tacho_LDL.hpp"
#include "Tacho_LDL_External.hpp"
#include "Tacho_LDL_Internal.hpp"
#include "Tacho_LDL_OnDevice.hpp"

#include "Tacho_Scale2x2_BlockInverseDiagonals.hpp"
#include "Tacho_Scale2x2_BlockInverseDiagonals_OnDevice.hpp"

#include "Tacho_Gemv.hpp"
#include "Tacho_Gemv_OnDevice.hpp"

#include "Tacho_Trsv.hpp"
#include "Tacho_Trsv_OnDevice.hpp"

#include "Tacho_ApplyPermutation.hpp"
#include "Tacho_ApplyPermutation_OnDevice.hpp"

template<typename value_type>
int driver_ldl (const int m, const bool verbose) {
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

    Kokkos::DefaultExecutionSpace exec_instance;

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
    int handle_blas(0), handle_lapack(0); // dummy
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

      //value_type * A_ptr = A.data();
      //value_type * x_ptr = x.data();
      //int * p_ptr = p.data();
        
      timer.reset();
      Kokkos::View<value_type*,Kokkos::LayoutLeft,device_type> W;
      const int lwork = Tacho::LDL<Tacho::Uplo::Lower,Tacho::Algo::OnDevice>
        ::invoke(handle_lapack, A, p, W);
      W = Kokkos::View<value_type*,Kokkos::LayoutLeft,device_type>("W", lwork);

      printf("LDLt lwork %d\n", lwork);
      Tacho::LDL<Tacho::Uplo::Lower,Tacho::Algo::OnDevice>
        ::invoke(handle_lapack, A, p, W);
      // {
      //   using policy_type = Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>;
      //   policy_type policy(1, 1, 1);
      //   Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const typename policy_type::member_type &member) {
      //       Tacho::LDL<Tacho::Uplo::Lower,Tacho::Algo::Internal>
      //         ::invoke(member, A, p, W);
      //     });
      // }
      Kokkos::fence();
      
      Tacho::LDL<Tacho::Uplo::Lower,Tacho::Algo::OnDevice>
        ::modify(exec_instance, A, p, D);
      Kokkos::fence();

      const double t = timer.seconds();
        
      auto AA = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A);
      auto pp = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), p);
      auto DD = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), D);
      {
        printf("L  = \n");
        for (int i=0;i<m;++i) {
          for (int j=0;j<m;++j)
            printf("%e ", j<=i ? A(i,j) : zero);
          printf("\n");
        }
        printf("D, lapack piv, flame piv, perm, peri = \n");
        for (int i=0;i<m;++i) {
          printf("% e % e ", DD(i,0), DD(i,1));
          printf("%4d %4d %4d %4d\n", pp(i),pp(i+m),pp(i+2*m),pp(i+3*m));
        }
        
      }
      printf("LDL time %e\n", t);

      Kokkos::fence();
        
      /// solve 
      {
        Kokkos::fence();
        auto perm = Kokkos::subview(p, Kokkos::pair<int,int>(2*m, 3*m));
        auto peri = Kokkos::subview(p, Kokkos::pair<int,int>(3*m, 4*m));
        /// copy and transpose
        Tacho::ApplyPermutation<Tacho::Side::Left,Tacho::Trans::NoTranspose,Tacho::Algo::OnDevice>
          ::invoke(exec_instance, b, perm, x);        

        {
          auto x_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), x);
          printf("x (permuted b) = \n");
          for (int i=0;i<m;++i) {
            printf("%e\n", x_host(i,0));
          }
        }
        
        /// trsv
        Tacho::Trsv<Tacho::Uplo::Lower,Tacho::Trans::NoTranspose,Tacho::Algo::OnDevice>
          ::invoke(handle_blas, Tacho::Diag::Unit(), A, x);
        Kokkos::fence();
        {
          auto x_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), x);
          printf("x (solve first) = \n");
          for (int i=0;i<m;++i) {
            printf("%e\n", x_host(i,0));
          }
        }

        Tacho::Scale2x2_BlockInverseDiagonals<Tacho::Side::Left,Tacho::Algo::OnDevice>
          ::invoke(exec_instance, p, D, x);
        Kokkos::fence();
        {
          auto x_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), x);
          printf("x (inverse scale diagoanl) = \n");
          for (int i=0;i<m;++i) {
            printf("%e\n", x_host(i,0));
          }
        }

        Tacho::Trsv<Tacho::Uplo::Lower,Tacho::Trans::Transpose,Tacho::Algo::OnDevice>
          ::invoke(handle_blas, Tacho::Diag::Unit(), A, x);
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
        Tacho::ApplyPermutation<Tacho::Side::Left,Tacho::Trans::NoTranspose,Tacho::Algo::OnDevice>
          ::invoke(exec_instance, x, peri, z);        

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
