#ifndef __TACHO_EXAMPLE_DEVICE_DENSE_LDL_HPP__
#define __TACHO_EXAMPLE_DEVICE_DENSE_LDL_HPP__

#include "Kokkos_Random.hpp"

#include "Tacho_Util.hpp"
#include "Tacho_Blas_External.hpp"
#include "Tacho_Lapack_External.hpp"

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
      auto AA = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), Aback);
      {
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
          x(i) = i+1;
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
#if defined (KOKKOS_ENABLE_CUDA) 
      int lwork(0);
      Tacho::Lapack<value_type>::sytrf_buffersize(handle_lapack, 
                                                  m, 
                                                  A_ptr, m, 
                                                  &lwork);
      printf("LDLt lwork %d\n", lwork);
      Kokkos::View<value_type*,Kokkos::LayoutLeft,device_type> W("W", lwork);
      value_type * W_ptr = W.data();
      
      Tacho::Lapack<value_type>::sytrf(handle_lapack,
                                       CUBLAS_FILL_MODE_LOWER,
                                       m, 
                                       A_ptr, m,
                                       p_ptr,
                                       W_ptr, lwork,
                                       dev.data());
#else
      int lwork(m*32);
      printf("LDLt lwork %d\n", lwork);
      Kokkos::View<value_type*,Kokkos::LayoutLeft,device_type> W("W", lwork);
      value_type * W_ptr = W.data();
      Tacho::Lapack<value_type>::sytrf('L',
                                       m,
                                       A_ptr, m,
                                       p_ptr,
                                       W_ptr, lwork,
                                       dev.data());

        printf("LDL org  = \n");
        for (int i=0;i<m;++i) {
          for (int j=0;j<m;++j)
            std::cout << A(i,j) << " ";
          printf("\n");
        }

      /// convert lapack to flame and perm vector
      {
        int 
          * flame_pivots = p_ptr + m, 
          * perm = flame_pivots + m, 
          * peri = perm + m;
        for (int i=0;i<m;++i) perm[i] = i;
        for (int i=0,cnt=0;i<m;++i) {
          if (p_ptr[i] < 0) {
            if (++cnt%2) {
              p_ptr[i] = 0; /// invalidate the pivot
              flame_pivots[i] = 0;
              D(i, 0) = A(i,i);     D(i, 1) = A(i+1,i);
              A(i,i) = one;
            } else {
              const int fla_pivot = -p_ptr[i] - i - 1;
              flame_pivots[i] = fla_pivot;
              if (fla_pivot) {
                value_type * src = A_ptr + i;
                value_type * tgt = src + fla_pivot;
                printf("i %d, fla pivot %d\n", i, fla_pivot);
                for (int j=0;j<(i-1);++j) {
                  const int idx = j*m;
                  printf("- j %d, src %e, tgt %e\n", j, src[idx], tgt[idx]);
                  const value_type tmp = src[idx];
                  src[idx] = tgt[idx];
                  tgt[idx] = tmp;                  
                  printf("+ j %d, src %e, tgt %e\n", j, src[idx], tgt[idx]);
                }
              }              
              D(i, 0) = A(i,i-1);     D(i, 1) = A(i,i);
              A(i,i-1) = zero; A(i,i) = one;
            }
          } else {
            const int fla_pivot = p_ptr[i] - i - 1;
            flame_pivots[i] = fla_pivot;
            if (fla_pivot) {
              printf("i %d, fla pivot %d\n", i, fla_pivot);
              value_type * src = A_ptr + i;
              value_type * tgt = src + fla_pivot;
              for (int j=0;j<i;++j) {
                const int idx = j*m;
                printf("- j %d, src %e, tgt %e\n", j, src[idx], tgt[idx]);
                const value_type tmp = src[idx];
                src[idx] = tgt[idx];
                tgt[idx] = tmp;     
                printf("- j %d, src %e, tgt %e\n", j, src[idx], tgt[idx]);             
              }
            }

            D(i, 0) = A(i,i); 
            A(i,i) = one;
          }

          /// apply pivots to perm vector
          if (flame_pivots[i]) {
            const int tmp = perm[i], pidx = i+flame_pivots[i];
            perm[i] = perm[pidx];
            perm[pidx] = tmp;
          }
        }
        for (int i=0;i<m;++i) peri[perm[i]] = i;
      }

      // extract and apply pivots

#endif
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
            x(i) = b(perm(i));
          });
        {
          auto x_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), x);
          printf("x (permuted b) = \n");
          for (int i=0;i<m;++i) {
            printf("%e\n", x_host(i));
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
            printf("%e\n", x_host(i));
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
                x0 = x(i-1),
                x1 = x(i);
              
              printf("i %d; %e %e; %e %e\n", a11/det, -a10/det, -a10/det, a00/det);
              x(i-1) = ( a11*x0 - a10*x1)/det;
              x(i  ) = (-a10*x0 + a00*x1)/det;
            } else {
              const value_type
                a00 = D(i,0);
              x(i) /= a00;
            }
          });
        Kokkos::fence();
        {
          auto x_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), x);
          printf("x (inverse scale diagoanl) = \n");
          for (int i=0;i<m;++i) {
            printf("%e\n", x_host(i));
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
            printf("%e\n", x_host(i));
          }
        }
        

        /// permute back to z
        Kokkos::View<value_type*,device_type> z("z",m);
        Kokkos::parallel_for
          (Kokkos::RangePolicy<typename device_type::execution_space>(0,m),
           KOKKOS_LAMBDA(const int i) {           
            z(i) = x(peri(i));
          });
        
        {
          auto x_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), z);
          printf("x (permute) = \n");
          for (int i=0;i<m;++i) {
            printf("%e\n", x_host(i));
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
