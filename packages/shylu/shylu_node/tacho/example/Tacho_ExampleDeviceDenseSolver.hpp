#include "Kokkos_Random.hpp"

#include "Tacho_Util.hpp"
#include "Tacho_Blas_External.hpp"
#include "Tacho_Lapack_External.hpp"
#include "Tacho_CommandLineParser.hpp" 

using ordinal_type = Tacho::ordinal_type;

template<typename value_type>
int driver (int argc, char *argv[]) {
  int nthreads = 1;
  bool verbose = true;

  int max_iter = 1;
  int m = 10;

  Tacho::CommandLineParser opts("This example program measure the Tacho on Kokkos::OpenMP");

  opts.set_option<int>("kokkos-threads", "Number of threads", &nthreads);
  opts.set_option<bool>("verbose", "Flag for verbose printing", &verbose);
  opts.set_option<int>("m", "Dense problem size", &m);

  const bool r_parse = opts.parse(argc, argv);
  if (r_parse) return 0; // print help return

  Kokkos::initialize(argc, argv);
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
      if (verbose) printf("cublas/cusolver handle create begin\n");
      {
        const int status = cublasCreate(&handle_blas); 
        if (status) printf("Nonzero error from cublasCreate %d\n", status);
      }
      {
        const int status = cusolverDnCreate(&handle_lapack); 
        if (status) printf("Nonzero error from cusolverDnCreate %d\n", status);
      }
      if (verbose) printf("cublas/cusolver handle create end\n");
    }
#endif

    Kokkos::View<value_type**,Kokkos::LayoutLeft,device_type> Arand("Arand", m, m);
    {
      if (verbose) printf("test problem randomization\n");
      Kokkos::Random_XorShift64_Pool<device_type> random(13718);
      Kokkos::fill_random(Arand, random, one);
    }

    Kokkos::View<value_type**,Kokkos::LayoutLeft,device_type> A("A", m, m);
    {
      if (verbose) printf("test problem symmetrization\n");
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
      printf("symmetrization time %e\n", t);
    } 

    Kokkos::View<value_type**,Kokkos::LayoutLeft,device_type> Aback("Aback", m, m);    
    {
      if (verbose) printf("test problem backup\n");
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

    /// run Cholesky
    for (int iter=0;iter<max_iter;++iter) {
      Kokkos::deep_copy(A, Aback);
      {
        value_type * A_ptr = A.data();
        int dev(0);
        
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
        
        Kokkos::View<int,device_type> dev_view("dev");
        Tacho::Lapack<value_type>::potrf(handle_lapack,
                                         CUBLAS_FILL_MODE_UPPER,
                                         m, 
                                         A_ptr, m,
                                         W_ptr, lwork,
                                         dev_view.data());
        auto dev_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), dev_view);
        dev = dev_host();
#else
        Tacho::Lapack<value_type>::potrf('U',
                                         m,
                                         A_ptr, m,
                                         &dev);
#endif
        Kokkos::fence();
        const double t = timer.seconds();
        printf("Cholesky time %e\n", t);
        if (dev) printf("Cholesky returns non-zero dev info %d\n", dev);
      }
    }

    /// run LDLt
    for (int iter=0;iter<max_iter;++iter) {
      Kokkos::deep_copy(A, Aback);
      Kokkos::View<int*,Kokkos::LayoutLeft,device_type> ipiv("pivot", m);
      {
        value_type * A_ptr = A.data();
        int * ipiv_ptr = ipiv.data();
        int dev(0);
        
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
        
        Kokkos::View<int,device_type> dev_view("dev");
        Tacho::Lapack<value_type>::sytrf(handle_lapack,
                                         CUBLAS_FILL_MODE_LOWER,
                                         m, 
                                         A_ptr, m,
                                         ipiv_ptr,
                                         W_ptr, lwork,
                                         dev_view.data());
        auto dev_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), dev_view);
        dev = dev_host();
#else
        int lwork(m*32);
        printf("LDLt lwork %d\n", lwork);
        Kokkos::View<value_type*,Kokkos::LayoutLeft,device_type> W("W", lwork);
        value_type * W_ptr = W.data();
        Tacho::Lapack<value_type>::sytrf('L',
                                         m,
                                         A_ptr, m,
                                         ipiv_ptr,
                                         W_ptr, lwork,
                                         &dev);
#endif
        Kokkos::fence();
        const double t = timer.seconds();
        auto AA = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A);
        auto pp = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), ipiv);
        if (0) {
          printf("LDL  = \n");
          for (int i=0;i<m;++i) {
            for (int j=0;j<m;++j)
              std::cout << AA(i,j) << " ";
            printf("\n");
          }
          printf("piv = \n");
          for (int i=0;i<m;++i) {
            printf("%d\n", pp(i));
          }
        }
        
        printf("LDLt time %e\n", t);
        if (dev) printf("LDLt returns non-zero dev info %d\n", dev);
      }
    }


    /// run LU
    for (int iter=0;iter<4;++iter) {
      Kokkos::deep_copy(A, Aback);
      Kokkos::View<int*,Kokkos::LayoutLeft,device_type> ipiv("pivot", m);
      {
        value_type * A_ptr = A.data();
        int * ipiv_ptr = ipiv.data();
        int dev(0);
        
        timer.reset();
#if defined (KOKKOS_ENABLE_CUDA) 
        int lwork(0);
        Tacho::Lapack<value_type>::getrf_buffersize(handle_lapack, 
                                                    m, m,
                                                    A_ptr, m, 
                                                    &lwork);
        printf("LU lwork %d\n", lwork);
        Kokkos::View<value_type*,Kokkos::LayoutLeft,device_type> W("W", lwork);
        value_type * W_ptr = W.data();
        
        Kokkos::View<int,device_type> dev_view("dev");
        Tacho::Lapack<value_type>::getrf(handle_lapack,
                                         m, m,
                                         A_ptr, m,
                                         W_ptr,
                                         ipiv_ptr,
                                         dev_view.data());
        auto dev_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), dev_view);
        dev = dev_host();
#else
        Tacho::Lapack<value_type>::getrf(m, m,
                                         A_ptr, m,
                                         ipiv_ptr,
                                         &dev);
#endif
        Kokkos::fence();
        const double t = timer.seconds();
        printf("LU time %e\n", t);
        if (dev) printf("LU returns non-zero dev info %d\n", dev);
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
  Kokkos::finalize();

  return r_val;
}
