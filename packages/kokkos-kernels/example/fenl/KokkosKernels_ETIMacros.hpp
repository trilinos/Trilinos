#ifdef MAKE_BUILD
#ifdef KOKKOS_ENABLE_CUDA
  #define KOKKOSKERNELS_ETI_MANGLING_TYPEDEFS()  \
                        typedef Kokkos::Device<Kokkos::Cuda, Kokkos::Cuda::memory_space> Kokkos_Device0Kokkos_Cuda_Kokkos_CudaSpace0; \
        typedef Kokkos::complex<double> Kokkos_complex0double0; \
        typedef long long longlong;
#else
  #ifdef KOKKOS_ENABLE_OPENMP
    #define KOKKOSKERNELS_ETI_MANGLING_TYPEDEFS()  \
                        typedef Kokkos::Device<Kokkos::OpenMP, Kokkos::OpenMP::memory_space> Kokkos_Device0Kokkos_OpenMP_Kokkos_HostSpace0; \
        typedef Kokkos::complex<double> Kokkos_complex0double0; \
        typedef long long longlong;
  #else
    #ifdef KOKKOS_ENABLE_THREADS
      #define KOKKOSKERNELS_ETI_MANGLING_TYPEDEFS()  \
                        typedef Kokkos::Device<Kokkos::Threads, Kokkos::Threads::memory_space> Kokkos_Device0Kokkos_Threads_Kokkos_HostSpace0; \
        typedef Kokkos::complex<double> Kokkos_complex0double0; \
        typedef long long longlong;
    #else
      #define KOKKOSKERNELS_ETI_MANGLING_TYPEDEFS()  \
                        typedef Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace> Kokkos_Device0Kokkos_OpenMP_Kokkos_HostSpace0; \
        typedef Kokkos::complex<double> Kokkos_complex0double0; \
        typedef long long longlong;
    #endif
  #endif
#endif

#endif

