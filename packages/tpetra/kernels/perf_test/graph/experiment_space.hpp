//#define CUDACONFIG
#define OPENMP


#ifdef CUDACONFIG
#define REPEAT 5
typedef Kokkos::Cuda MyExecSpace;
typedef Kokkos::Cuda MyMemorySpace;
typedef Kokkos::CudaUVMSpace MyEdgeMemorySpace;

typedef Kokkos::CudaHostPinnedSpace TemporaryWorkSpace;
typedef Kokkos::CudaUVMSpace PersistentWorkSpace;

#ifdef CUDAPINNEDSPACE
  typedef Kokkos::CudaHostPinnedSpace MyEdgeMemorySpace;
#else

#endif

#endif

#ifdef OPENMP
#define REPEAT 10
typedef Kokkos::OpenMP MyExecSpace;
typedef Kokkos::OpenMP MyMemorySpace;
typedef Kokkos::OpenMP MyEdgeMemorySpace;

typedef Kokkos::OpenMP TemporaryWorkSpace;
typedef Kokkos::OpenMP PersistentWorkSpace;
#endif



typedef int idx;
typedef double wt;
typedef unsigned int color_type;
typedef unsigned long long int color_type_eb;


typedef Kokkos::View<idx *, MyMemorySpace> idx_array_type;
typedef Kokkos::View<idx *, MyEdgeMemorySpace> idx_edge_array_type;
typedef Kokkos::View<wt *, MyEdgeMemorySpace> value_array_type;

typedef Kokkos::View<color_type_eb *, MyMemorySpace> color_eb_array_type;
typedef Kokkos::View<color_type * , MyMemorySpace> color_array_type;

typedef Kokkos::View<idx *, MyMemorySpace::array_layout, Kokkos::Serial, Kokkos::MemoryUnmanaged> um_array_type;
typedef Kokkos::View<idx *, MyMemorySpace::array_layout, MyEdgeMemorySpace, Kokkos::MemoryUnmanaged> um_edge_array_type;


typedef Kokkos::View<wt *, MyMemorySpace::array_layout, MyEdgeMemorySpace, Kokkos::MemoryUnmanaged> wt_um_edge_array_type;


typedef Kokkos::View<wt *, MyMemorySpace::array_layout, Kokkos::Serial, Kokkos::MemoryUnmanaged> um_value_array_type;


