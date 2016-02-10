#include <KokkosKernels_Utils.hpp>
#include <cstdlib>
#include <iostream>

#include "KokkosKernels_GraphHelpers.hpp"

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





int main (int argc, char ** argv){
  if (argc < 2){
    std::cerr << "Usage:" << argv[0] << " input_bin_file" << std::endl;
    exit(1);
  }


  Kokkos::initialize(argc, argv);
  MyExecSpace::print_configuration(std::cout);
  idx nv = 0, ne = 0;
  idx *xadj, *adj;
  wt *ew;

  KokkosKernels::Experimental::Graph::Utils::read_graph_bin<idx, wt> (&nv, &ne, &xadj, &adj, &ew, argv[1]);

  std::cout << "nv:" << nv << " ne:" << ne << std::endl;

  um_array_type _xadj (xadj, nv + 1);
  um_edge_array_type _adj (adj, ne);

  idx_array_type kok_xadj ("xadj", nv + 1);
  idx_edge_array_type kok_adj("adj", ne);

  idx_array_type sym_xadj;
  idx_edge_array_type sym_adj;

  Kokkos::deep_copy (kok_xadj, _xadj);
  Kokkos::deep_copy (kok_adj, _adj);

  wt_um_edge_array_type _mtx_vals (ew, ne);
  value_array_type kok_mtx_vals ("MTX_VALS", ne);
  Kokkos::deep_copy (kok_mtx_vals, _mtx_vals);

  delete [] xadj;
  delete [] adj;
  delete [] ew;

  std::cout << "Symetrizing Graph" << std::endl;

  Kokkos::Impl::Timer timer;

  KokkosKernels::Experimental::Util::symmetrize_graph_symbolic_hashmap<
  idx_array_type, idx_edge_array_type, idx_array_type, idx_edge_array_type, MyExecSpace>
    (nv, kok_xadj, kok_adj,sym_xadj, sym_adj);
  Kokkos::fence();

  double t = timer.seconds();
  std::cout << "Time to symmetrize:" << t << std::endl;
  KokkosKernels::Experimental::Util::print_1Dview(kok_xadj);
  KokkosKernels::Experimental::Util::print_1Dview(kok_adj);

  std::cout << "Symetric Graph" << std::endl;
  KokkosKernels::Experimental::Util::print_1Dview(sym_xadj);

  KokkosKernels::Experimental::Util::print_1Dview(sym_adj);

  Kokkos::finalize();



  return 0;
}
