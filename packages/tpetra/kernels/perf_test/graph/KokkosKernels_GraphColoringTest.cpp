#include <KokkosKernels_GraphColor.hpp>
#include <KokkosKernels_Handle.hpp>

#include <cstdlib>
#include <iostream>

#include <random>       // std::default_random_engine
#include <algorithm>    // std::shuffle
#include <vector>

#include <KokkosKernels_GraphHelpers.hpp>

int main (int argc, char ** argv){
  if (argc < 2){
    std::cerr << "Usage:" << argv[0] << " input_bin_file" << std::endl;
    exit(1);
  }
  Kokkos::initialize(argc, argv);


  const int numColoringAlgos = 6;
  using namespace KokkosKernels::Experimental::Graph;
  const ColoringAlgorithm ColoringAlgorithms[] = {COLORING_DEFAULT, COLORING_SERIAL, COLORING_VB, COLORING_VBBIT, COLORING_VBCS, COLORING_EB, COLORING_EB};
  const std::string ColoringAlgorithmNames[] = {"COLORING_DEFAULT", "COLORING_SERIAL", "COLORING_VB", "COLORING_VBBIT", "COLORING_VBCS", "COLORING_EB"};


  typedef int idx;
  typedef double wt;


  idx nr = 0, ne = 0;
  idx *xadj, *adj;
  wt *ew;

  KokkosKernels::Experimental::Graph::Utils::read_graph_bin<idx, wt> (
      &nr, &ne, &xadj, &adj, &ew, argv[1]);
  delete [] ew;

  idx nc = 0;
  for (idx i = 0; i < ne; ++i){
    if (adj[i] > nc) nc = adj[i];
  }

  nc += 1;

#if defined( KOKKOS_HAVE_SERIAL )


  {


    std::cout << "SERIAL TEST 1" << std::endl;
    typedef int lno_t;
    typedef double scalar_t;
    //typedef unsigned int color_t;

    typedef Kokkos::Serial ExecSpace;
    typedef Kokkos::Serial RowMemorySpace;
    typedef Kokkos::Serial NonzeroMemorySpace;

    ExecSpace::print_configuration(std::cout);

    typedef Kokkos::Serial::memory_space TempWorkSpace;
    typedef Kokkos::Serial::memory_space PersistentWorkSpace;

    typedef Kokkos::View<lno_t *, RowMemorySpace> row_index_view_type;
    typedef Kokkos::View<lno_t *, NonzeroMemorySpace> nonzero_index_view_type;
    typedef Kokkos::View<scalar_t *, NonzeroMemorySpace> nonzero_scalar_view_type;

    //typedef Kokkos::View<color_t * , RowMemorySpace> color_view_type;

    typedef Kokkos::View<idx *, RowMemorySpace::array_layout, Kokkos::Serial, Kokkos::MemoryUnmanaged> um_array_type;
    typedef Kokkos::View<idx *, RowMemorySpace::array_layout, NonzeroMemorySpace, Kokkos::MemoryUnmanaged> um_edge_array_type;

    um_array_type _xadj (xadj, nr + 1);
    um_edge_array_type _adj (adj, ne);


    row_index_view_type kok_xadj ("xadj", nr + 1);
    nonzero_index_view_type kok_adj("adj", ne);
    KokkosKernels::Experimental::Util::copy_vector<um_array_type, row_index_view_type, ExecSpace>(nr+1, _xadj, kok_xadj );
    KokkosKernels::Experimental::Util::copy_vector<um_edge_array_type, nonzero_index_view_type, ExecSpace>(ne, _adj, kok_adj );



    typedef KokkosKernels::Experimental::KokkosKernelsHandle
        <row_index_view_type,nonzero_index_view_type, nonzero_scalar_view_type, ExecSpace, TempWorkSpace,PersistentWorkSpace > KernelHandle;
    KernelHandle kkh;

    for (int i = 0; i < numColoringAlgos; ++i){
      std::cout << "\t:" << "Running " << ColoringAlgorithmNames[i] << std::endl;
      kkh.create_graph_coloring_handle(ColoringAlgorithms[i]);
      KokkosKernels::Experimental::Graph::graph_color_symbolic
      <KernelHandle,row_index_view_type,nonzero_index_view_type> (&kkh,nr, nc, kok_xadj, kok_adj);

      std::cout << "\t:" << ColoringAlgorithmNames[i] <<
          " Time:" << kkh.get_graph_coloring_handle()->get_overall_coloring_time() << " "
          "Num colors:" << kkh.get_graph_coloring_handle()->get_num_colors() << " "
          "Num Phases:" << kkh.get_graph_coloring_handle()->get_num_phases() << std::endl;
      std::cout << "\t"; KokkosKernels::Experimental::Util::print_1Dview(kkh.get_graph_coloring_handle()->get_vertex_colors());
    }

    /*
      typedef typename KernelHandle::row_index_temp_work_view_type row_view_type;
      typedef typename KernelHandle::row_index_persistent_work_view_type edge_view_type;
      row_view_type sym_row_map;
      edge_view_type sym_entries;

      KokkosKernels::Experimental::Util::symmetrize_graph_symbolic
        < row_index_view_type, nonzero_index_view_type,row_view_type, edge_view_type,ExecSpace>
        (nr, kok_xadj, kok_adj, sym_row_map,sym_entries );
      kok_xadj = row_index_view_type("s");
      kok_adj = nonzero_index_view_type("s");

      for (int i = 0; i < numColoringAlgos; ++i){
        std::cout << "\t:" << "On Temp Memory Space Running " << ColoringAlgorithmNames[i] << std::endl;
        kkh.create_graph_coloring_handle(ColoringAlgorithms[i]);
        KokkosKernels::Experimental::Graph::graph_color_symbolic
            <KernelHandle,row_view_type,edge_view_type> (&kkh,nr, nc, sym_row_map, sym_entries);

        std::cout << "\t:" << ColoringAlgorithmNames[i] <<
            " Time:" << kkh.get_graph_coloring_handle()->get_overall_coloring_time() << " "
            "Num colors:" << kkh.get_graph_coloring_handle()->get_num_colors() << " "
            "Num Phases:" << kkh.get_graph_coloring_handle()->get_num_phases() << std::endl;
        std::cout << "\t"; KokkosKernels::Experimental::Util::print_1Dview(kkh.get_graph_coloring_handle()->get_vertex_colors());
      }
     */
  }
#endif
#if defined( KOKKOS_HAVE_OPENMP )


  {
    //OPENMP TEST1
    {

      std::cout << "OPENMP TEST 1" << std::endl;
      typedef int lno_t;
      typedef double scalar_t;
      //typedef unsigned int color_t;

      typedef Kokkos::OpenMP ExecSpace;
      typedef Kokkos::OpenMP RowMemorySpace;
      typedef Kokkos::OpenMP::memory_space NonzeroMemorySpace;

      ExecSpace::print_configuration(std::cout);

      typedef Kokkos::OpenMP TempWorkSpace;
      typedef Kokkos::OpenMP PersistentWorkSpace;

      typedef Kokkos::View<lno_t *, RowMemorySpace> row_index_view_type;
      typedef Kokkos::View<lno_t *, NonzeroMemorySpace> nonzero_index_view_type;
      typedef Kokkos::View<scalar_t *, NonzeroMemorySpace> nonzero_scalar_view_type;

      //typedef Kokkos::View<color_t * , RowMemorySpace> color_view_type;

      typedef Kokkos::View<idx *, RowMemorySpace::array_layout, Kokkos::Serial, Kokkos::MemoryUnmanaged> um_array_type;
      typedef Kokkos::View<idx *, RowMemorySpace::array_layout, NonzeroMemorySpace, Kokkos::MemoryUnmanaged> um_edge_array_type;

      um_array_type _xadj (xadj, nr + 1);
      um_edge_array_type _adj (adj, ne);


      row_index_view_type kok_xadj ("xadj", nr + 1);
      nonzero_index_view_type kok_adj("adj", ne);
      KokkosKernels::Experimental::Util::copy_vector<um_array_type, row_index_view_type, ExecSpace>(nr+1, _xadj, kok_xadj );
      KokkosKernels::Experimental::Util::copy_vector<um_edge_array_type, nonzero_index_view_type, ExecSpace>(ne, _adj, kok_adj );



      typedef KokkosKernels::Experimental::KokkosKernelsHandle
              <row_index_view_type,nonzero_index_view_type, nonzero_scalar_view_type, ExecSpace, TempWorkSpace,PersistentWorkSpace > KernelHandle;
      KernelHandle kkh;

      for (int i = 0; i < numColoringAlgos; ++i){
        std::cout << "\t:" << "Running " << ColoringAlgorithmNames[i] << std::endl;
        kkh.create_graph_coloring_handle(ColoringAlgorithms[i]);
        KokkosKernels::Experimental::Graph::graph_color_symbolic
              <KernelHandle,row_index_view_type,nonzero_index_view_type> (&kkh,nr, nc, kok_xadj, kok_adj);

        std::cout << "\t:" << ColoringAlgorithmNames[i] <<
            " Time:" << kkh.get_graph_coloring_handle()->get_overall_coloring_time() << " "
            "Num colors:" << kkh.get_graph_coloring_handle()->get_num_colors() << " "
            "Num Phases:" << kkh.get_graph_coloring_handle()->get_num_phases() << std::endl;
        std::cout << "\t"; KokkosKernels::Experimental::Util::print_1Dview(kkh.get_graph_coloring_handle()->get_vertex_colors());
      }

      /*
      typedef typename KernelHandle::row_index_temp_work_view_type row_view_type;
      typedef typename KernelHandle::row_index_persistent_work_view_type edge_view_type;
      row_view_type sym_row_map;
      edge_view_type sym_entries;

      KokkosKernels::Experimental::Util::symmetrize_graph_symbolic
        < row_index_view_type, nonzero_index_view_type,row_view_type, edge_view_type,ExecSpace>
        (nr, kok_xadj, kok_adj, sym_row_map,sym_entries );
      kok_xadj = row_index_view_type("s");
      kok_adj = nonzero_index_view_type("s");

      for (int i = 0; i < numColoringAlgos; ++i){
        std::cout << "\t:" << "On Temp Memory Space Running " << ColoringAlgorithmNames[i] << std::endl;
        kkh.create_graph_coloring_handle(ColoringAlgorithms[i]);
        KokkosKernels::Experimental::Graph::graph_color_symbolic
            <KernelHandle,row_view_type,edge_view_type> (&kkh,nr, nc, sym_row_map, sym_entries);

        std::cout << "\t:" << ColoringAlgorithmNames[i] <<
            " Time:" << kkh.get_graph_coloring_handle()->get_overall_coloring_time() << " "
            "Num colors:" << kkh.get_graph_coloring_handle()->get_num_colors() << " "
            "Num Phases:" << kkh.get_graph_coloring_handle()->get_num_phases() << std::endl;
        std::cout << "\t"; KokkosKernels::Experimental::Util::print_1Dview(kkh.get_graph_coloring_handle()->get_vertex_colors());
      }
      */
    }

    //OPENMP TEST2
    {

      std::cout << "OPENMP TEST 2" << std::endl;
      typedef unsigned long lno_t;
      typedef float scalar_t;
      //typedef int color_t;

      typedef Kokkos::OpenMP ExecSpace;
      typedef Kokkos::OpenMP RowMemorySpace;
      typedef Kokkos::OpenMP::memory_space NonzeroMemorySpace;
      ExecSpace::print_configuration(std::cout);

      typedef Kokkos::OpenMP TempWorkSpace;
      typedef Kokkos::OpenMP PersistentWorkSpace;

      typedef Kokkos::View<lno_t *, Kokkos::LayoutLeft, RowMemorySpace> row_index_view_type;
      typedef Kokkos::View<lno_t *, Kokkos::LayoutRight, NonzeroMemorySpace> nonzero_index_view_type;
      typedef Kokkos::View<scalar_t *, NonzeroMemorySpace> nonzero_scalar_view_type;

      typedef Kokkos::View<const lno_t *, Kokkos::LayoutLeft, RowMemorySpace> const_row_index_view_type;
      typedef Kokkos::View<const lno_t *, Kokkos::LayoutRight, NonzeroMemorySpace> const_nonzero_index_view_type;


      //typedef Kokkos::View<color_t * , Kokkos::LayoutStride, RowMemorySpace> color_view_type;

      typedef Kokkos::View<idx *, RowMemorySpace::array_layout, Kokkos::Serial, Kokkos::MemoryUnmanaged> um_array_type;
      typedef Kokkos::View<idx *, RowMemorySpace::array_layout, NonzeroMemorySpace, Kokkos::MemoryUnmanaged> um_edge_array_type;

      um_array_type _xadj (xadj, nr + 1);
      um_edge_array_type _adj (adj, ne);


      row_index_view_type kok_xadj ("xadj", nr + 1);
      nonzero_index_view_type kok_adj("adj", ne);
      KokkosKernels::Experimental::Util::copy_vector<um_array_type, row_index_view_type, ExecSpace>(nr+1, _xadj, kok_xadj );
      KokkosKernels::Experimental::Util::copy_vector<um_edge_array_type, nonzero_index_view_type, ExecSpace>(ne, _adj, kok_adj );

      typedef KokkosKernels::Experimental::KokkosKernelsHandle
              <const_row_index_view_type,const_nonzero_index_view_type, nonzero_scalar_view_type, ExecSpace, TempWorkSpace,PersistentWorkSpace > KernelHandle;
      KernelHandle kkh;

      for (int i = 0; i < numColoringAlgos; ++i){
        std::cout << "\t:" << "Running " << ColoringAlgorithmNames[i] << std::endl;
        kkh.create_graph_coloring_handle(ColoringAlgorithms[i]);
        KokkosKernels::Experimental::Graph::graph_color_symbolic
              <KernelHandle,const_row_index_view_type,const_nonzero_index_view_type> (&kkh,nr, nc, kok_xadj, kok_adj);

        std::cout << "\t:" << ColoringAlgorithmNames[i] <<
            " Time:" << kkh.get_graph_coloring_handle()->get_overall_coloring_time() << " "
            "Num colors:" << kkh.get_graph_coloring_handle()->get_num_colors() << " "
            "Num Phases:" << kkh.get_graph_coloring_handle()->get_num_phases() << std::endl;
        std::cout << "\t"; KokkosKernels::Experimental::Util::print_1Dview(kkh.get_graph_coloring_handle()->get_vertex_colors());
      }

      /*
      typedef typename KernelHandle::row_index_temp_work_view_type row_view_type;
      typedef typename KernelHandle::row_index_persistent_work_view_type edge_view_type;
      row_view_type sym_row_map;
      edge_view_type sym_entries;

      KokkosKernels::Experimental::Util::symmetrize_graph_symbolic
        < row_index_view_type, nonzero_index_view_type,row_view_type, edge_view_type,ExecSpace>
        (nr, kok_xadj, kok_adj, sym_row_map,sym_entries );
      kok_xadj = row_index_view_type("s");
      kok_adj = nonzero_index_view_type("s");

      for (int i = 0; i < numColoringAlgos; ++i){
        std::cout << "\t:" << "On Temp Memory Space Running " << ColoringAlgorithmNames[i] << std::endl;
        kkh.create_graph_coloring_handle(ColoringAlgorithms[i]);
        KokkosKernels::Experimental::Graph::graph_color_symbolic
            <KernelHandle,row_view_type,edge_view_type> (&kkh,nr, nc, sym_row_map, sym_entries);

        std::cout << "\t:" << ColoringAlgorithmNames[i] <<
            " Time:" << kkh.get_graph_coloring_handle()->get_overall_coloring_time() << " "
            "Num colors:" << kkh.get_graph_coloring_handle()->get_num_colors() << " "
            "Num Phases:" << kkh.get_graph_coloring_handle()->get_num_phases() << std::endl;
        std::cout << "\t"; KokkosKernels::Experimental::Util::print_1Dview(kkh.get_graph_coloring_handle()->get_vertex_colors());
      }
      */
    }

  }
#endif

#if defined( KOKKOS_HAVE_CUDA )

  {

    //CUDA TEST1
    {

      std::cout << "Cuda TEST 1" << std::endl;
      typedef int lno_t;
      typedef double scalar_t;
      typedef unsigned int color_t;

      typedef Kokkos::Cuda ExecSpace;
      typedef Kokkos::Cuda RowMemorySpace;
      typedef Kokkos::Cuda::memory_space NonzeroMemorySpace;
      ExecSpace::print_configuration(std::cout);

      typedef Kokkos::Cuda TempWorkSpace;
      typedef Kokkos::Cuda PersistentWorkSpace;

      typedef Kokkos::View<lno_t *, RowMemorySpace> row_index_view_type;
      typedef Kokkos::View<lno_t *, NonzeroMemorySpace> nonzero_index_view_type;
      typedef Kokkos::View<scalar_t *, NonzeroMemorySpace> nonzero_scalar_view_type;

      typedef Kokkos::View<color_t * , RowMemorySpace> color_view_type;

      typedef Kokkos::View<idx *, RowMemorySpace::array_layout, Kokkos::Serial, Kokkos::MemoryUnmanaged> um_array_type;
      typedef Kokkos::View<idx *, RowMemorySpace::array_layout, Kokkos::Serial, Kokkos::MemoryUnmanaged> um_edge_array_type;


      um_array_type _xadj (xadj, nr + 1);
      um_edge_array_type _adj (adj, ne);


      row_index_view_type kok_xadj ("xadj", nr + 1);
      nonzero_index_view_type kok_adj("adj", ne);
      Kokkos::deep_copy (kok_xadj , _xadj);
      Kokkos::deep_copy (kok_adj , _adj);
      //KokkosKernels::Experimental::Util::copy_vector<um_array_type, row_index_view_type, ExecSpace>(nr+1, _xadj, kok_xadj );
      //KokkosKernels::Experimental::Util::copy_vector<um_edge_array_type, nonzero_index_view_type, ExecSpace>(ne, _adj, kok_adj );

      typedef KokkosKernels::Experimental::KokkosKernelsHandle
              <row_index_view_type,nonzero_index_view_type, nonzero_scalar_view_type, ExecSpace, TempWorkSpace,PersistentWorkSpace > KernelHandle;
      KernelHandle kkh;

      for (int i = 0; i < numColoringAlgos; ++i){
        std::cout << "\t:" << "Running " << ColoringAlgorithmNames[i] << std::endl;
        kkh.create_graph_coloring_handle(ColoringAlgorithms[i]);
        KokkosKernels::Experimental::Graph::graph_color_symbolic
              <KernelHandle,row_index_view_type,nonzero_index_view_type> (&kkh,nr, nc, kok_xadj, kok_adj);

        std::cout << "\t:" << ColoringAlgorithmNames[i] <<
            " Time:" << kkh.get_graph_coloring_handle()->get_overall_coloring_time() << " "
            "Num colors:" << kkh.get_graph_coloring_handle()->get_num_colors() << " "
            "Num Phases:" << kkh.get_graph_coloring_handle()->get_num_phases() << std::endl;
        std::cout << "\t"; KokkosKernels::Experimental::Util::print_1Dview(kkh.get_graph_coloring_handle()->get_vertex_colors());
      }
      typedef typename KernelHandle::row_lno_temp_work_view_t row_view_type;

      typedef typename KernelHandle::row_lno_persistent_work_view_t edge_view_type;


      row_view_type sym_row_map;
      edge_view_type sym_entries;

      /*
      KokkosKernels::Experimental::Util::symmetrize_graph_symbolic
        < row_index_view_type, nonzero_index_view_type,row_view_type, edge_view_type,ExecSpace>
        (nr, kok_xadj, kok_adj, sym_row_map,sym_entries );
      kok_xadj = row_index_view_type("s");
      kok_adj = nonzero_index_view_type("s");

      for (int i = 0; i < numColoringAlgos; ++i){
        std::cout << "\t:" << "On Temp Memory Space Running " << ColoringAlgorithmNames[i] << std::endl;
        kkh.create_graph_coloring_handle(ColoringAlgorithms[i]);
        KokkosKernels::Experimental::Graph::graph_color_symbolic
            <KernelHandle,row_view_type,edge_view_type> (&kkh,nr, nc, sym_row_map, sym_entries);

        std::cout << "\t:" << ColoringAlgorithmNames[i] <<
            " Time:" << kkh.get_graph_coloring_handle()->get_overall_coloring_time() << " "
            "Num colors:" << kkh.get_graph_coloring_handle()->get_num_colors() << " "
            "Num Phases:" << kkh.get_graph_coloring_handle()->get_num_phases() << std::endl;
        std::cout << "\t"; KokkosKernels::Experimental::Util::print_1Dview(kkh.get_graph_coloring_handle()->get_vertex_colors());
      }
      */
    }
    //CUDA TEST2
    {

      std::cout << "Cuda TEST 2" << std::endl;
      typedef unsigned long lno_t;
      typedef float scalar_t;
      typedef int color_t;

      typedef Kokkos::Cuda ExecSpace;
      typedef Kokkos::Cuda RowMemorySpace;
      typedef Kokkos::Cuda::memory_space NonzeroMemorySpace;
      ExecSpace::print_configuration(std::cout);

      typedef Kokkos::Cuda TempWorkSpace;
      typedef Kokkos::Cuda PersistentWorkSpace;

      typedef Kokkos::View<lno_t *, Kokkos::LayoutLeft, RowMemorySpace> row_index_view_type;
      typedef Kokkos::View<lno_t *, Kokkos::LayoutRight, NonzeroMemorySpace> nonzero_index_view_type;
      typedef Kokkos::View<scalar_t *, NonzeroMemorySpace> nonzero_scalar_view_type;

      typedef Kokkos::View<color_t * , Kokkos::LayoutStride, RowMemorySpace> color_view_type;

      typedef Kokkos::View<idx *, RowMemorySpace::array_layout, Kokkos::Serial, Kokkos::MemoryUnmanaged> um_array_type;
      typedef Kokkos::View<idx *, RowMemorySpace::array_layout, NonzeroMemorySpace, Kokkos::MemoryUnmanaged> um_edge_array_type;

      um_array_type _xadj (xadj, nr + 1);
      um_edge_array_type _adj (adj, ne);

      typedef Kokkos::View<const lno_t *, Kokkos::LayoutLeft, RowMemorySpace> const_row_index_view_type;
      typedef Kokkos::View<const lno_t *, Kokkos::LayoutRight, NonzeroMemorySpace> const_nonzero_index_view_type;



      typedef Kokkos::View<idx *, RowMemorySpace> tmp_view_type;
      typedef Kokkos::View<idx *, NonzeroMemorySpace> tmp_edge_view_type;
      tmp_view_type __xadj ("xadj", nr + 1);
      tmp_edge_view_type __adj("adj", ne);

      Kokkos::deep_copy (__xadj , _xadj);
      Kokkos::deep_copy (__adj , _adj);

      row_index_view_type kok_xadj ("xadj", nr + 1);
      nonzero_index_view_type kok_adj("adj", ne);

      KokkosKernels::Experimental::Util::copy_vector<tmp_view_type, row_index_view_type, ExecSpace>(nr+1, __xadj, kok_xadj );
      KokkosKernels::Experimental::Util::copy_vector<tmp_edge_view_type, nonzero_index_view_type, ExecSpace>(ne, __adj, kok_adj );
      __xadj = tmp_view_type("");
      __adj = tmp_edge_view_type("");

      //KokkosKernels::Experimental::Util::copy_vector<um_array_type, row_index_view_type, ExecSpace>(nr+1, _xadj, kok_xadj );
      //KokkosKernels::Experimental::Util::copy_vector<um_edge_array_type, nonzero_index_view_type, ExecSpace>(ne, _adj, kok_adj );

      typedef KokkosKernels::Experimental::KokkosKernelsHandle
              <const_row_index_view_type,const_nonzero_index_view_type, nonzero_scalar_view_type, ExecSpace, TempWorkSpace,PersistentWorkSpace > KernelHandle;
      KernelHandle kkh;

      for (int i = 0; i < numColoringAlgos; ++i){
        std::cout << "\t:" << "Running " << ColoringAlgorithmNames[i] << std::endl;
        kkh.create_graph_coloring_handle(ColoringAlgorithms[i]);
        KokkosKernels::Experimental::Graph::graph_color_symbolic
              <KernelHandle,const_row_index_view_type,const_nonzero_index_view_type> (&kkh,nr, nc, kok_xadj, kok_adj);

        std::cout << "\t:" << ColoringAlgorithmNames[i] <<
            " Time:" << kkh.get_graph_coloring_handle()->get_overall_coloring_time() << " "
            "Num colors:" << kkh.get_graph_coloring_handle()->get_num_colors() << " "
            "Num Phases:" << kkh.get_graph_coloring_handle()->get_num_phases() << std::endl;
        std::cout << "\t"; KokkosKernels::Experimental::Util::print_1Dview(kkh.get_graph_coloring_handle()->get_vertex_colors());
      }
      /*
      typedef typename KernelHandle::row_index_temp_work_view_type row_view_type;
      typedef typename KernelHandle::row_index_persistent_work_view_type edge_view_type;
      row_view_type sym_row_map;
      edge_view_type sym_entries;

      KokkosKernels::Experimental::Util::symmetrize_graph_symbolic
        < const_row_index_view_type, const_nonzero_index_view_type,row_view_type, edge_view_type,ExecSpace>
        (nr, kok_xadj, kok_adj, sym_row_map,sym_entries );
      kok_xadj = row_index_view_type("s");
      kok_adj = nonzero_index_view_type("s");

      for (int i = 0; i < numColoringAlgos; ++i){
        std::cout << "\t:" << "On Temp Memory Space Running " << ColoringAlgorithmNames[i] << std::endl;
        kkh.create_graph_coloring_handle(ColoringAlgorithms[i]);
        KokkosKernels::Experimental::Graph::graph_color_symbolic
            <KernelHandle,row_view_type,edge_view_type> (&kkh,nr, nc, sym_row_map, sym_entries);

        std::cout << "\t:" << ColoringAlgorithmNames[i] <<
            " Time:" << kkh.get_graph_coloring_handle()->get_overall_coloring_time() << " "
            "Num colors:" << kkh.get_graph_coloring_handle()->get_num_colors() << " "
            "Num Phases:" << kkh.get_graph_coloring_handle()->get_num_phases() << std::endl;
        std::cout << "\t"; KokkosKernels::Experimental::Util::print_1Dview(kkh.get_graph_coloring_handle()->get_vertex_colors());
      }
      */
    }


#if defined( KOKKOS_USE_CUDA_UVM )

    //CUDA TEST3
    {
      std::cout << "Cuda TEST 3" << std::endl;
      typedef unsigned long lno_t;
      typedef float scalar_t;
      typedef int color_t;

      typedef Kokkos::Cuda ExecSpace;
      typedef Kokkos::CudaUVMSpace RowMemorySpace;
      typedef Kokkos::Cuda NonzeroMemorySpace;
      ExecSpace::print_configuration(std::cout);

      typedef Kokkos::CudaUVMSpace TempWorkSpace;
      typedef Kokkos::CudaHostPinnedSpace PersistentWorkSpace;

      typedef Kokkos::View<lno_t *, Kokkos::LayoutLeft, RowMemorySpace> row_index_view_type;
      typedef Kokkos::View<lno_t *, Kokkos::LayoutRight, NonzeroMemorySpace> nonzero_index_view_type;
      typedef Kokkos::View<scalar_t *, NonzeroMemorySpace> nonzero_scalar_view_type;

      typedef Kokkos::View<color_t * , Kokkos::LayoutStride, RowMemorySpace> color_view_type;

      typedef Kokkos::View<idx *,NonzeroMemorySpace::array_layout, Kokkos::Serial, Kokkos::MemoryUnmanaged> um_array_type;
      typedef Kokkos::View<idx *,NonzeroMemorySpace::array_layout,  Kokkos::Serial, Kokkos::MemoryUnmanaged> um_edge_array_type;

      um_array_type _xadj (xadj, nr + 1);
      um_edge_array_type _adj (adj, ne);


      //row_index_view_type kok_xadj ("xadj", nr + 1);
      //nonzero_index_view_type kok_adj("adj", ne);



      typedef Kokkos::View<idx *, RowMemorySpace> tmp_view_type;
      typedef Kokkos::View<idx *, NonzeroMemorySpace> tmp_edge_view_type;
      tmp_view_type __xadj ("xadj", nr + 1);
      tmp_edge_view_type __adj("adj", ne);

      Kokkos::deep_copy (__xadj , _xadj);
      Kokkos::deep_copy (__adj , _adj);

      row_index_view_type kok_xadj ("xadj", nr + 1);
      nonzero_index_view_type kok_adj("adj", ne);

      KokkosKernels::Experimental::Util::copy_vector<tmp_view_type, row_index_view_type, ExecSpace>(nr+1, __xadj, kok_xadj );
      KokkosKernels::Experimental::Util::copy_vector<tmp_edge_view_type, nonzero_index_view_type, ExecSpace>(ne, __adj, kok_adj );
      __xadj = tmp_view_type("");
      __adj = tmp_edge_view_type("");



      typedef KokkosKernels::Experimental::KokkosKernelsHandle
              <row_index_view_type,nonzero_index_view_type, nonzero_scalar_view_type, ExecSpace, TempWorkSpace,PersistentWorkSpace > KernelHandle;
      KernelHandle kkh;

      for (int i = 0; i < numColoringAlgos; ++i){
        std::cout << "\t:" << "Running " << ColoringAlgorithmNames[i] << std::endl;
        kkh.create_graph_coloring_handle(ColoringAlgorithms[i]);
        KokkosKernels::Experimental::Graph::graph_color_symbolic
              <KernelHandle,row_index_view_type,nonzero_index_view_type> (&kkh,nr, nc, kok_xadj, kok_adj);

        std::cout << "\t:" << ColoringAlgorithmNames[i] <<
            " Time:" << kkh.get_graph_coloring_handle()->get_overall_coloring_time() << " "
            "Num colors:" << kkh.get_graph_coloring_handle()->get_num_colors() << " "
            "Num Phases:" << kkh.get_graph_coloring_handle()->get_num_phases() << std::endl;
        std::cout << "\t"; KokkosKernels::Experimental::Util::print_1Dview(kkh.get_graph_coloring_handle()->get_vertex_colors());
      }
      /*
      typedef typename KernelHandle::row_index_temp_work_view_type row_view_type;
      typedef typename KernelHandle::row_index_persistent_work_view_type edge_view_type;
      row_view_type sym_row_map;
      edge_view_type sym_entries;

      KokkosKernels::Experimental::Util::symmetrize_graph_symbolic
        < row_index_view_type, nonzero_index_view_type,row_view_type, edge_view_type,ExecSpace>
        (nr, kok_xadj, kok_adj, sym_row_map,sym_entries );
      kok_xadj = row_index_view_type("s");
      kok_adj = nonzero_index_view_type("s");

      for (int i = 0; i < numColoringAlgos; ++i){
        std::cout << "\t:" << "On Temp Memory Space Running " << ColoringAlgorithmNames[i] << std::endl;
        kkh.create_graph_coloring_handle(ColoringAlgorithms[i]);
        KokkosKernels::Experimental::Graph::graph_color_symbolic
            <KernelHandle,row_view_type,edge_view_type> (&kkh,nr, nc, sym_row_map, sym_entries);

        std::cout << "\t:" << ColoringAlgorithmNames[i] <<
            " Time:" << kkh.get_graph_coloring_handle()->get_overall_coloring_time() << " "
            "Num colors:" << kkh.get_graph_coloring_handle()->get_num_colors() << " "
            "Num Phases:" << kkh.get_graph_coloring_handle()->get_num_phases() << std::endl;
        std::cout << "\t"; KokkosKernels::Experimental::Util::print_1Dview(kkh.get_graph_coloring_handle()->get_vertex_colors());
      }
      */
    }

#endif
  }
#endif
  delete [] xadj;
  delete [] adj;


  Kokkos::finalize();

  return 0;
}
