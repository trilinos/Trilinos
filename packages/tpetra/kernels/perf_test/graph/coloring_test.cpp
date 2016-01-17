#include <GraphColor.hpp>

#include <KokkosKernelsGraphHelpers.hpp>
#include <KokkosKernelsHandle.hpp>

#include <cstdlib>
#include <iostream>

#include <random>       // std::default_random_engine
#include <algorithm>    // std::shuffle
#include <vector>
#include "experiment_space.hpp"


int main (int argc, char ** argv){
  if (argc < 2){
    std::cerr << "Usage:" << argv[0] << " input_bin_file" << std::endl;
    exit(1);
  }

  Kokkos::initialize(argc, argv);
  MyExecSpace::print_configuration(std::cout);

  idx nr = 0, ne = 0;
  idx *xadj, *adj;
  //idx *half_srcs, *half_dsts;
  wt *ew;

  KokkosKernels::Experimental::Graph::Utils::read_graph_bin<idx, wt> (
      &nr, &ne, &xadj, &adj, &ew, argv[1]);
  delete [] ew;

  idx nc = 0;
  for (idx i = 0; i < ne; ++i){
    if (adj[i] > nc) nc = adj[i];
  }
  nc += 1;

  um_array_type _xadj (xadj, nr + 1);
  um_edge_array_type _adj (adj, ne);
  idx_array_type kok_xadj ("xadj", nr + 1);
  idx_edge_array_type kok_adj("adj", ne);
  Kokkos::deep_copy (kok_xadj, _xadj);
  Kokkos::deep_copy (kok_adj, _adj);
  delete [] xadj;
  delete [] adj;

  bool is_symmetric = true; //check symmetry.

  if (nr > nc){
    std::cerr << "Num Rows > Num Cols is not handled in graph coloring." << std::endl;
    return 1;
  }
  if (!is_symmetric || nr < nc){
    //symetrize the graph rows.
    //if num_cols > num_rows we ignore those, because we are coloring rows in this case.

    idx_array_type sym_row_map;
    idx_edge_array_type sym_entries;
    KokkosKernels::Experimental::Util::symmetrize_graph_symbolic
      < idx_array_type, idx_edge_array_type,
        idx_array_type, idx_edge_array_type,
        MyExecSpace>
      (nr, kok_xadj, kok_adj, sym_row_map,sym_entries );
    kok_xadj = sym_row_map;
    kok_adj = sym_entries;
  }


  //KokkosKernels::Experimental::Util::print_1Dview(kok_xadj);
  //KokkosKernels::Experimental::Util::print_1Dview(kok_adj);

  typedef KokkosKernels::Experimental::KokkosKernelsHandle
        <idx_array_type,idx_edge_array_type, value_array_type,
        MyExecSpace, TemporaryWorkSpace,PersistentWorkSpace > KernelHandle;

  KernelHandle kkh;

  //kkh.set_row_map(kok_xadj);
  //kkh.set_entries(kok_adj);

  //kkh.set_values();

  kkh.create_graph_coloring_handle(KokkosKernels::Experimental::Graph::COLORING_EB);

  KokkosKernels::Experimental::Graph::graph_color_symbolic<KernelHandle> (&kkh,nr, nc, kok_xadj, kok_adj);

  std::cout << "EB    " <<
      "Time:" << kkh.get_graph_coloring_handle()->get_overall_coloring_time() << " "
      "Num colors:" << kkh.get_graph_coloring_handle()->get_num_colors() << " "
      "Num Phases:" << kkh.get_graph_coloring_handle()->get_num_phases() << std::endl;
  KokkosKernels::Experimental::Util::print_1Dview(kkh.get_graph_coloring_handle()->get_vertex_colors());

  KokkosKernels::Experimental::Graph::graph_color_symbolic<KernelHandle> (&kkh,nr, nc, kok_xadj, kok_adj);


  std::cout << "EBS   " <<
      "Time:" << kkh.get_graph_coloring_handle()->get_overall_coloring_time() << " "
      "Num colors:" << kkh.get_graph_coloring_handle()->get_num_colors() << " "
      "Num Phases:" << kkh.get_graph_coloring_handle()->get_num_phases() << std::endl;

  KokkosKernels::Experimental::Util::print_1Dview(kkh.get_graph_coloring_handle()->get_vertex_colors());

  kkh.destroy_graph_coloring_handle();

  kkh.create_graph_coloring_handle(KokkosKernels::Experimental::Graph::COLORING_SERIAL);

  KokkosKernels::Experimental::Graph::graph_color_symbolic<KernelHandle> (&kkh,nr, nc, kok_xadj, kok_adj);


  std::cout << "SEQ   " <<
      "Time:" << kkh.get_graph_coloring_handle()->get_overall_coloring_time() << " "
      "Num colors:" << kkh.get_graph_coloring_handle()->get_num_colors() << " "
      "Num Phases:" << kkh.get_graph_coloring_handle()->get_num_phases() << std::endl;
  KokkosKernels::Experimental::Util::print_1Dview(kkh.get_graph_coloring_handle()->get_vertex_colors());

  kkh.destroy_graph_coloring_handle();

  kkh.create_graph_coloring_handle(KokkosKernels::Experimental::Graph::COLORING_VB);

  KokkosKernels::Experimental::Graph::graph_color_symbolic<KernelHandle> (&kkh,nr, nc, kok_xadj, kok_adj);


  std::cout << "VB    " <<
      "Time:" << kkh.get_graph_coloring_handle()->get_overall_coloring_time() << " "
      "Num colors:" << kkh.get_graph_coloring_handle()->get_num_colors() << " "
      "Num Phases:" << kkh.get_graph_coloring_handle()->get_num_phases() << std::endl;
  KokkosKernels::Experimental::Util::print_1Dview(kkh.get_graph_coloring_handle()->get_vertex_colors());
  kkh.destroy_graph_coloring_handle();


  kkh.create_graph_coloring_handle(KokkosKernels::Experimental::Graph::COLORING_VBBIT);
  //kkh.get_graph_coloring_handle()->set_tictoc(true);
  KokkosKernels::Experimental::Graph::graph_color_symbolic<KernelHandle> (&kkh,nr, nc, kok_xadj, kok_adj);


  std::cout << "VBBIT " <<
      "Time:" << kkh.get_graph_coloring_handle()->get_overall_coloring_time() << " "
      "Num colors:" << kkh.get_graph_coloring_handle()->get_num_colors() << " "
      "Num Phases:" << kkh.get_graph_coloring_handle()->get_num_phases() << std::endl;
  KokkosKernels::Experimental::Util::print_1Dview(kkh.get_graph_coloring_handle()->get_vertex_colors());
  kkh.destroy_graph_coloring_handle();

  kkh.create_graph_coloring_handle(KokkosKernels::Experimental::Graph::COLORING_VBCS);
  //kkh.get_graph_coloring_handle()->set_tictoc(true);
  KokkosKernels::Experimental::Graph::graph_color_symbolic<KernelHandle> (&kkh, nr, nc, kok_xadj, kok_adj);


  std::cout << "VBCS  " <<
      "Time:" << kkh.get_graph_coloring_handle()->get_overall_coloring_time() << " "
      "Num colors:" << kkh.get_graph_coloring_handle()->get_num_colors() << " "
      "Num Phases:" << kkh.get_graph_coloring_handle()->get_num_phases() << std::endl;
  KokkosKernels::Experimental::Util::print_1Dview(kkh.get_graph_coloring_handle()->get_vertex_colors());


  KokkosKernels::Experimental::Graph::graph_color_symbolic<KernelHandle> (&kkh,nr, nc, kok_xadj, kok_adj);


  std::cout << "VBCSS " <<
      "Time:" << kkh.get_graph_coloring_handle()->get_overall_coloring_time() << " "
      "Num colors:" << kkh.get_graph_coloring_handle()->get_num_colors() << " "
      "Num Phases:" << kkh.get_graph_coloring_handle()->get_num_phases() << std::endl;
  KokkosKernels::Experimental::Util::print_1Dview(kkh.get_graph_coloring_handle()->get_vertex_colors());

  kkh.destroy_graph_coloring_handle();

  Kokkos::finalize();

  return 0;
}
