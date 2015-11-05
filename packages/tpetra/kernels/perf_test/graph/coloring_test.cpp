
#include <KokkosKernelsGraphHelpers.hpp>
#include <KokkosKernelsHandle.hpp>
#include <GraphColor.hpp>
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

  idx nv = 0, ne = 0;
  idx *xadj, *adj;
  //idx *half_srcs, *half_dsts;
  wt *ew;

  Experimental::KokkosKernels::Graph::Utils::read_graph_bin<idx, wt> (
      &nv, &ne, &xadj, &adj, &ew, argv[1]);
  delete [] ew;


  um_array_type _xadj (xadj, nv + 1);
  um_edge_array_type _adj (adj, ne);
  idx_array_type kok_xadj ("xadj", nv + 1);
  idx_edge_array_type kok_adj("adj", ne);
  Kokkos::deep_copy (kok_xadj, _xadj);
  Kokkos::deep_copy (kok_adj, _adj);
  delete [] xadj;
  delete [] adj;


  typedef Experimental::KokkosKernels::KokkosKernelsHandle
        <idx_array_type,idx_edge_array_type, value_array_type,
        MyExecSpace, TemporaryWorkSpace,PersistentWorkSpace > KernelHandle;

  KernelHandle kkh;

  kkh.set_row_map(kok_xadj);
  kkh.set_entries(kok_adj);
  //kkh.set_values();

  kkh.create_graph_coloring_handle(Experimental::KokkosKernels::Graph::COLORING_EB);

  Experimental::KokkosKernels::Graph::graph_color_solve<KernelHandle> (&kkh);

  std::cout << "EB    " <<
      "Time:" << kkh.get_graph_coloring_handle()->get_overall_coloring_time() << " "
      "Num colors:" << kkh.get_graph_coloring_handle()->get_num_colors() << " "
      "Num Phases:" << kkh.get_graph_coloring_handle()->get_num_phases() << std::endl;

  Experimental::KokkosKernels::Graph::graph_color_solve<KernelHandle> (&kkh);

  std::cout << "EBS   " <<
      "Time:" << kkh.get_graph_coloring_handle()->get_overall_coloring_time() << " "
      "Num colors:" << kkh.get_graph_coloring_handle()->get_num_colors() << " "
      "Num Phases:" << kkh.get_graph_coloring_handle()->get_num_phases() << std::endl;

  kkh.destroy_graph_coloring_handle();

  kkh.create_graph_coloring_handle(Experimental::KokkosKernels::Graph::COLORING_SERIAL);

  Experimental::KokkosKernels::Graph::graph_color_solve<KernelHandle> (&kkh);


  std::cout << "SEQ   " <<
      "Time:" << kkh.get_graph_coloring_handle()->get_overall_coloring_time() << " "
      "Num colors:" << kkh.get_graph_coloring_handle()->get_num_colors() << " "
      "Num Phases:" << kkh.get_graph_coloring_handle()->get_num_phases() << std::endl;

  kkh.destroy_graph_coloring_handle();

  kkh.create_graph_coloring_handle(Experimental::KokkosKernels::Graph::COLORING_VB);

  Experimental::KokkosKernels::Graph::graph_color_solve<KernelHandle> (&kkh);


  std::cout << "VB    " <<
      "Time:" << kkh.get_graph_coloring_handle()->get_overall_coloring_time() << " "
      "Num colors:" << kkh.get_graph_coloring_handle()->get_num_colors() << " "
      "Num Phases:" << kkh.get_graph_coloring_handle()->get_num_phases() << std::endl;
  kkh.destroy_graph_coloring_handle();


  kkh.create_graph_coloring_handle(Experimental::KokkosKernels::Graph::COLORING_VBBIT);
  //kkh.get_graph_coloring_handle()->set_tictoc(true);
  Experimental::KokkosKernels::Graph::graph_color_solve<KernelHandle> (&kkh);


  std::cout << "VBBIT " <<
      "Time:" << kkh.get_graph_coloring_handle()->get_overall_coloring_time() << " "
      "Num colors:" << kkh.get_graph_coloring_handle()->get_num_colors() << " "
      "Num Phases:" << kkh.get_graph_coloring_handle()->get_num_phases() << std::endl;
  kkh.destroy_graph_coloring_handle();

  kkh.create_graph_coloring_handle(Experimental::KokkosKernels::Graph::COLORING_VBCS);
  //kkh.get_graph_coloring_handle()->set_tictoc(true);
  Experimental::KokkosKernels::Graph::graph_color_solve<KernelHandle> (&kkh);


  std::cout << "VBCS  " <<
      "Time:" << kkh.get_graph_coloring_handle()->get_overall_coloring_time() << " "
      "Num colors:" << kkh.get_graph_coloring_handle()->get_num_colors() << " "
      "Num Phases:" << kkh.get_graph_coloring_handle()->get_num_phases() << std::endl;


  Experimental::KokkosKernels::Graph::graph_color_solve<KernelHandle> (&kkh);


  std::cout << "VBCSS " <<
      "Time:" << kkh.get_graph_coloring_handle()->get_overall_coloring_time() << " "
      "Num colors:" << kkh.get_graph_coloring_handle()->get_num_colors() << " "
      "Num Phases:" << kkh.get_graph_coloring_handle()->get_num_phases() << std::endl;

  kkh.destroy_graph_coloring_handle();

  Kokkos::finalize();

  return 0;
}
