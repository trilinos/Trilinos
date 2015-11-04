#include "pcg.hpp"
#include "experiment_space.hpp"
#include "KokkosKernelsGraphHelpers.hpp"
#include <iostream>

#define LOOPCOUNT 1//1000
#define RESCALCLOOP 1//50
#define MAXVAL 1

#define numExperimentsWithCT1 3
int colorings [] {4, 0, -1};
std::string names[] { "CUSPARSE", "EB", "CG"};
value_array_type create_x_vector(idx nv, wt max_value = 1.0){
  wt *x_vector = new wt[nv];
  for (idx i = 0; i < nv; ++i){
    wt r = static_cast <wt> (rand()) / static_cast <wt> (RAND_MAX / max_value);
    x_vector[i] = r;
  }
  um_value_array_type um_x(x_vector, nv);
  value_array_type kok_x("XVECTOR",nv);
  Kokkos::deep_copy (kok_x, um_x);
  delete [] x_vector;
  return kok_x;
}

struct SPMV{
  idx_array_type _kok_xadj;
  idx_edge_array_type _kok_adj;
  value_array_type _kok_mtx_vals;
  value_array_type _x_vector;
  value_array_type _b_vector;
  SPMV(idx_array_type kok_xadj,
      idx_edge_array_type kok_adj,
      value_array_type kok_mtx_vals,
      value_array_type x_vector,
      value_array_type b_vector):
        _kok_xadj(kok_xadj), _kok_adj(kok_adj),
        _kok_mtx_vals(kok_mtx_vals), _x_vector(x_vector), _b_vector(b_vector){}

  KOKKOS_INLINE_FUNCTION
  void operator()(const idx ii) const {
    wt result = 0;

    idx rb = _kok_xadj[ii];
    idx re = _kok_xadj[ii + 1];

    for (idx i = rb; i < re; ++i){
      idx col = _kok_adj[i];
      wt val = _kok_mtx_vals[i];
      result += val * _x_vector[col];
    }
    _b_vector[ii] = result;
  }
};

value_array_type create_b_vector(
    idx nv,
    idx_array_type kok_xadj,
    idx_edge_array_type kok_adj,
    value_array_type kok_mtx_vals,
    value_array_type x_vector){


  value_array_type b_vector ("B VECTOR", nv);
  Kokkos::parallel_for (Kokkos::RangePolicy<MyExecSpace> (0, nv) ,
      SPMV(
          kok_xadj,
          kok_adj,
          kok_mtx_vals,
          x_vector,
          b_vector));
  return b_vector;
}

void run_experiment(
    int repeatcount, idx nv, idx ne,
    idx_array_type kok_xadj,
    idx_edge_array_type kok_adj,
    value_array_type kok_mtx_vals
){

  value_array_type kok_x_original = create_x_vector(nv, MAXVAL);
  value_array_type kok_b_vector = create_b_vector(
      nv,
      kok_xadj,
      kok_adj,
      kok_mtx_vals,
      kok_x_original);

  Experimental::KokkosKernels::Example::CrsMatrix<wt, idx, MyExecSpace> A(nv ,ne, kok_xadj, kok_adj, kok_mtx_vals);


  //create X vector
  value_array_type kok_x_vector("kok_x_vector", nv);


  double gs_time = 0;
  const unsigned cg_iteration_limit = 100000;
  const double   cg_iteration_tolerance     = 1e-7 ;

  Experimental::KokkosKernels::Example::CGSolveResult cg_result ;
  typedef Experimental::KokkosKernels::KokkosKernelsHandle
        <idx_array_type,idx_edge_array_type, value_array_type,
        MyExecSpace, TemporaryWorkSpace,PersistentWorkSpace > KernelHandle;

  KernelHandle kh;
  kh.set_row_map(A.graph.row_map);
  kh.set_entries(A.graph.entries);
  kh.set_values(A.coeff);

  Kokkos::Impl::Timer timer1;
  Experimental::KokkosKernels::Example::pcgsolve(
        kh
      , A
      , kok_b_vector
      , kok_x_vector
      , cg_iteration_limit
      , cg_iteration_tolerance
      , & cg_result
      , true
  );
  Kokkos::fence();
  gs_time = timer1.seconds();


  std::cout  << " cg_iteration[" << cg_result.iteration << "]"
      << " matvec_time[" << cg_result.matvec_time << "]"
      << " cg_residual[" << cg_result.norm_res << "]"
      << " cg_iter_time[" << cg_result.iter_time << "]"
      << " precond_time[" << cg_result.precond_time << "]"
      << " precond_init_time[" << cg_result.precond_init_time << "]"
      << " precond_time/iter[" << cg_result.precond_time / (cg_result.iteration  + 1) << "]"
      << " GSTIME[" << gs_time<< "]"
      << std::endl ;



  kh.create_graph_coloring_handle(Experimental::KokkosKernels::Graph::COLORING_VB);
  Experimental::KokkosKernels::Graph::graph_color_solve<KernelHandle> (&kh);

  kok_x_vector = value_array_type("kok_x_vector", nv);
  timer1.reset();
  Experimental::KokkosKernels::Example::pcgsolve(
        kh
      , A
      , kok_b_vector
      , kok_x_vector
      , cg_iteration_limit
      , cg_iteration_tolerance
      , & cg_result
      , true
  );
  Kokkos::fence();

  gs_time = timer1.seconds();
  std::cout  << "\n\nCOLORING_VB PRECALL:\n cg_iteration[" << cg_result.iteration << "]"
      << " matvec_time[" << cg_result.matvec_time << "]"
      << " cg_residual[" << cg_result.norm_res << "]"
      << " cg_iter_time[" << cg_result.iter_time << "]"
      << " precond_time[" << cg_result.precond_time << "]"
      << " precond_init_time[" << cg_result.precond_init_time << "]"
      << " precond_time/iter[" << cg_result.precond_time / (cg_result.iteration  + 1) << "]"
      << " GSTIME[" << gs_time<< "]"
      << " numColor[" << kh.get_graph_coloring_handle()->get_num_colors()<<"]"
      << std::endl ;

  kh.destroy_graph_coloring_handle();
  kh.create_graph_coloring_handle(Experimental::KokkosKernels::Graph::COLORING_EB);
  Experimental::KokkosKernels::Graph::graph_color_solve<KernelHandle> (&kh);

  kok_x_vector = value_array_type("kok_x_vector", nv);
  timer1.reset();
  Experimental::KokkosKernels::Example::pcgsolve(
        kh
      , A
      , kok_b_vector
      , kok_x_vector
      , cg_iteration_limit
      , cg_iteration_tolerance
      , & cg_result
      , true
  );
  Kokkos::fence();

  gs_time = timer1.seconds();
  std::cout  << "\n\nCOLORING_EB PRECALL:\n cg_iteration[" << cg_result.iteration << "]"
      << " matvec_time[" << cg_result.matvec_time << "]"
      << " cg_residual[" << cg_result.norm_res << "]"
      << " cg_iter_time[" << cg_result.iter_time << "]"
      << " precond_time[" << cg_result.precond_time << "]"
      << " precond_init_time[" << cg_result.precond_init_time << "]"
      << " precond_time/iter[" << cg_result.precond_time / (cg_result.iteration  + 1) << "]"
      << " GSTIME[" << gs_time<< "]"
      << " numColor[" << kh.get_graph_coloring_handle()->get_num_colors()<<"]"
      << std::endl ;

  kh.destroy_graph_coloring_handle();
  kh.destroy_gs_handle();
  kh.create_gs_handle(Experimental::KokkosKernels::Graph::GS_PERMUTED);


  kok_x_vector = value_array_type("kok_x_vector", nv);
  timer1.reset();
  Experimental::KokkosKernels::Example::pcgsolve(
        kh
      , A
      , kok_b_vector
      , kok_x_vector
      , cg_iteration_limit
      , cg_iteration_tolerance
      , & cg_result
      , true
  );

  Kokkos::fence();
  gs_time = timer1.seconds();
  std::cout  << "\n\nPERMUTED:\n cg_iteration[" << cg_result.iteration << "]"
      << " matvec_time[" << cg_result.matvec_time << "]"
      << " cg_residual[" << cg_result.norm_res << "]"
      << " cg_iter_time[" << cg_result.iter_time << "]"
      << " precond_time[" << cg_result.precond_time << "]"
      << " precond_init_time[" << cg_result.precond_init_time << "]"
      << " precond_time/iter[" << cg_result.precond_time / (cg_result.iteration  + 1) << "]"
      << " GSTIME[" << gs_time<< "]"
      << std::endl ;


  kh.destroy_graph_coloring_handle();
  kh.destroy_gs_handle();
  kh.create_gs_handle(Experimental::KokkosKernels::Graph::GS_TEAM);

  kok_x_vector = value_array_type("kok_x_vector", nv);
  timer1.reset();
  Experimental::KokkosKernels::Example::pcgsolve(
        kh
      , A
      , kok_b_vector
      , kok_x_vector
      , cg_iteration_limit
      , cg_iteration_tolerance
      , & cg_result
      , true
  );
  Kokkos::fence();

  gs_time = timer1.seconds();
  std::cout  << "\n\nGSTEAM:\n cg_iteration[" << cg_result.iteration << "]"
      << " matvec_time[" << cg_result.matvec_time << "]"
      << " cg_residual[" << cg_result.norm_res << "]"
      << " cg_iter_time[" << cg_result.iter_time << "]"
      << " precond_time[" << cg_result.precond_time << "]"
      << " precond_init_time[" << cg_result.precond_init_time << "]"
      << " precond_time/iter[" << cg_result.precond_time / (cg_result.iteration  + 1) << "]"
      << " GSTIME[" << gs_time<< "]"
      << std::endl ;




}





int main (int argc, char ** argv){
  if (argc < 2){
    std::cerr << "Usage:" << argv[0] << " input_bin_file" << std::endl;
    exit(1);
  }


  Kokkos::initialize(argc, argv);
  MyExecSpace::print_configuration(std::cout);
  idx nv = 0, ne = 0;
  idx *xadj, *adj, *half_srcs, *half_dsts;
  wt *ew;

  Experimental::KokkosKernels::Graph::Utils::read_graph_bin<idx, wt> (&nv, &ne, &xadj, &adj, &ew, argv[1]);

  Experimental::KokkosKernels::Graph::Utils::md_malloc<idx>(&half_srcs, ne/2);
  Experimental::KokkosKernels::Graph::Utils::md_malloc<idx>(&half_dsts, ne/2);
  Experimental::KokkosKernels::Graph::Utils::convert_crs_to_lower_triangle_edge_list<idx>(nv, xadj, adj, half_srcs, half_dsts);

  std::cout << "nv:" << nv << " ne:" << ne << std::endl;

  um_array_type _xadj (xadj, nv + 1);
  um_edge_array_type _adj (adj, ne);

  um_edge_array_type _half_srcs (half_srcs, ne/2);
  um_edge_array_type _half_dests (half_dsts, ne/2);



  idx_array_type kok_xadj ("xadj", nv + 1);
  idx_edge_array_type kok_adj("adj", ne);

  Kokkos::deep_copy (kok_xadj, _xadj);
  Kokkos::deep_copy (kok_adj, _adj);

  idx_edge_array_type kok_half_srcs("adj", ne/2), kok_half_dsts("adj", ne/2);
  Kokkos::deep_copy (kok_half_srcs, _half_srcs);
  Kokkos::deep_copy (kok_half_dsts, _half_dests);


  wt_um_edge_array_type _mtx_vals (ew, ne);
  value_array_type kok_mtx_vals ("MTX_VALS", ne);
  Kokkos::deep_copy (kok_mtx_vals, _mtx_vals);

  delete [] xadj;
  delete [] adj;
  delete [] half_srcs;
  delete [] half_dsts;
  delete [] ew;

  //fill_experiments(nv, ne, kok_xadj, kok_adj, kok_half_srcs, kok_half_dsts);
  run_experiment(REPEAT, nv, ne, kok_xadj, kok_adj, kok_mtx_vals);
  //free_experiments();

  Kokkos::finalize();

  return 0;
}
