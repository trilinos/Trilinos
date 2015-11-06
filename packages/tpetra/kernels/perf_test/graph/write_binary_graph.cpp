#include <KokkosKernelsGraphHelpers.hpp>
#include <cstdlib>
#include <iostream>

#include "experiment_space.hpp"

int main (int argc, char ** argv){
  if (argc != 3){
    std::cerr << "Usage:" << argv[0] << " matrixfile.mtx output_bin_file" << std::endl;
    exit(1);
  }

  idx nv = 0, ne = 0;
  idx *xadj, *adj;
  wt *ew;
  Experimental::KokkosKernels::Graph::Utils::read_mtx<idx,wt>(argv[1], &nv, &ne, &xadj, &adj, &ew, true, false);

  Experimental::KokkosKernels::Graph::Utils::write_graph_bin<idx, wt> (nv, ne, xadj, adj, ew, argv[2]);

  delete [] xadj;
  delete [] adj;
  delete [] ew;
}
