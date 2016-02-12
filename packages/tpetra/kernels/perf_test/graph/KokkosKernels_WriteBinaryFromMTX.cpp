#include <KokkosKernels_GraphHelpers.hpp>
#include <cstdlib>
#include <iostream>


typedef int idx;
typedef double wt;

int main (int argc, char ** argv){
  if (argc != 3){
    std::cerr << "Usage:" << argv[0] << " matrixfile.mtx output_bin_file" << std::endl;
    exit(1);
  }

  idx nv = 0, ne = 0;
  idx *xadj, *adj;
  wt *ew;
  KokkosKernels::Experimental::Graph::Utils::read_mtx<idx,wt>(argv[1], &nv, &ne, &xadj, &adj, &ew, false, false);

  KokkosKernels::Experimental::Graph::Utils::write_graph_bin<idx, wt> (nv, ne, xadj, adj, ew, argv[2]);

  delete [] xadj;
  delete [] adj;
  delete [] ew;
}
