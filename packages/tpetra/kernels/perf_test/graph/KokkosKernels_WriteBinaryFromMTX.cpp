#include <KokkosKernels_GraphHelpers.hpp>
#include <cstdlib>
#include <iostream>


typedef int idx;
typedef double wt;

int main (int argc, char ** argv){


  bool symmetrize = false, remove_diagonal = false, transpose = false;
  char *in_mtx = NULL, *out_bin = NULL;
  for ( int i = 1 ; i < argc ; ++i ) {
    if ( 0 == strcasecmp( argv[i] , "symmetrize" ) ) {
      symmetrize = true;
    }
    else if ( 0 == strcasecmp( argv[i] , "remove_diagonal" ) ) {
      remove_diagonal = true;
    }
    else if ( 0 == strcasecmp( argv[i] , "transpose" ) ) {
      transpose = true;
    }
    else if ( 0 == strcasecmp( argv[i] , "in_mtx" ) ) {
      in_mtx = argv[++i];
    }
    else if ( 0 == strcasecmp( argv[i] , "out_bin" ) ) {
      out_bin = argv[++i];
    }
    else {
      std::cerr << "Usage:" << argv[0]
                << " in_mtx matrixfile.mtx out_bin output_bin_file [symmetrize] [remove_diagonal] [transpose]" << std::endl;
      exit(1);
    }
  }
  if (in_mtx == NULL || out_bin == NULL){
    std::cerr << "Usage:" << argv[0]
              << " in_mtx matrixfile.mtx out_bin output_bin_file [symmetrize] [remove_diagonal] [transpose]" << std::endl;
    exit(1);
  }

  idx nv = 0, ne = 0;
  idx *xadj, *adj;
  wt *ew;
  KokkosKernels::Experimental::Graph::Utils::read_mtx<idx,wt>
      (in_mtx, &nv, &ne, &xadj, &adj, &ew, symmetrize, remove_diagonal, transpose);

  KokkosKernels::Experimental::Graph::Utils::write_graph_bin<idx, wt> (nv, ne, xadj, adj, ew, out_bin);

  delete [] xadj;
  delete [] adj;
  delete [] ew;
}
