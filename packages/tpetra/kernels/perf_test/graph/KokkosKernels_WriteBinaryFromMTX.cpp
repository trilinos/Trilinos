/*
//@HEADER
// ************************************************************************
//
//          KokkosKernels: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER
*/
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
