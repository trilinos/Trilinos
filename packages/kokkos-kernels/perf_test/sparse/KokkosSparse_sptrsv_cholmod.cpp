/*
//@HEADER
// ************************************************************************
//               KokkosKernels 0.9: Linear Algebra and Graph Kernels
//                 Copyright 2017 Sandia Corporation
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
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#include "Kokkos_Random.hpp"
#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosSparse_spmv.hpp"
#include "KokkosSparse_sptrsv.hpp"

#if defined( KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA )         && \
  (!defined(KOKKOS_ENABLE_CUDA) || (8000 <= CUDA_VERSION)) && \
    defined(KOKKOSKERNELS_INST_DOUBLE)

#if defined(KOKKOSKERNELS_ENABLE_TPL_CHOLMOD) && \
    defined(KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV)

#include "cholmod.h"
// auxiliary functions in perf_test (e.g., pivoting, printing)
#include "KokkosSparse_sptrsv_aux.hpp"

using namespace KokkosKernels;
using namespace KokkosKernels::Experimental;
using namespace KokkosSparse;
using namespace KokkosSparse::Experimental;
using namespace KokkosSparse::PerfTest::Experimental;

enum {CUSPARSE, SUPERNODAL_NAIVE, SUPERNODAL_ETREE};


/* ========================================================================================= */
template<typename scalar_type>
void print_factor_cholmod(cholmod_factor *L, cholmod_common *cm) {

  scalar_type *Lx;
  int *mb, *colptr, *rowind, *nb;
  int nsuper, j1, j2, i1, i2, psx, nsrow, nscol, i, ii, jj, s,
      nsrow2, ps2;

  /* ---------------------------------------------------------------------- */
  /* get inputs */
  /* ---------------------------------------------------------------------- */

  nsuper = L->nsuper;      // # of supernodal columns
  mb = (int*)(L->pi);      // mb[s+1] - mb[s] = total number of rows in all the s-th supernodes (diagonal+off-diagonal)
  nb = (int*)(L->super);
  colptr = (int*)(L->px);
  rowind = (int*)(L->s);               // rowind
  Lx = (scalar_type*)(L->x);                // data

  printf( " >> print factor(n=%ld, nsuper=%d) <<\n",L->n,nsuper );
  for (s = 0 ; s < nsuper ; s++) {
    j1 = nb [s];
    j2 = nb [s+1];
    nscol = j2 - j1;      // number of columns in the s-th supernode column
    printf( " nb[%d] = %d\n",s,nscol );
  }
  for (s = 0 ; s < nsuper ; s++) {
    j1 = nb [s];
    j2 = nb [s+1];
    nscol = j2 - j1;      // number of columns in the s-th supernode column

    i1 = mb [s];
    i2 = mb [s+1];
    nsrow  = i2 - i1;    // "total" number of rows in all the supernodes (diagonal+off-diagonal)
    nsrow2 = nsrow - nscol;  // "total" number of rows in all the off-diagonal supernodes
    ps2    = i1 + nscol;    // offset into rowind

    psx = colptr [s];           // offset into data,   Lx[s][s]

    /* print diagonal block */
    for (ii = 0; ii < nscol; ii++) {
      for (jj = 0; jj <= ii; jj++)
        std::cout << j1+ii+1 << " " << j1+jj+1 << " " <<  Lx[psx + (ii + jj*nsrow)] << std::endl;
    }

    /* print off-diagonal blocks */
    for (ii = 0; ii < nsrow2; ii++) {
      i = rowind [ps2 + ii] ;
      for (jj = 0; jj < nscol; jj++)
        std::cout << i+1 << " " << j1+jj+1 << " " << Lx[psx + (nscol+ii + jj*nsrow)] << std::endl;
    }
  }
}


template <typename crsmat_t>
void print_factor_cholmod(crsmat_t *L) {
  using graph_t       = typename crsmat_t::StaticCrsGraphType;
  using values_view_t = typename crsmat_t::values_type::non_const_type;
  using scalar_type   = typename values_view_t::value_type;

  graph_t  graph = L->graph;
  const int      *colptr = graph.row_map.data ();
  const int      *rowind = graph.entries.data ();
  const scalar_type *Lx     = L->values.data ();

  printf( "\n -- print cholmod factor in crs (numCols = %d) --\n",L->numCols () );
  for (int j = 0; j < L->numCols (); j++) {
    for (int k = colptr[j]; k < colptr[j+1]; k++) {
      std::cout << rowind[k] << " " <<  j << " " << Lx[k] << std::endl;
    }
  }
}


/* ========================================================================================= */
template<typename scalar_type>
cholmod_factor* factor_cholmod(const int nrow, const int nnz, scalar_type *nzvals, int *rowptr, int *colind, cholmod_common *Comm, int **etree) {

  // Start Cholmod
  cholmod_common *cm = Comm;
  cholmod_start (cm);
  cm->supernodal = CHOLMOD_SUPERNODAL;

  // Manually, initialize the matrix
  cholmod_sparse A;
  A.stype = 1;   // symmetric
  A.sorted = 0;
  A.packed = 1;
  A.itype = CHOLMOD_INT;
  A.xtype = CHOLMOD_REAL;
  A.dtype = CHOLMOD_DOUBLE;

  A.nrow = nrow;
  A.ncol = nrow;
  A.nzmax = nnz;

  A.p = rowptr;
  A.x = nzvals;
  A.i = colind;

  // Symbolic factorization
  cholmod_factor *L;
  L = cholmod_analyze (&A, cm);
  if (cm->status != CHOLMOD_OK) {
    printf( " ** cholmod_analyze returned with status = %d **",cm->status );
  }

  // Numerical factorization
  if (!cholmod_factorize (&A, L, cm)) {
    printf( " ** cholmod_factorize returned FALSE **\n" );
  }
  if (cm->status != CHOLMOD_OK) {
    printf( " ** cholmod_factorize returned with status = %d, minor = %ld **",cm->status,L->minor );
    int i;
    int *Perm = (int*)(L->Perm);
    for (i = 0; i < (int)(L->n); i++) printf( "%d %d\n",i,Perm[i] );
  }
  switch (cm->selected) {
    case CHOLMOD_NATURAL: printf( "  > NATURAL ordering (%d)\n", CHOLMOD_NATURAL ); break;
    case CHOLMOD_AMD:     printf( "  > AMD ordering (%d)\n",     CHOLMOD_AMD     ); break;
    case CHOLMOD_METIS:   printf( "  > METIS ordering (%d)\n",   CHOLMOD_METIS   ); break;
    case CHOLMOD_NESDIS:  printf( "  > NESDIS ordering (%d)\n",  CHOLMOD_NESDIS  ); break;
  }
  //print_factor_cholmod<scalar_type>(L, cm);
  compute_etree_cholmod(&A, cm, etree);

  return L;
}



/* ========================================================================================= */
template<typename scalar_type>
int test_sptrsv_perf(std::vector<int> tests, std::string& filename, int loop) {

  using STS = Kokkos::Details::ArithTraits<scalar_type>;
  using ordinal_type = int;
  using size_type    = int;

  // Default spaces
  using execution_space = Kokkos::DefaultExecutionSpace;
  using memory_space = typename execution_space::memory_space;

  // Host spaces
  using host_execution_space = Kokkos::DefaultHostExecutionSpace;
  using host_memory_space = typename host_execution_space::memory_space;

  //
  using host_crsmat_t = KokkosSparse::CrsMatrix<scalar_type, ordinal_type, host_execution_space, void, size_type>;
  using      crsmat_t = KokkosSparse::CrsMatrix<scalar_type, ordinal_type,      execution_space, void, size_type>;

  //
  using graph_t = typename crsmat_t::StaticCrsGraphType;

  //
  using host_scalar_view_t = Kokkos::View<scalar_type*, host_memory_space>;
  using      scalar_view_t = Kokkos::View<scalar_type*,      memory_space>;

  //
  using KernelHandle = KokkosKernels::Experimental::KokkosKernelsHandle <size_type, ordinal_type, scalar_type,
    execution_space, memory_space, memory_space >;

  const scalar_type ZERO (0.0);
  const scalar_type ONE (1.0);

  // tolerance
  scalar_type tol = STS::epsilon();

  int num_failed = 0;
  std::cout << std::endl;
  std::cout << "Execution space: " << execution_space::name () << std::endl;
  std::cout << "Memory space   : " << memory_space::name () << std::endl;
  std::cout << std::endl;
  if (!filename.empty())
  {
    // ==============================================
    // read the matrix ** on host **
    std::cout << " CHOLMOD Tester Begin: Read matrix filename " << filename << std::endl;
    host_crsmat_t Mtx = KokkosKernels::Impl::read_kokkos_crst_matrix<host_crsmat_t>(filename.c_str()); //in_matrix
    auto  graph_host  = Mtx.graph; // in_graph
    const size_type nrows = graph_host.numRows();
    auto row_map_host = graph_host.row_map;
    auto entries_host = graph_host.entries;
    auto values_host  = Mtx.values;
    //print_factor_cholmod(&Mtx);

    Kokkos::Timer timer;
    // ==============================================
    // call CHOLMOD on the host    
    cholmod_common cm;
    cholmod_factor *L = NULL;
    int *etree;
    timer.reset();
    std::cout << " > call CHOLMOD for factorization" << std::endl;
    L = factor_cholmod<scalar_type> (nrows, Mtx.nnz(), values_host.data(), const_cast<int*> (row_map_host.data()), entries_host.data(),
                                  &cm, &etree);
    std::cout << "   Factorization Time: " << timer.seconds() << std::endl << std::endl;
    int* iperm = (int*)(L->Perm);
    int*  perm = new int[nrows];
    for (int i = 0; i < nrows; i++) {
      perm[iperm[i]] = i;
    }

    // ==============================================
    // Run all requested algorithms
    for ( auto test : tests ) {
      std::cout << "\ntest = " << test << std::endl;

      KernelHandle khL, khU;
      switch(test) {
        case SUPERNODAL_NAIVE:
        case SUPERNODAL_ETREE:
        {
          // ==============================================
          // Create handles for U and U^T solves
          if (test == SUPERNODAL_NAIVE) {
            std::cout << " > create handle for SUPERNODAL_NAIVE" << std::endl << std::endl;
            khL.create_sptrsv_handle (SPTRSVAlgorithm::SUPERNODAL_NAIVE, nrows, true);
            khU.create_sptrsv_handle (SPTRSVAlgorithm::SUPERNODAL_NAIVE, nrows, false);
          } else {
            std::cout << " > create handle for SUPERNODAL_ETREE" << std::endl << std::endl;
            khL.create_sptrsv_handle (SPTRSVAlgorithm::SUPERNODAL_ETREE, nrows, true);
            khU.create_sptrsv_handle (SPTRSVAlgorithm::SUPERNODAL_ETREE, nrows, false);
          }

          // ==============================================
          // set etree (required)
          khL.set_sptrsv_etree (etree);
          khU.set_sptrsv_etree (etree);

          // ==============================================
          // set permutation
          khL.set_sptrsv_perm (perm);
          khU.set_sptrsv_perm (perm);

          // ==============================================
          // Do symbolic analysis
          sptrsv_symbolic<scalar_type, ordinal_type, size_type> (&khL, &khU, L, &cm);

          // ==============================================
          // Do numerical compute
          sptrsv_compute<scalar_type, ordinal_type, size_type> (&khL, &khU, L, &cm);

          // ==============================================
          // Create the known solution and set to all 1's ** on host **
          host_scalar_view_t sol_host("sol_host", nrows);
          //Kokkos::deep_copy(sol_host, ONE);
          Kokkos::Random_XorShift64_Pool<host_execution_space> random(13718);
          Kokkos::fill_random(sol_host, random, scalar_type(1));

          // ==============================================
          // Create the rhs ** on host **
          // A*sol 
          host_scalar_view_t rhs_host("rhs_host", nrows);
          KokkosSparse::spmv( "N", ONE, Mtx, sol_host, ZERO, rhs_host);

          // ==============================================
          // apply forward-pivot on the host
          host_scalar_view_t tmp_host ("temp", nrows);
          forwardP_supernode<scalar_type> (nrows, perm, 1, rhs_host.data(), nrows, tmp_host.data(), nrows);

          // ==============================================
          // copy rhs to the default host/device
          scalar_view_t rhs ("rhs", nrows);
          scalar_view_t sol ("sol", nrows);
          Kokkos::deep_copy (rhs, tmp_host);

          // ==============================================
          // do L solve
         // numeric (only rhs is modified) on the default device/host space
          timer.reset();
           sptrsv_solve (&khL, sol, rhs);
          Kokkos::fence();
          std::cout << "   Solve Time   : " << timer.seconds() << std::endl;

          // ==============================================
          // do L^T solve
          // numeric (only rhs is modified) on the default device/host space
          timer.reset();
           sptrsv_solve (&khU, rhs, sol);
          Kokkos::fence ();
          std::cout << "   Solve Time   : " << timer.seconds() << std::endl;
 

          // ==============================================
          // apply backward-pivot
          // > copy solution to host
          Kokkos::deep_copy(tmp_host, rhs);
          backwardP_supernode<scalar_type>(nrows, perm, 1, tmp_host.data(), nrows, sol_host.data(), nrows);


          // ==============================================
          // Error Check ** on host **
          Kokkos::fence();
          if (!check_errors(tol, Mtx, rhs_host, sol_host)) {
            num_failed ++;
          }

          // try again?
          {
            Kokkos::deep_copy(sol_host, ONE);
            KokkosSparse::spmv( "N", ONE, Mtx, sol_host, ZERO, rhs_host);
            forwardP_supernode<scalar_type> (nrows, perm, 1, rhs_host.data(), nrows, tmp_host.data(), nrows);
            Kokkos::deep_copy (rhs, tmp_host);
             sptrsv_solve (&khL, sol, rhs);
             sptrsv_solve (&khU, rhs, sol);
            Kokkos::fence();
            Kokkos::deep_copy(tmp_host, rhs);
            backwardP_supernode<scalar_type>(nrows, perm, 1, tmp_host.data(), nrows, sol_host.data(), nrows);

            if (!check_errors(tol, Mtx, rhs_host, sol_host)) {
              num_failed ++;
            }
          }
          std::cout << std::endl;

          // Benchmark
          // L-solve
          double min_time = 1.0e32;
          double max_time = 0.0;
          double ave_time = 0.0;
          Kokkos::fence();
          for(int i=0;i<loop;i++) {
            timer.reset();
            sptrsv_solve (&khL, sol, rhs);
            Kokkos::fence();
            double time = timer.seconds();
            ave_time += time;
            if(time>max_time) max_time = time;
            if(time<min_time) min_time = time;
            //std::cout << time << std::endl;
          }
          std::cout << " L-solve: loop = " << loop << std::endl;
          std::cout << "  LOOP_AVG_TIME:  " << ave_time/loop << std::endl;
          std::cout << "  LOOP_MAX_TIME:  " << max_time << std::endl;
          std::cout << "  LOOP_MIN_TIME:  " << min_time << std::endl << std::endl;

          // U-solve
          min_time = 1.0e32;
          max_time = 0.0;
          ave_time = 0.0;
          Kokkos::fence();
          for(int i=0;i<loop;i++) {
            timer.reset();
            sptrsv_solve (&khU, rhs, sol);
            Kokkos::fence();
            double time = timer.seconds();
            ave_time += time;
            if(time>max_time) max_time = time;
            if(time<min_time) min_time = time;
            //std::cout << time << std::endl;
          }
          std::cout << " U-solve: loop = " << loop << std::endl;
          std::cout << "  LOOP_AVG_TIME:  " << ave_time/loop << std::endl;
          std::cout << "  LOOP_MAX_TIME:  " << max_time << std::endl;
          std::cout << "  LOOP_MIN_TIME:  " << min_time << std::endl << std::endl;
        }
        break;

        case CUSPARSE:
        {
          // ==============================================
          // read CHOLMOD factor int crsMatrix on the host (cholmodMat_host) and copy to default host/device (cholmodMtx)
          timer.reset();
          std::cout << " > Read Cholmod factor into KokkosSparse::CrsMatrix (invert diagonabl, and copy to device) " << std::endl;
          bool cusparse = true;
          bool invert_diag = false;
          auto graph = read_cholmod_graphL<graph_t>(cusparse, L, &cm);
          auto cholmodMtx = read_cholmod_factor<crsmat_t, graph_t> (cusparse, invert_diag, L, &cm, graph);
          std::cout << "   Conversion Time: " << timer.seconds() << std::endl << std::endl;

          bool col_majorL = true;
          bool col_majorU = false;
          if (!check_cusparse(Mtx, col_majorL, cholmodMtx, col_majorU, cholmodMtx, perm, perm, tol, loop)) {
            num_failed ++;
          }
        }
        break;

        default:
          std::cout << " > Testing only Cholmod < " << std::endl;
          exit(0);
      }
    }
  }
  std::cout << std::endl << std::endl;

  return num_failed;
}


void print_help_sptrsv() {
  printf("Options:\n");
  printf("  --test [OPTION] : Use different kernel implementations\n");
  printf("                    Options:\n");
  printf("                    cholmod_naive, cholmod_etree\n\n");
  printf("  -f [file]       : Read in Matrix Market formatted text file 'file'.\n");
  printf("  --offset [O]    : Subtract O from every index.\n");
  printf("                    Useful in case the matrix market file is not 0 based.\n\n");
  printf("  -rpt [K]        : Number of Rows assigned to a thread.\n");
  printf("  -ts [T]         : Number of threads per team.\n");
  printf("  -vl [V]         : Vector-length (i.e. how many Cuda threads are a Kokkos 'thread').\n");
  printf("  --loop [LOOP]   : How many spmv to run to aggregate average time. \n");
}


int main(int argc, char **argv)
{
  std::vector<int> tests;
  std::string filename;

  int loop = 1;

  if(argc == 1)
  {
    print_help_sptrsv();
    return 0;
  }

  for(int i=0;i<argc;i++)
  {
    if((strcmp(argv[i],"--test")==0)) {
      i++;
      if((strcmp(argv[i],"cholmod-naive")==0)) {
        tests.push_back( SUPERNODAL_NAIVE );
      }
      if((strcmp(argv[i],"cholmod-etree")==0)) {
        tests.push_back( SUPERNODAL_ETREE );
      }
      if((strcmp(argv[i],"cusparse")==0)) {
        tests.push_back( CUSPARSE );
      }
      continue;
    }
    if((strcmp(argv[i],"-f")==0)) {
      filename = argv[++i];
      continue;
    }
    if((strcmp(argv[i],"--loop")==0)) {
      loop = atoi(argv[++i]);
      continue;
    }
    if((strcmp(argv[i],"--help")==0) || (strcmp(argv[i],"-h")==0)) {
      print_help_sptrsv();
      return 0;
    }
  }

  std::cout << std::endl;
  for (size_t i = 0; i < tests.size(); ++i) {
    std::cout << "tests[" << i << "] = " << tests[i] << std::endl;
  }

  {
    // Cholmod may not support single, yet
    //int total_errors = test_sptrsv_perf<float>(tests, filename, loop);
    // Kokkos::IO may not read complex?
    //int total_errors = test_sptrsv_perf<Kokkos::complex<double>>(tests, filename, loop);
    Kokkos::ScopeGuard kokkosScope (argc, argv);
    int total_errors = test_sptrsv_perf<double>(tests, filename, loop);
    if(total_errors == 0)
      std::cout << "Kokkos::SPTRSV Test: Passed" << std::endl << std::endl;
    else
      std::cout << "Kokkos::SPTRSV Test: Failed" << std::endl << std::endl;
  }

  return 0;
}
#else // defined(KOKKOSKERNELS_ENABLE_TPL_CHOLMOD)
int main(int argc, char **argv)
{
  std::cout << std::endl << "** CHOLMOD NOT ENABLED **" << std::endl << std::endl;
  return 0;
}
#endif

#else // defined( KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA ) && (!defined(KOKKOS_ENABLE_CUDA) || ( 8000 <= CUDA_VERSION ))
int main() {
#if !defined(KOKKOSKERNELS_INST_DOUBLE)
  std::cout << " Only supported with double precision" << std::endl;
#endif
#if !defined( KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA )
  std::cout << " KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA **not** defined" << std::endl;
#endif
#if defined(KOKKOS_ENABLE_CUDA)
  std::cout << " KOKKOS_ENABLE_CUDA defined" << std::endl;
  #if !defined(KOKKOS_ENABLE_CUDA_LAMBDA)
  std::cout << " KOKKOS_ENABLE_CUDA_LAMBDA not defined\n" << std::endl;
  #endif
  std::cout << " CUDA_VERSION = " << CUDA_VERSION << std::endl;
#endif
  return 0;
}
#endif
