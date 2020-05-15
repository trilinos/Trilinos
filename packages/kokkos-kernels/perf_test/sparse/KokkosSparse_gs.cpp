/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
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
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
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

#include <Kokkos_Core.hpp>
#include <Kokkos_Timer.hpp>
#include <KokkosKernels_Handle.hpp>
#include <KokkosSparse_gauss_seidel.hpp>
#include <KokkosSparse_spmv.hpp>
#include <KokkosKernels_IOUtils.hpp>
#include <KokkosBlas1_fill.hpp>
#include <KokkosBlas1_nrm2.hpp>
#include <KokkosKernels_config.h>
#include "KokkosKernels_default_types.hpp"
#include <iostream>
#include <vector>
#include <string>

using std::cout;
using std::string;

template<typename device_t>
void runGS(string matrixPath, string devName, bool symmetric, bool twostage, bool classic)
{
  typedef default_scalar scalar_t;
  typedef default_lno_t lno_t;
  typedef default_size_type size_type;
  typedef typename device_t::execution_space exec_space;
  typedef typename device_t::memory_space mem_space;
  typedef KokkosKernels::Experimental::KokkosKernelsHandle<size_type, lno_t, scalar_t, exec_space, mem_space, mem_space> KernelHandle;
  typedef typename KokkosSparse::CrsMatrix<scalar_t, lno_t, device_t, void, size_type> crsmat_t;
  //typedef typename crsmat_t::StaticCrsGraphType graph_t;
  typedef typename crsmat_t::values_type::non_const_type scalar_view_t;
  //typedef typename graph_t::row_map_type::non_const_type lno_view_t;
  //typedef typename graph_t::entries_type::non_const_type lno_nnz_view_t;
  crsmat_t A = KokkosKernels::Impl::read_kokkos_crst_matrix<crsmat_t>(matrixPath.c_str());
  lno_t nrows = A.numRows();
  lno_t ncols = A.numCols();
  if(nrows != ncols)
  {
    throw std::runtime_error("Gauss_Seidel only works for square matrices"); 
  }
  //size_type nnz = A.nnz();
  KernelHandle kh;
  //use a random RHS - uniformly distributed over (-5, 5)
  scalar_view_t b("b", nrows);
  {
    srand(54321);
    auto bhost = Kokkos::create_mirror_view(b);
    for(lno_t i = 0; i < nrows; i++)
    {
      bhost(i) = 10.0 * rand() / RAND_MAX - 5.0;
    }
    Kokkos::deep_copy(b, bhost);
  }
  kh.set_team_work_size(16);
  kh.set_dynamic_scheduling(true);
  //initial LHS is 0
  scalar_view_t x("x", nrows);
  //Data to dump to CSV
  //cluster size sequence: 1, 2, 3, 5, 8, 12, 18, ... up to 1000, or < 96 clusters (whichever comes first)
  std::vector<int> clusterSizes;
  //how long symbolic/numeric phases take (the graph reuse case isn't that interesting since numeric doesn't do much)
  std::vector<double> symbolicTimes;
  std::vector<double> numericTimes;
  std::vector<double> applyTimes;
  std::vector<double> scaledRes;
  //the initial residual norm (for x = 0) is bnorm
  double bnorm = KokkosBlas::nrm2(b);
  Kokkos::Timer timer;
  clusterSizes.push_back(1);
  clusterSizes.push_back(4);
  clusterSizes.push_back(8);
  clusterSizes.push_back(16);
  //for(int clusterSize = 1; clusterSize <= nrows / 64 && clusterSize <= 1000; clusterSize = ceil(1.5 * clusterSize))
  //  clusterSizes.push_back(clusterSize);
  for(int clusterSize : clusterSizes)
  {
    //cluster size of 1 is standard multicolor GS

    if(twostage || classic) {
      // Two-stage or Classical GS
      if (classic) {
        std::cout << "\n\n***** RUNNING CLASSICAL SGS (two-stage with inner triangular solve)\n";
      } else {
        std::cout << "\n\n***** RUNNING TWO-STAGE SGS\n";
      }
      //this constructor is for two-stage
      kh.create_gs_handle(KokkosSparse::GS_TWOSTAGE);
      kh.set_gs_twostage(!classic, nrows);
    } else if(clusterSize == 1)
    {
      std::cout << "\n\n***** RUNNING POINT COLORING SGS\n";
      //this constructor is for point coloring
      kh.create_gs_handle(KokkosSparse::GS_DEFAULT);
    }
    else
    {
      std::cout << "\n\n***** RUNNING CLUSTER SGS, cluster size = " << clusterSize << "\n";
      //this constructor is for cluster (block) coloring
      //kh.create_gs_handle(KokkosSparse::CLUSTER_CUTHILL_MCKEE, clusterSize);
      //kh.create_gs_handle(KokkosSparse::CLUSTER_DEFAULT, clusterSize);
      //kh.create_gs_handle(KokkosSparse::CLUSTER_DO_NOTHING, clusterSize);
      kh.create_gs_handle(KokkosSparse::CLUSTER_BALLOON, clusterSize);
    }
    timer.reset();
    KokkosSparse::Experimental::gauss_seidel_symbolic//<KernelHandle, lno_view_t, lno_nnz_view_t>
      (&kh, nrows, nrows, A.graph.row_map, A.graph.entries, symmetric);
    symbolicTimes.push_back(timer.seconds());
    std::cout << "\n*** symbolic time: " << symbolicTimes.back() << '\n';
    timer.reset();
    KokkosSparse::Experimental::gauss_seidel_numeric//<KernelHandle, lno_view_t, lno_nnz_view_t, scalar_view_t>
      (&kh, nrows, nrows, A.graph.row_map, A.graph.entries, A.values, symmetric);
    numericTimes.push_back(timer.seconds());
    std::cout << "\n*** numeric time: " << numericTimes.back() << '\n';
    timer.reset();
    //Last two parameters are damping factor (should be 1) and sweeps
    KokkosSparse::Experimental::symmetric_gauss_seidel_apply
      (&kh, nrows, nrows, A.graph.row_map, A.graph.entries, A.values, x, b, true, true, 1.0, 1);
    applyTimes.push_back(timer.seconds());
    std::cout << "\n*** apply time: " << applyTimes.back() << '\n';
    //Now, compute the 2-norm of residual 
    scalar_view_t res("Ax-b", nrows);
    Kokkos::deep_copy(res, b);
    typedef Kokkos::Details::ArithTraits<scalar_t> KAT;
    scalar_t alpha = KAT::one();
    scalar_t beta = -KAT::one();
    KokkosSparse::spmv<scalar_t, crsmat_t, scalar_view_t, scalar_t, scalar_view_t>
      ("N", alpha, A, x, beta, res);
    double resnorm = KokkosBlas::nrm2(res);
    //note: this still works if the solution diverges
    scaledRes.push_back(resnorm / bnorm);
    kh.destroy_gs_handle();
  }
  string csvName = "gs_perf_" + devName + ".csv";
  std::cout << "Writing results to " << csvName << "\n";
  FILE* csvDump = fopen(csvName.c_str(), "w");
  fprintf(csvDump, "ClusterSize,Symbolic,Numeric,Apply,Residual\n");
  for(size_t i = 0; i < clusterSizes.size(); i++)
  {
    fprintf(csvDump, "%d,%.5e,%.5e,%.5e,%.5e\n",
        clusterSizes[i], symbolicTimes[i], numericTimes[i], applyTimes[i], scaledRes[i]);
  }
  fclose(csvDump);
}

int main(int argc, char** argv)
{
  //Expect two args: matrix name and device flag.
  if(argc != 3 && argc != 4 && argc != 5)
  {
    std::cout << "Usage: ./sparse_gs matrix.mtx [--device] [--symmetric]\n\n";
    std::cout << "device can be \"serial\", \"openmp\", \"cuda\" or \"threads\".\n";
    std::cout << "If device is not given, the default device is used.\n";
    std::cout << "Add the --symmetric flag if the matrix is known to be symmetric.\n";
    return 0;
  }
  string device;
  string matrixPath;
  bool sym = false;
  bool twostage = false;
  bool classic = false;
  for(int i = 1; i < argc; i++)
  {
    if(!strcmp(argv[i], "--twostage"))
      twostage = true;
    else if(!strcmp(argv[i], "--classic"))
      classic = true;
    else if(!strcmp(argv[i], "--symmetric"))
      sym = true;
    else if(!strcmp(argv[i], "--serial"))
      device = "serial";
    else if(!strcmp(argv[i], "--openmp"))
      device = "openmp";
    else if(!strcmp(argv[i], "--threads"))
      device = "threads";
    else if(!strcmp(argv[i], "--cuda"))
      device = "cuda";
    else
      matrixPath = argv[i];
  }
  //No device given, so use the default one
  if(!device.length())
  {
    #ifdef KOKKOS_ENABLE_SERIAL
    if(std::is_same<Kokkos::DefaultExecutionSpace, Kokkos::Serial>::value)
      device = "serial";
    #endif
    #ifdef KOKKOS_ENABLE_OPENMP
    if(std::is_same<Kokkos::DefaultExecutionSpace, Kokkos::OpenMP>::value)
      device = "openmp";
    #endif
    #ifdef KOKKOS_ENABLE_CUDA
    if(std::is_same<Kokkos::DefaultExecutionSpace, Kokkos::Cuda>::value)
      device = "cuda";
    #endif
    #ifdef KOKKOS_ENABLE_THREADS
    if(std::is_same<Kokkos::DefaultExecutionSpace, Kokkos::Threads>::value)
      device = "threads";
    #endif
  }
  Kokkos::initialize();
  //Kokkos::ScopeGuard kokkosScope (argc, argv);

  bool run = false;
  #ifdef KOKKOS_ENABLE_SERIAL
  if(device == "serial")
  {
    runGS<Kokkos::Serial>(matrixPath, device, sym, twostage, classic);
    run = true;
  }
  #endif
  #ifdef KOKKOS_ENABLE_OPENMP
  if(device == "openmp")
  {
    runGS<Kokkos::OpenMP>(matrixPath, device, sym, twostage, classic);
    run = true;
  }
  #endif
  #ifdef KOKKOS_ENABLE_THREADS
  if(device == "threads")
  {
    runGS<Kokkos::Threads>(matrixPath, device, sym, twostage, classic);
    run = true;
  }
  #endif
  #ifdef KOKKOS_ENABLE_CUDA
  if(device == "cuda")
  {
    runGS<Kokkos::Cuda>(matrixPath, device, sym, twostage, classic);
    run = true;
  }
  #endif
  if(!run)
  {
    std::cerr << "Error: device " << device << " was requested but it's not enabled.\n";
    return 1;
  }
  Kokkos::finalize();
  return 0;
}

