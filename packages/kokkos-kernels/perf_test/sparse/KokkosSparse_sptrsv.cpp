/*
//@HEADER
// ************************************************************************
//
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

#include <cstdio>

#include <ctime>
#include <cstring>
#include <cstdlib>
#include <limits>
#include <limits.h>
#include <cmath>
#include <unordered_map>

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
#include <cusparse.h>
#endif

#include <Kokkos_Core.hpp>
#include <matrix_market.hpp>

#include "KokkosKernels_SparseUtils.hpp"
#include "KokkosSparse_sptrsv.hpp"
#include "KokkosSparse_spmv.hpp"
#include "KokkosSparse_CrsMatrix.hpp"
#include <KokkosKernels_IOUtils.hpp>

#if defined( KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA ) && (!defined(KOKKOS_ENABLE_CUDA) || ( 8000 <= CUDA_VERSION ))
using namespace KokkosSparse;
using namespace KokkosSparse::Experimental;
using namespace KokkosKernels;
using namespace KokkosKernels::Experimental;

enum {DEFAULT, CUSPARSE, LVLSCHED_RP, LVLSCHED_TP1/*, LVLSCHED_TP2*/};



template<typename Scalar>
int test_sptrsv_perf(std::vector<int> tests, std::string& lfilename, std::string& ufilename, int team_size, int vector_length, int idx_offset, int loop) {

  typedef Scalar scalar_t;
  typedef int lno_t;
  typedef int size_type;
  typedef Kokkos::DefaultExecutionSpace execution_space;
  typedef typename execution_space::memory_space memory_space;

  typedef KokkosSparse::CrsMatrix<scalar_t, lno_t, execution_space, void, size_type> crsmat_t;
  typedef typename crsmat_t::StaticCrsGraphType graph_t;

  typedef Kokkos::View< scalar_t*, memory_space >     ValuesType;

  typedef KokkosKernels::Experimental::KokkosKernelsHandle <size_type, lno_t, scalar_t,
    execution_space, memory_space, memory_space > KernelHandle;

  scalar_t ZERO = scalar_t(0);
  scalar_t ONE = scalar_t(1);


// Read lmtx
// Run all requested algorithms
// Read umtx
// Run all requested algorithms

// LOWERTRI
  std::cout << "\n\n" << std::endl;
  if (!lfilename.empty())
  {
    std::cout << "Lower Tri Begin: Read matrix filename " << lfilename << std::endl;
    crsmat_t triMtx = KokkosKernels::Impl::read_kokkos_crst_matrix<crsmat_t>(lfilename.c_str()); //in_matrix
    graph_t  graph  = triMtx.graph; // in_graph
    const size_type nrows = graph.numRows();

    // Create the rhs and lhs_known solution
    ValuesType known_lhs("known_lhs", nrows);
    // Create known solution lhs set to all 1's
    Kokkos::deep_copy(known_lhs, ONE);

    // Solution to find
    ValuesType lhs("lhs", nrows);

    // A*known_lhs generates rhs: rhs is dense, use spmv
    ValuesType rhs("rhs", nrows);

    std::cout << "SPMV" << std::endl;
    KokkosSparse::spmv( "N", ONE, triMtx, known_lhs, ZERO, rhs);


    auto row_map = graph.row_map;
    auto entries = graph.entries;
    auto values  = triMtx.values;



#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
  //std::cout << "  cusparse: create handle" << std::endl;
  cusparseStatus_t status;
  cusparseHandle_t handle = 0;
  status = cusparseCreate(&handle);
  if (CUSPARSE_STATUS_SUCCESS != status)
    std::cout << "handle create status error name " << (status) << std::endl;
  cusparseSetPointerMode(handle, CUSPARSE_POINTER_MODE_HOST);
  cusparseMatDescr_t descr = 0;
  csrsv2Info_t info = 0;
  int pBufferSize;
  void *pBuffer = 0;
  int structural_zero;
  int numerical_zero;
  const double alpha = 1.;
  const cusparseSolvePolicy_t policy = CUSPARSE_SOLVE_POLICY_USE_LEVEL;
  const cusparseOperation_t trans = CUSPARSE_OPERATION_NON_TRANSPOSE;
  
  // step 1: create a descriptor which contains
  // - matrix L is lower triangular
  //   (L may not have all diagonal elements.)
  status = cusparseCreateMatDescr(&descr);
  if (CUSPARSE_STATUS_SUCCESS != status)
    std::cout << "matdescr create status error name " << (status) << std::endl;
  //cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ONE);
  cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO);
  cusparseSetMatFillMode(descr, CUSPARSE_FILL_MODE_LOWER);
  cusparseSetMatDiagType(descr, CUSPARSE_DIAG_TYPE_NON_UNIT);
  cusparseSetMatType(descr,CUSPARSE_MATRIX_TYPE_GENERAL);
  //cusparseSetMatDiagType(descr, CUSPARSE_DIAG_TYPE_UNIT);
  
  // step 2: create a empty info structure
  //std::cout << "  cusparse: create csrsv2info" << std::endl;
  status = cusparseCreateCsrsv2Info(&info);
  if (CUSPARSE_STATUS_SUCCESS != status)
    std::cout << "csrsv2info create status error name " << (status) << std::endl;
  
  // step 3: query how much memory used in csrsv2, and allocate the buffer
        int nnz = triMtx.nnz();
  cusparseDcsrsv2_bufferSize(handle, trans, nrows, nnz, descr,
      values.data(), row_map.data(), entries.data(), info, &pBufferSize);
  // pBuffer returned by cudaMalloc is automatically aligned to 128 bytes.
  cudaMalloc((void**)&pBuffer, pBufferSize);
#endif


  for ( auto test : tests ) {
    std::cout << "\ntest = " << test << std::endl;

    KernelHandle kh;
    bool is_lower_tri = true;

    std::cout << "Create handle" << std::endl;
    switch(test) {
      case LVLSCHED_RP:
        kh.create_sptrsv_handle(SPTRSVAlgorithm::SEQLVLSCHD_RP, nrows, is_lower_tri);
        kh.get_sptrsv_handle()->print_algorithm();
        break;
      case LVLSCHED_TP1:
        kh.create_sptrsv_handle(SPTRSVAlgorithm::SEQLVLSCHD_TP1, nrows, is_lower_tri);
        kh.get_sptrsv_handle()->print_algorithm();
        break;
/*
      case LVLSCHED_TP2:
        kh.create_sptrsv_handle(SPTRSVAlgorithm::SEQLVLSCHED_TP2, nrows, is_lower_tri);
        kh.get_sptrsv_handle()->print_algorithm();
        break;
*/
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
      case CUSPARSE:
        std::cout << "CUSPARSE: No kk interface added yet" << std::endl;
        //cusparse_matvec(A, x, y, rows_per_thread, team_size, vector_length);
        break;
#endif
      default:
        kh.create_sptrsv_handle(SPTRSVAlgorithm::SEQLVLSCHD_TP1, nrows, is_lower_tri);
        kh.get_sptrsv_handle()->print_algorithm();
    }


    // Init run to clear cache etc.
    Kokkos::Timer timer;
    if (test != CUSPARSE) {
    timer.reset();
    sptrsv_symbolic( &kh, row_map, entries );
    std::cout << "LTRI Symbolic Time: " << timer.seconds() << std::endl;

    //std::cout << "TriSolve Solve" << std::endl;
    timer.reset();
    sptrsv_solve( &kh, row_map, entries, values, rhs, lhs );
    Kokkos::fence();
    std::cout << "LTRI Solve Time: " << timer.seconds() << std::endl;
  
    }
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
// step 4: perform analysis
    else {
#if 0
      double *dvalues = (double *)(values.data());
      int *drow_map = (int *)(row_map.data());
      int *dentries = (int *)(entries.data());
      double *dlhs = (double *)(lhs.data());
      double *drhs = (double *)(rhs.data());
#endif
      //int nnz = triMtx.nnz();
      //std::cout << "  cusparse path: analysis" << std::endl;
      //status = cusparseDcsrsv2_analysis(handle, trans, nrows, nnz, descr, (double*)dvalues, (int *)drow_map, (int *)dentries, info, policy, pBuffer);
      timer.reset();
      status = cusparseDcsrsv2_analysis(handle, trans, nrows, triMtx.nnz(), descr, values.data(), row_map.data(), entries.data(), info, policy, pBuffer);
      std::cout << "LTRI Cusparse Symbolic Time: " << timer.seconds() << std::endl;
      if (CUSPARSE_STATUS_SUCCESS != status)
        std::cout << "analysis status error name " << (status) << std::endl;
// L has unit diagonal, so no structural zero is reported.

      //std::cout << "  cusparse path: analysis" << std::endl;
      status = cusparseXcsrsv2_zeroPivot(handle, info, &structural_zero);
      if (CUSPARSE_STATUS_ZERO_PIVOT == status){
         printf("L(%d,%d) is missing\n", structural_zero, structural_zero);
      }

// step 5: solve L*y = x
      //std::cout << "  cusparse path: solve" << std::endl;
      //status = cusparseDcsrsv2_solve(handle, trans, nrows, nnz, &alpha, descr, (double*)dvalues, (int *)drow_map, (int *)dentries, info, (double*)drhs, (double*)dlhs, policy, pBuffer);
      timer.reset();
      status = cusparseDcsrsv2_solve(handle, trans, nrows, triMtx.nnz(), &alpha, descr, values.data(), row_map.data(), entries.data(), info, rhs.data(), lhs.data(), policy, pBuffer);
      Kokkos::fence();
      std::cout << "LTRI Cusparse Solve Time: " << timer.seconds() << std::endl;
      if (CUSPARSE_STATUS_SUCCESS != status)
        std::cout << "solve status error name " << (status) << std::endl;
// L has unit diagonal, so no numerical zero is reported.
      status = cusparseXcsrsv2_zeroPivot(handle, info, &numerical_zero);
      if (CUSPARSE_STATUS_ZERO_PIVOT == status){
         printf("L(%d,%d) is zero\n", numerical_zero, numerical_zero);
      }
    }
#endif
    // Error Check
    scalar_t sum = 0.0;
    Kokkos::fence();
    Kokkos::parallel_reduce( Kokkos::RangePolicy<execution_space>(0, lhs.extent(0)), 
      KOKKOS_LAMBDA ( const lno_t i, scalar_t &tsum ) {
        tsum += lhs(i);
      }, sum);
  
    if ( sum != lhs.extent(0) ) {
      std::cout << "Lower Tri Solve FAILURE: sum = " << sum << std::endl;
      auto hsoln = Kokkos::create_mirror_view(lhs);
      Kokkos::deep_copy(hsoln, lhs);
      for ( size_t i = 0; i < hsoln.extent(0); ++i ) {
        std::cout << "lhs(" << i << ") = " << hsoln(i) << std::endl;
      }
      return 1;
    }
    else {
     std::cout << "Lower Tri Solve SUCCESS!" << std::endl;
    }

  
    // Benchmark
    Kokkos::fence();
    double min_time = 1.0e32;
    double max_time = 0.0;
    double ave_time = 0.0;

    for(int i=0;i<loop;i++) {
      timer.reset();
  
    if (test != CUSPARSE) {
      sptrsv_solve( &kh, row_map, entries, values, rhs, lhs );
    }
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
    else {
      cusparseDcsrsv2_solve(handle, trans, nrows, triMtx.nnz(), &alpha, descr, values.data(), row_map.data(), entries.data(), info, rhs.data(), lhs.data(), policy, pBuffer);
    }
#endif
  
      Kokkos::fence();
      double time = timer.seconds();
      ave_time += time;
      if(time>max_time) max_time = time;
      if(time<min_time) min_time = time;
    }
    std::cout << "LOOP_AVG_TIME:  " << ave_time/loop << std::endl;
    std::cout << "LOOP_MAX_TIME:  " << max_time << std::endl;
    std::cout << "LOOP_MIN_TIME:  " << min_time << std::endl;

  }
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
// step 6: free resources
    cudaFree(pBuffer);
    cusparseDestroyCsrsv2Info(info);
    cusparseDestroyMatDescr(descr);
    cusparseDestroy(handle);
#endif
  }

  std::cout << "\n\n" << std::endl;
// UPPERTRI
  if (!ufilename.empty())
  {
    std::cout << "Upper Tri Begin: Read matrix filename " << ufilename << std::endl;
    crsmat_t triMtx = KokkosKernels::Impl::read_kokkos_crst_matrix<crsmat_t>(ufilename.c_str()); //in_matrix
    graph_t  graph  = triMtx.graph; // in_graph
    const size_type nrows = graph.numRows();

    // Create the rhs and lhs_known solution
    ValuesType known_lhs("known_lhs", nrows);
    // Create known solution lhs set to all 1's
    Kokkos::deep_copy(known_lhs, ONE);

    // Solution to find
    ValuesType lhs("lhs", nrows);

    // A*known_lhs generates rhs: rhs is dense, use spmv
    ValuesType rhs("rhs", nrows);

    std::cout << "SPMV" << std::endl;
    KokkosSparse::spmv( "N", ONE, triMtx, known_lhs, ZERO, rhs);

    auto row_map = graph.row_map;
    auto entries = graph.entries;
    auto values  = triMtx.values;

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
  //std::cout << "  cusparse: create handle" << std::endl;
  cusparseStatus_t status;
  cusparseHandle_t handle = 0;
  status = cusparseCreate(&handle);
  if (CUSPARSE_STATUS_SUCCESS != status)
    std::cout << "handle create status error name " << (status) << std::endl;
  cusparseSetPointerMode(handle, CUSPARSE_POINTER_MODE_HOST);
  cusparseMatDescr_t descr = 0;
  csrsv2Info_t info = 0;
  int pBufferSize;
  void *pBuffer = 0;
  int structural_zero;
  int numerical_zero;
  const double alpha = 1.;
  const cusparseSolvePolicy_t policy = CUSPARSE_SOLVE_POLICY_USE_LEVEL;
  const cusparseOperation_t trans = CUSPARSE_OPERATION_NON_TRANSPOSE;
  
  // step 1: create a descriptor which contains
  //   (L may not have all diagonal elements.)
  status = cusparseCreateMatDescr(&descr);
  if (CUSPARSE_STATUS_SUCCESS != status)
    std::cout << "matdescr create status error name " << (status) << std::endl;
  //cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ONE);
  cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO);
  cusparseSetMatFillMode(descr, CUSPARSE_FILL_MODE_UPPER);
  cusparseSetMatDiagType(descr, CUSPARSE_DIAG_TYPE_NON_UNIT);
  cusparseSetMatType(descr,CUSPARSE_MATRIX_TYPE_GENERAL);
  //cusparseSetMatType(descr,CUSPARSE_MATRIX_TYPE_TRIANGULAR);
  //cusparseSetMatDiagType(descr, CUSPARSE_DIAG_TYPE_UNIT);
  
  // step 2: create a empty info structure
  //std::cout << "  cusparse: create csrsv2info" << std::endl;
  status = cusparseCreateCsrsv2Info(&info);
  if (CUSPARSE_STATUS_SUCCESS != status)
    std::cout << "csrsv2info create status error name " << (status) << std::endl;
  
  // step 3: query how much memory used in csrsv2, and allocate the buffer
        int nnz = triMtx.nnz();
  cusparseDcsrsv2_bufferSize(handle, trans, nrows, nnz, descr,
      values.data(), row_map.data(), entries.data(), info, &pBufferSize);
  // pBuffer returned by cudaMalloc is automatically aligned to 128 bytes.
  cudaMalloc((void**)&pBuffer, pBufferSize);
#endif


  for ( auto test : tests ) {
    std::cout << "\ntest = " << test << std::endl;

    KernelHandle kh;
    bool is_lower_tri = false;

    std::cout << "Create handle" << std::endl;
    switch(test) {
      case LVLSCHED_RP:
        kh.create_sptrsv_handle(SPTRSVAlgorithm::SEQLVLSCHD_RP, nrows, is_lower_tri);
        kh.get_sptrsv_handle()->print_algorithm();
        break;
      case LVLSCHED_TP1:
        kh.create_sptrsv_handle(SPTRSVAlgorithm::SEQLVLSCHD_TP1, nrows, is_lower_tri);
        kh.get_sptrsv_handle()->print_algorithm();
        break;
/*
      case LVLSCHED_TP2:
        kh.create_sptrsv_handle(SPTRSVAlgorithm::SEQLVLSCHED_TP2, nrows, is_lower_tri);
        kh.get_sptrsv_handle()->print_algorithm();
        break;
*/
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
      case CUSPARSE:
        std::cout << "CUSPARSE: No kk interface added yet" << std::endl;
        //cusparse_matvec(A, x, y, rows_per_thread, team_size, vector_length);
        break;
#endif
      default:
        kh.create_sptrsv_handle(SPTRSVAlgorithm::SEQLVLSCHD_TP1, nrows, is_lower_tri);
        kh.get_sptrsv_handle()->print_algorithm();
    }


    // Init run to clear cache etc.
    Kokkos::Timer timer;
    if (test != CUSPARSE) {
    timer.reset();
    sptrsv_symbolic( &kh, row_map, entries );
    std::cout << "UTRI Symbolic Time: " << timer.seconds() << std::endl;

    //std::cout << "TriSolve Solve" << std::endl;
    timer.reset();
    sptrsv_solve( &kh, row_map, entries, values, rhs, lhs );
    Kokkos::fence();
    std::cout << "UTRI Solve Time: " << timer.seconds() << std::endl;
  
    }
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
// step 4: perform analysis
    else {
#if 0
      double *dvalues = (double *)(values.data());
      int *drow_map = (int *)(row_map.data());
      int *dentries = (int *)(entries.data());
      double *dlhs = (double *)(lhs.data());
      double *drhs = (double *)(rhs.data());
#endif
      //int nnz = triMtx.nnz();
      //std::cout << "  cusparse path: analysis" << std::endl;
      //status = cusparseDcsrsv2_analysis(handle, trans, nrows, nnz, descr, (double*)dvalues, (int *)drow_map, (int *)dentries, info, policy, pBuffer);
      timer.reset();
      status = cusparseDcsrsv2_analysis(handle, trans, nrows, triMtx.nnz(), descr, values.data(), row_map.data(), entries.data(), info, policy, pBuffer);
      std::cout << "UTRI Cusparse Symbolic Time: " << timer.seconds() << std::endl;
      if (CUSPARSE_STATUS_SUCCESS != status)
        std::cout << "analysis status error name " << (status) << std::endl;
      // L has unit diagonal, so no structural zero is reported.

      status = cusparseXcsrsv2_zeroPivot(handle, info, &structural_zero);
      if (CUSPARSE_STATUS_ZERO_PIVOT == status){
         printf("L(%d,%d) is missing\n", structural_zero, structural_zero);
      }

// step 5: solve L*y = x
      //std::cout << "  cusparse path: solve" << std::endl;
      //status = cusparseDcsrsv2_solve(handle, trans, nrows, nnz, &alpha, descr, (double*)dvalues, (int *)drow_map, (int *)dentries, info, (double*)drhs, (double*)dlhs, policy, pBuffer);
      timer.reset();
      status = cusparseDcsrsv2_solve(handle, trans, nrows, triMtx.nnz(), &alpha, descr, values.data(), row_map.data(), entries.data(), info, rhs.data(), lhs.data(), policy, pBuffer);
      Kokkos::fence();
      std::cout << "UTRI Cusparse Solve Time: " << timer.seconds() << std::endl;
      if (CUSPARSE_STATUS_SUCCESS != status)
        std::cout << "solve status error name " << (status) << std::endl;
      // L has unit diagonal, so no numerical zero is reported.
      status = cusparseXcsrsv2_zeroPivot(handle, info, &numerical_zero);
      if (CUSPARSE_STATUS_ZERO_PIVOT == status){
         printf("L(%d,%d) is zero\n", numerical_zero, numerical_zero);
      }
    }
#endif
    // Error Check
    Kokkos::fence();
    scalar_t sum = 0.0;
    Kokkos::parallel_reduce( Kokkos::RangePolicy<execution_space>(0, lhs.extent(0)), 
      KOKKOS_LAMBDA ( const lno_t i, scalar_t &tsum ) {
        tsum += lhs(i);
      }, sum);
  
    if ( sum != scalar_t(lhs.extent(0)) ) {
      std::cout << "Upper Tri Solve FAILURE: sum = " << sum << std::endl;
      auto hsoln = Kokkos::create_mirror_view(lhs);
      Kokkos::deep_copy(hsoln, lhs);
      for ( size_t i = 0; i < hsoln.extent(0); ++i ) {
        std::cout << "lhs(" << i << ") = " << hsoln(i) << std::endl;
      }
      return 1;
    }
    else {
     std::cout << "Upper Tri Solve SUCCESS!" << std::endl;
    }
  
    // Benchmark
    Kokkos::fence();
    double min_time = 1.0e32;
    double max_time = 0.0;
    double ave_time = 0.0;

    for(int i=0;i<loop;i++) {
      timer.reset();
  
    if (test != CUSPARSE) {
      sptrsv_solve( &kh, row_map, entries, values, rhs, lhs );
    }
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
    else {
      cusparseDcsrsv2_solve(handle, trans, nrows, triMtx.nnz(), &alpha, descr, values.data(), row_map.data(), entries.data(), info, rhs.data(), lhs.data(), policy, pBuffer);
    }
#endif
  
      Kokkos::fence();
      double time = timer.seconds();
      ave_time += time;
      if(time>max_time) max_time = time;
      if(time<min_time) min_time = time;
    }

    std::cout << "LOOP_AVG_TIME:  " << ave_time/loop << std::endl;
    std::cout << "LOOP_MAX_TIME:  " << max_time << std::endl;
    std::cout << "LOOP_MIN_TIME:  " << min_time << std::endl;
  }
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
// step 6: free resources
    cudaFree(pBuffer);
    cusparseDestroyCsrsv2Info(info);
    cusparseDestroyMatDescr(descr);
    cusparseDestroy(handle);
#endif
  }

  return 0;
}


void print_help_sptrsv() {
  printf("Options:\n");
  printf("  --test [OPTION] : Use different kernel implementations\n");
  printf("                    Options:\n");
  printf("                      lvlrp, lvltp1, lvltp2\n\n");
  printf("                      cusparse           (Vendor Libraries)\n\n");
  printf("  -lf [file]       : Read in Matrix Market formatted text file 'file'.\n");
  printf("  -uf [file]       : Read in Matrix Market formatted text file 'file'.\n");
//  printf("  -s [N]          : generate a semi-random banded (band size 0.01xN) NxN matrix\n");
//  printf("                    with average of 10 entries per row.\n");
//  printf("  --schedule [SCH]: Set schedule for kk variant (static,dynamic,auto [ default ]).\n");
//  printf("  -fb [file]      : Read in binary Matrix files 'file'.\n");
//  printf("  --write-binary  : In combination with -f, generate binary files.\n");
  printf("  --offset [O]    : Subtract O from every index.\n");
  printf("                    Useful in case the matrix market file is not 0 based.\n\n");
  printf("  -rpt [K]        : Number of Rows assigned to a thread.\n");
  printf("  -ts [T]         : Number of threads per team.\n");
  printf("  -vl [V]         : Vector-length (i.e. how many Cuda threads are a Kokkos 'thread').\n");
  printf("  --loop [LOOP]       : How many spmv to run to aggregate average time. \n");
}


int main(int argc, char **argv)
{
 std::vector<int> tests;

 std::string lfilename;
 std::string ufilename;

 int vector_length = -1;
 int team_size = -1;
 int idx_offset = 0;
 int loop = 1;
// int schedule=AUTO;

 if(argc == 1) {
   print_help_sptrsv();
   return 0;
 }

 for(int i=0;i<argc;i++)
 {
  if((strcmp(argv[i],"--test")==0)) {
    i++;
    if((strcmp(argv[i],"lvlrp")==0)) {
      tests.push_back( LVLSCHED_RP );
    }
    if((strcmp(argv[i],"lvltp1")==0)) {
      tests.push_back( LVLSCHED_TP1 );
    }
/*
    if((strcmp(argv[i],"lvltp2")==0)) {
      tests.push_back( LVLSCHED_TP2 );
    }
*/
    if((strcmp(argv[i],"cusparse")==0)) {
      tests.push_back( CUSPARSE );
    }
    continue;
  }
  if((strcmp(argv[i],"-lf")==0)) {lfilename = argv[++i]; continue;}
  if((strcmp(argv[i],"-uf")==0)) {ufilename = argv[++i]; continue;}
  if((strcmp(argv[i],"-ts")==0)) {team_size=atoi(argv[++i]); continue;}
  if((strcmp(argv[i],"-vl")==0)) {vector_length=atoi(argv[++i]); continue;}
  if((strcmp(argv[i],"--offset")==0)) {idx_offset=atoi(argv[++i]); continue;}
  if((strcmp(argv[i],"--loop")==0)) {loop=atoi(argv[++i]); continue;}
/*
  if((strcmp(argv[i],"-lfb")==0)) {lfilename = argv[++i]; binaryfile = true; continue;}
  if((strcmp(argv[i],"-ufb")==0)) {ufilename = argv[++i]; binaryfile = true; continue;}
  if((strcmp(argv[i],"--schedule")==0)) {
    i++;
    if((strcmp(argv[i],"auto")==0))
      schedule = AUTO;
    if((strcmp(argv[i],"dynamic")==0))
      schedule = DYNAMIC;
    if((strcmp(argv[i],"static")==0))
      schedule = STATIC;
    continue;
  }
*/
  if((strcmp(argv[i],"--help")==0) || (strcmp(argv[i],"-h")==0)) {
    print_help_sptrsv();
    return 0;
  }
 }

 if (tests.size() == 0) {
   tests.push_back(DEFAULT);
 }
 for (size_t i = 0; i < tests.size(); ++i) {
    std::cout << "tests[" << i << "] = " << tests[i] << std::endl;
 }


 Kokkos::initialize(argc,argv);
 {
   int total_errors = test_sptrsv_perf<double>(tests,lfilename,ufilename,team_size,vector_length,idx_offset,loop);

   if(total_errors == 0)
   printf("Kokkos::SPTRSV Test: Passed\n");
   else
   printf("Kokkos::SPTRSV Test: Failed\n");


  }
  Kokkos::finalize();
  return 0;
}
#else
int main() {
  return 0;
}
#endif
