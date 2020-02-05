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

#include <Kokkos_Core.hpp>
#include <KokkosSparse_spmv.hpp>
#include <KokkosKernels_Test_Structured_Matrix.hpp>
#include <KokkosSparse_spmv_struct_impl.hpp>
#include <KokkosSparse_spmv_impl.hpp>

enum {STRUCT, UNSTR};
enum {AUTO, DYNAMIC, STATIC};

#if defined(KOKKOSKERNELS_INST_ORDINAL_INT)
  typedef int default_lno_t;
#elif defined(KOKKOSKERNELS_INST_ORDINAL_INT64_T)
  typedef int64_t default_lno_t;
#else
  #error "Expect int and/or int64_t to be enabled as ORDINAL (lno_t) types"
#endif
  //Prefer int as the default offset type, because cuSPARSE doesn't support size_t for rowptrs.
#if defined(KOKKOSKERNELS_INST_OFFSET_INT)
  typedef int default_size_type;
#elif defined(KOKKOSKERNELS_INST_OFFSET_SIZE_T)
  typedef size_t default_size_type;
#else
  #error "Expect size_t and/or int to be enabled as OFFSET (size_type) types"
#endif

void print_help() {
  printf("SPMV_struct benchmark code written by Luc Berger-Vergiat.\n");
  printf("Options:\n");
  printf("  --check-errors  : Determine if the result of spmv_struct is compared to serial unstructured spmv.\n");
  printf("  --compare       : Compare results efficiency of spmv_struct and spmv.\n");
  printf("  -l [LOOP]       : How many spmv to run to aggregate average time. \n");
  printf("  -nx             : How many nodes in x direction. \n");
  printf("  -ny             : How many nodes in y direction. \n");
  printf("  -nz             : How many nodes in z direction. \n");
  printf("  -st             : The stencil type used for discretization: 1 -> FD, 2 -> FE.\n");
  printf("  -dim            : Number of spacial dimensions used in the problem: 1, 2 or 3\n");
  printf("  -ws             : Worksets. \n");
  printf("  -ts             : Team Size. \n");
  printf("  -vl             : Vector length. \n");
}

int main(int argc, char **argv)
{
  typedef double Scalar;

  int  nx = 100;
  int  ny = 100;
  int  nz = 100;
  int  stencil_type = 1;
  int  numDimensions = 2;
  int  numVecs = 1;
  bool check_errors = false;
  bool compare = false;
  int  loop = 100;
  int  vl = -1;
  int  ts = -1;
  int  ws = -1;

  if(argc == 1) {
    print_help();
    return 0;
  }

  for(int i=0;i<argc;i++)
    {
      if((strcmp(argv[i],"-nx" )==0)) {nx=atoi(argv[++i]); continue;}
      if((strcmp(argv[i],"-ny" )==0)) {ny=atoi(argv[++i]); continue;}
      if((strcmp(argv[i],"-nz" )==0)) {nz=atoi(argv[++i]); continue;}
      if((strcmp(argv[i],"-st" )==0)) {stencil_type=atoi(argv[++i]); continue;}
      if((strcmp(argv[i],"-dim")==0)) {numDimensions=atoi(argv[++i]); continue;}
      if((strcmp(argv[i],"-mv" )==0)) {numVecs=atoi(argv[++i]); continue;}
      if((strcmp(argv[i],"-l"  )==0)) {loop=atoi(argv[++i]); continue;}
      if((strcmp(argv[i],"-vl" )==0)) {vl=atoi(argv[++i]); continue;}
      if((strcmp(argv[i],"-ts" )==0)) {ts=atoi(argv[++i]); continue;}
      if((strcmp(argv[i],"-ws" )==0)) {ws=atoi(argv[++i]); continue;}
      if((strcmp(argv[i],"--check-errors")==0)) {check_errors=true; continue;}
      if((strcmp(argv[i],"--compare")==0)) {compare=true; continue;}
      if((strcmp(argv[i],"--help")==0) || (strcmp(argv[i],"-h")==0)) {
        print_help();
        return 0;
      }
    }

  if(vl < 0) {
    vl = 2;
  } else if(ts < 0) {
    ts = 256 / vl;
    if(ws < 0) {
      if(numDimensions == 1) {
	ws = (nx - 2 + ts - 1) / ts;
      } else if(numDimensions == 2) {
	ws = ((nx - 2)*(ny - 2) + ts - 1) / ts;
      } else if(numDimensions == 3) {
	ws = ((nx - 2)*(ny - 2)*(nz - 2) + ts - 1) / ts;
      }
    }
  }

  Kokkos::initialize(argc,argv);
  {
    typedef default_lno_t lno_t;
    typedef default_size_type size_type;
    typedef KokkosSparse::CrsMatrix<Scalar,lno_t,Kokkos::DefaultExecutionSpace,void,size_type> matrix_type;
    typedef typename Kokkos::View<Scalar**,Kokkos::LayoutLeft> mv_type;
    // typedef typename Kokkos::View<Scalar*,Kokkos::LayoutLeft,Kokkos::MemoryRandomAccess > mv_random_read_type;
    typedef typename mv_type::HostMirror h_mv_type;

    int leftBC = 1, rightBC = 1, frontBC = 1, backBC = 1, bottomBC = 1, topBC = 1;

    Kokkos::View<lno_t*, Kokkos::HostSpace> structure("Spmv Structure", numDimensions);
    Kokkos::View<lno_t*[3], Kokkos::HostSpace> mat_structure("Matrix Structure", numDimensions);
    if(numDimensions == 1) {
      structure(0) = nx;
      mat_structure(0, 0) = nx;
    } else if(numDimensions == 2) {
      structure(0) = nx;
      structure(1) = ny;
      mat_structure(0, 0) = nx;
      mat_structure(1, 0) = ny;
      if(leftBC   == 1) { mat_structure(0, 1) = 1; }
      if(rightBC  == 1) { mat_structure(0, 2) = 1; }
      if(bottomBC == 1) { mat_structure(1, 1) = 1; }
      if(topBC    == 1) { mat_structure(1, 2) = 1; }
    } else if(numDimensions == 3) {
      structure(0) = nx;
      structure(1) = ny;
      structure(2) = nz;
      mat_structure(0, 0) = nx;
      mat_structure(1, 0) = ny;
      mat_structure(2, 0) = nz;
      if(leftBC   == 1) { mat_structure(0, 1) = 1; }
      if(rightBC  == 1) { mat_structure(0, 2) = 1; }
      if(frontBC  == 1) { mat_structure(1, 1) = 1; }
      if(backBC   == 1) { mat_structure(1, 2) = 1; }
      if(bottomBC == 1) { mat_structure(2, 1) = 1; }
      if(topBC    == 1) { mat_structure(2, 2) = 1; }
    }

    std::string discrectization_stencil;
    if(stencil_type == 1) {
      discrectization_stencil = "FD";
    } else if(stencil_type == 2) {
      discrectization_stencil = "FE";
    }

    matrix_type A;
    if(numDimensions == 1) {
      A = Test::generate_structured_matrix1D<matrix_type>(mat_structure);
    } else if(numDimensions == 2) {
      A = Test::generate_structured_matrix2D<matrix_type>(discrectization_stencil,
							  mat_structure);
    } else if(numDimensions == 3) {
      A = Test::generate_structured_matrix3D<matrix_type>(discrectization_stencil,
							  mat_structure);
    }

    mv_type x("X", A.numCols(), numVecs);
    // mv_random_read_type t_x(x);
    mv_type y("Y", A.numRows(), numVecs);
    h_mv_type h_x = Kokkos::create_mirror_view(x);
    h_mv_type h_y = Kokkos::create_mirror_view(y);
    h_mv_type h_y_compare;

    for(int rowIdx = 0; rowIdx < A.numCols(); ++rowIdx) {
      for(int vecIdx = 0; vecIdx < numVecs; ++vecIdx) {
        h_x(rowIdx, vecIdx) = static_cast<Scalar>(1.0*(rand()%40)-20.0);
        h_y(rowIdx, vecIdx) = static_cast<Scalar>(1.0*(rand()%40)-20.0);
      }
    }

    if(check_errors) {
      h_y_compare = Kokkos::create_mirror(y);
      typename matrix_type::StaticCrsGraphType::HostMirror h_graph = Kokkos::create_mirror(A.graph);
      typename matrix_type::values_type::HostMirror h_values = Kokkos::create_mirror_view(A.values);

      // Error Check Gold Values
      for(int rowIdx = 0; rowIdx < A.numRows(); ++rowIdx) {
        int start = h_graph.row_map(rowIdx);
        int end = h_graph.row_map(rowIdx + 1);

        for(int vecIdx = 0; vecIdx < numVecs; ++vecIdx) {
          h_y_compare(rowIdx, vecIdx) = 0;
        }
        for(int entryIdx = start; entryIdx < end; ++entryIdx) {
          // Scalar tmp_val = h_graph.entries(entryIdx) + i;
          int idx = h_graph.entries(entryIdx);
          for(int vecIdx = 0; vecIdx < numVecs; ++vecIdx) {
            h_y_compare(rowIdx, vecIdx) += h_values(entryIdx)*h_x(idx, vecIdx);
          }
        }
      }
    }

    Kokkos::deep_copy(x, h_x);
    Kokkos::deep_copy(y, h_y);
    Kokkos::View<Scalar**, Kokkos::LayoutLeft, Kokkos::DefaultExecutionSpace> x1("X1", A.numCols(), numVecs);
    Kokkos::View<Scalar**, Kokkos::LayoutLeft, Kokkos::DefaultExecutionSpace> y1("Y1", A.numRows(), numVecs);
    Kokkos::deep_copy(x1, h_x);
    // typename KokkosSparse::CrsMatrix<Scalar,int,Kokkos::DefaultExecutionSpace,void,int>::values_type y1("Y1", A.numRows(), numVecs);

    {
      Kokkos::Profiling::pushRegion("Structured spmv test");
      // Benchmark
      double min_time = 1.0e32;
      double max_time = 0.0;
      double ave_time = 0.0;
      for(int i=0;i<loop;i++) {
	Kokkos::Timer timer;
	KokkosSparse::Experimental::spmv_struct("N", stencil_type, structure, 1.0, A, x1, 1.0, y1);
	Kokkos::fence();
	double time = timer.seconds();
	ave_time += time;
	if(time>max_time) max_time = time;
	if(time<min_time) min_time = time;
      }

      // Performance Output
      double matrix_size = 1.0*((A.nnz()*(sizeof(Scalar) + sizeof(int)) + A.numRows()*sizeof(int)))/1024/1024;
      double vector_size = 2.0*A.numRows()*sizeof(Scalar)/1024/1024;
      double vector_readwrite = (A.nnz() + A.numCols())*sizeof(Scalar)/1024/1024;

      double problem_size = matrix_size+vector_size;
      printf("Type NNZ NumRows NumCols ProblemSize(MB) AveBandwidth(GB/s) MinBandwidth(GB/s) MaxBandwidth(GB/s) AveGFlop MinGFlop MaxGFlop aveTime(ms) maxTime(ms) minTime(ms)\n");
      printf("Struct %zu %zu %zu %6.2lf ( %6.2lf %6.2lf %6.2lf ) ( %6.3lf %6.3lf %6.3lf ) ( %6.3lf %6.3lf %6.3lf )\n",
	     (size_t) A.nnz(), (size_t) A.numRows(), (size_t) A.numCols(), problem_size,
	     (matrix_size+vector_readwrite)/ave_time*loop/1024, (matrix_size+vector_readwrite)/max_time/1024, (matrix_size+vector_readwrite)/min_time/1024,
	     2.0*A.nnz()*loop/ave_time/1e9, 2.0*A.nnz()/max_time/1e9, 2.0*A.nnz()/min_time/1e9,
	     ave_time/loop*1000, max_time*1000, min_time*1000);
      Kokkos::Profiling::popRegion();
    }

    if(compare) {
      Kokkos::Profiling::pushRegion("Unstructured spmv test");
      // Benchmark
      double min_time = 1.0e32;
      double max_time = 0.0;
      double ave_time = 0.0;
      for(int i=0;i<loop;i++) {
	Kokkos::Timer timer;
	KokkosSparse::spmv("N", 1.0, A, x1, 1.0, y1);
	Kokkos::fence();
	double time = timer.seconds();
	ave_time += time;
	if(time>max_time) max_time = time;
	if(time<min_time) min_time = time;
      }

      // Performance Output
      double matrix_size = 1.0*((A.nnz()*(sizeof(Scalar)+sizeof(int)) + A.numRows()*sizeof(int)))/1024/1024;
      double vector_size = 2.0*A.numRows()*sizeof(Scalar)/1024/1024;
      double vector_readwrite = (A.nnz() + A.numCols())*sizeof(Scalar)/1024/1024;

      double problem_size = matrix_size+vector_size;
      printf("Unstr %zu %zu %zu %6.2lf ( %6.2lf %6.2lf %6.2lf ) ( %6.3lf %6.3lf %6.3lf ) ( %6.3lf %6.3lf %6.3lf )\n",
	     (size_t) A.nnz(), (size_t) A.numRows(), (size_t) A.numCols(), problem_size,
	     (matrix_size+vector_readwrite)/ave_time*loop/1024, (matrix_size+vector_readwrite)/max_time/1024,(matrix_size+vector_readwrite)/min_time/1024,
	     2.0*A.nnz()*loop/ave_time/1e9, 2.0*A.nnz()/max_time/1e9, 2.0*A.nnz()/min_time/1e9,
	     ave_time/loop*1000, max_time*1000, min_time*1000);
      Kokkos::Profiling::popRegion();
    }

    if(check_errors) {
      // Error Check
      Kokkos::deep_copy(h_y, y1);
      Scalar error = 0;
      Scalar sum = 0;
      for(int rowIdx = 0; rowIdx < A.numRows(); ++rowIdx) {
        for(int vecIdx = 0; vecIdx < numVecs; ++vecIdx) {
          error += (h_y_compare(rowIdx, vecIdx) - h_y(rowIdx, vecIdx))*(h_y_compare(rowIdx, vecIdx) - h_y(rowIdx, vecIdx));
          sum += h_y_compare(rowIdx, vecIdx)*h_y_compare(rowIdx, vecIdx);
        }
      }

      int num_errors = 0;
      double total_error = 0;
      double total_sum = 0;
      num_errors += (error/(sum==0?1:sum))>1e-5?1:0;
      total_error += error;
      total_sum += sum;

      if(total_error == 0) {
	printf("Kokkos::MultiVector Test: Passed\n");
      } else {
	printf("Kokkos::MultiVector Test: Failed\n");
      }
    }

  }
  Kokkos::finalize();
}
