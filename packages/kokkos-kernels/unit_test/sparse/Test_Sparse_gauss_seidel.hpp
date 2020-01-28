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

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>
#include "KokkosKernels_Handle.hpp"
#include "KokkosKernels_IOUtils.hpp"
//#include <Kokkos_Sparse_CrsMatrix.hpp>
#include <KokkosSparse_spmv.hpp>
#include <KokkosBlas1_dot.hpp>
#include <KokkosBlas1_axpby.hpp>
#include <KokkosBlas1_nrm2.hpp>
#include <cstdlib>
#include <iostream>
#include <complex>
#include <map>
#include <vector>
#include "KokkosSparse_gauss_seidel.hpp"
#include "KokkosSparse_partitioning_impl.hpp"
#include "impl/KokkosSparse_sor_sequential_impl.hpp"

#ifndef kokkos_complex_double
#define kokkos_complex_double Kokkos::complex<double>
#define kokkos_complex_float Kokkos::complex<float>
#endif

using namespace KokkosKernels;
using namespace KokkosKernels::Experimental;
using namespace KokkosSparse;
using namespace KokkosSparse::Experimental;
namespace Test {

template <typename crsMat_t, typename vec_t, typename device>
int run_gauss_seidel(
    crsMat_t input_mat,
    GSAlgorithm gs_algorithm,
    vec_t x_vector,
    vec_t y_vector,
    bool is_symmetric_graph,
    int apply_type = 0, // 0 for symmetric, 1 for forward, 2 for backward.
    int cluster_size = 1,
    ClusteringAlgorithm cluster_algorithm = CLUSTER_DEFAULT)
{
  typedef typename crsMat_t::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type lno_view_t;
  typedef typename graph_t::entries_type lno_nnz_view_t;
  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;

  typedef typename lno_view_t::value_type size_type;
  typedef typename lno_nnz_view_t::value_type lno_t;
  typedef typename scalar_view_t::value_type scalar_t;

  typedef KokkosKernelsHandle
      <size_type,lno_t, scalar_t,
      typename device::execution_space, typename device::memory_space,typename device::memory_space > KernelHandle;

  KernelHandle kh;
  kh.set_team_work_size(16);
  kh.set_dynamic_scheduling(true);
  if(gs_algorithm == GS_CLUSTER)
    kh.create_gs_handle(cluster_algorithm, cluster_size);
  else
    kh.create_gs_handle(GS_DEFAULT);

  const size_t num_rows_1 = input_mat.numRows();
  const size_t num_cols_1 = input_mat.numCols();
  const int apply_count = 100;

  gauss_seidel_symbolic
    (&kh, num_rows_1, num_cols_1, input_mat.graph.row_map, input_mat.graph.entries, is_symmetric_graph);
  gauss_seidel_numeric
    (&kh, num_rows_1, num_cols_1, input_mat.graph.row_map, input_mat.graph.entries, input_mat.values, is_symmetric_graph);

  scalar_t omega(0.9);

  switch (apply_type){
  case 0:
    symmetric_gauss_seidel_apply
      (&kh, num_rows_1, num_cols_1, input_mat.graph.row_map, input_mat.graph.entries, input_mat.values, x_vector, y_vector, false, true, omega, apply_count);
    break;
  case 1:
    forward_sweep_gauss_seidel_apply
    (&kh, num_rows_1, num_cols_1, input_mat.graph.row_map, input_mat.graph.entries, input_mat.values, x_vector, y_vector, false, true, omega, apply_count);
    break;
  case 2:
    backward_sweep_gauss_seidel_apply
    (&kh, num_rows_1, num_cols_1, input_mat.graph.row_map, input_mat.graph.entries, input_mat.values, x_vector, y_vector, false, true, omega, apply_count);
    break;
  default:
    symmetric_gauss_seidel_apply
    (&kh, num_rows_1, num_cols_1, input_mat.graph.row_map, input_mat.graph.entries, input_mat.values, x_vector, y_vector, false, true, omega, apply_count);
    break;
  }
  kh.destroy_gs_handle();
  return 0;
}

template<typename vec_t>
vec_t create_x_vector(vec_t& kok_x, double max_value = 10.0) {
  typedef typename vec_t::value_type scalar_t;
  auto h_x = Kokkos::create_mirror_view (kok_x);
  for (size_t j = 0; j < h_x.extent(1); ++j){
    for (size_t i = 0; i < h_x.extent(0); ++i){
      scalar_t r =
          static_cast <scalar_t> (rand()) /
          static_cast <scalar_t> (RAND_MAX / max_value);
      h_x.access(i, j) = r;
    }
  }
  Kokkos::deep_copy (kok_x, h_x);
  return kok_x;
}

template <typename crsMat_t, typename vector_t>
vector_t create_y_vector(crsMat_t crsMat, vector_t x_vector){
  vector_t y_vector (Kokkos::ViewAllocateWithoutInitializing("Y VECTOR"),
      crsMat.numRows());
  KokkosSparse::spmv("N", 1, crsMat, x_vector, 0, y_vector);
  return y_vector;
}

template <typename crsMat_t, typename vector_t>
vector_t create_y_vector_mv(crsMat_t crsMat, vector_t x_vector){
  vector_t y_vector (Kokkos::ViewAllocateWithoutInitializing("Y VECTOR"),
      crsMat.numRows(), x_vector.extent(1));
  KokkosSparse::spmv("N", 1, crsMat, x_vector, 0, y_vector);
  return y_vector;
}
}

template<typename scalar_t, typename lno_t, typename size_type, typename device, typename crsMat_t>
crsMat_t symmetrize(crsMat_t A)
{
  typedef typename crsMat_t::StaticCrsGraphType graph_t;
  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;
  typedef typename graph_t::row_map_type::non_const_type lno_view_t;
  typedef typename graph_t::entries_type::non_const_type lno_nnz_view_t;
  auto host_rowmap = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A.graph.row_map);
  auto host_entries = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A.graph.entries);
  auto host_values = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A.values);
  lno_t numRows = A.numRows();
  //symmetrize as input_mat + input_mat^T, to still have a diagonally dominant matrix
  typedef std::map<lno_t, scalar_t> Row;
  std::vector<Row> symRows(numRows);
  for(lno_t r = 0; r < numRows; r++)
  {
    auto& row = symRows[r];
    for(size_type i = host_rowmap(r); i < host_rowmap(r + 1); i++)
    {
      lno_t c = host_entries(i);
      auto& col = symRows[c];
      auto it = row.find(c);
      if(it == row.end())
        row[c] = host_values(i);
      else
        row[c] += host_values(i);
      it = col.find(r);
      if(it == col.end())
        col[r] = host_values(i);
      else
        col[r] += host_values(i);
    }
  }
  //Count entries
  Kokkos::View<size_type*, Kokkos::LayoutLeft, Kokkos::HostSpace> new_host_rowmap("Rowmap", numRows + 1);
  size_t accum = 0;
  for(lno_t r = 0; r <= numRows; r++)
  {
    new_host_rowmap(r) = accum;
    if(r < numRows)
      accum += symRows[r].size();
  }
  //Allocate new entries/values
  Kokkos::View<lno_t*, Kokkos::LayoutLeft, Kokkos::HostSpace> new_host_entries("Entries", accum);
  Kokkos::View<scalar_t*, Kokkos::LayoutLeft, Kokkos::HostSpace> new_host_values("Values", accum);
  for(lno_t r = 0; r < numRows; r++)
  {
    auto rowIt = symRows[r].begin();
    for(size_type i = new_host_rowmap(r); i < new_host_rowmap(r + 1); i++)
    {
      new_host_entries(i) = rowIt->first;
      new_host_values(i) = rowIt->second;
      rowIt++;
    }
  }
  lno_view_t new_rowmap("Rowmap", numRows + 1);
  lno_nnz_view_t new_entries("Entries", accum);
  scalar_view_t new_values("Values", accum);
  Kokkos::deep_copy(new_rowmap, new_host_rowmap);
  Kokkos::deep_copy(new_entries, new_host_entries);
  Kokkos::deep_copy(new_values, new_host_values);
  return crsMat_t("SymA", numRows, numRows, accum, new_values, new_rowmap, new_entries);
}

template <typename scalar_t, typename lno_t, typename size_type, typename device>
void test_gauss_seidel_rank1(lno_t numRows, size_type nnz, lno_t bandwidth, lno_t row_size_variance, bool symmetric)
{
  using namespace Test;
  typedef typename KokkosSparse::CrsMatrix<scalar_t, lno_t, device, void, size_type> crsMat_t;
  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;
  typedef typename Kokkos::Details::ArithTraits<scalar_t>::mag_type mag_t;
  srand(245);
  lno_t numCols = numRows;
  crsMat_t input_mat = KokkosKernels::Impl::kk_generate_diagonally_dominant_sparse_matrix<crsMat_t>(numRows, numCols, nnz, row_size_variance, bandwidth);
  if(symmetric)
  {
    //Symmetrize on host, rather than relying on the parallel versions (those can be tested for symmetric=false)
    input_mat = symmetrize<scalar_t, lno_t, size_type, device, crsMat_t>(input_mat);
  }
  lno_t nv = input_mat.numRows();
  scalar_view_t solution_x(Kokkos::ViewAllocateWithoutInitializing("X (correct)"), nv);
  create_x_vector(solution_x);
  mag_t initial_norm_res = KokkosBlas::nrm2(solution_x);
  scalar_view_t y_vector = create_y_vector(input_mat, solution_x);
  //GS_DEFAULT is GS_TEAM on CUDA and GS_PERMUTED on other spaces, and the behavior
  //of each algorithm _should be_ the same on every execution space, which is why
  //we just test GS_DEFAULT.
  int apply_count = 3;  //test symmetric, forward, backward
  scalar_view_t x_vector(Kokkos::ViewAllocateWithoutInitializing("x vector"), nv);
  const scalar_t alpha = 1.0;
  //*** Point-coloring version ****
  for (int apply_type = 0; apply_type < apply_count; ++apply_type)
  {
    Kokkos::Impl::Timer timer1;
    Kokkos::deep_copy(x_vector, scalar_t(0.0));
    run_gauss_seidel<crsMat_t, scalar_view_t, device>(input_mat, GS_DEFAULT, x_vector, y_vector, symmetric, apply_type);
    //double gs = timer1.seconds();
    //KokkosKernels::Impl::print_1Dview(x_vector);
    KokkosBlas::axpby(alpha, solution_x, -alpha, x_vector);
    mag_t result_norm_res = KokkosBlas::nrm2(x_vector);
    EXPECT_LT(result_norm_res, initial_norm_res);
  }
  //*** Cluster-coloring version ****
  int clusterSizes[3] = {2, 5, 34};
  for(int csize = 0; csize < 3; csize++)
  {
    for(int algo = 0; algo < (int) NUM_CLUSTERING_ALGORITHMS; algo++)
    {
      for(int apply_type = 0; apply_type < apply_count; ++apply_type)
      {
        Kokkos::Impl::Timer timer1;
        //Zero out X before solving
        Kokkos::deep_copy(x_vector, scalar_t(0.0));
        run_gauss_seidel<crsMat_t, scalar_view_t, device>(
            input_mat, GS_CLUSTER, x_vector, y_vector, symmetric, apply_type, clusterSizes[csize], (ClusteringAlgorithm) algo);
        KokkosBlas::axpby(alpha, solution_x, -alpha, x_vector);
        mag_t result_norm_res = KokkosBlas::nrm2(x_vector);
        EXPECT_LT(result_norm_res, initial_norm_res);
      }
    }
  }
}

template <typename scalar_t, typename lno_t, typename size_type, typename device>
void test_gauss_seidel_rank2(lno_t numRows, size_type nnz, lno_t bandwidth, lno_t row_size_variance, lno_t numVecs, bool symmetric)
{
  using namespace Test;
  srand(245);
  typedef typename KokkosSparse::CrsMatrix<scalar_t, lno_t, device, void, size_type> crsMat_t;
  typedef Kokkos::View<scalar_t**, Kokkos::LayoutLeft, device> scalar_view2d_t;
  typedef Kokkos::View<scalar_t**, Kokkos::LayoutLeft, Kokkos::HostSpace> host_scalar_view2d_t;
  typedef typename Kokkos::Details::ArithTraits<scalar_t>::mag_type mag_t;

  lno_t numCols = numRows;
  crsMat_t input_mat = KokkosKernels::Impl::kk_generate_diagonally_dominant_sparse_matrix<crsMat_t>(numRows,numCols,nnz,row_size_variance, bandwidth);
  if(symmetric)
  {
    //Symmetrize on host, rather than relying on the parallel versions (those can be tested for symmetric=false)
    input_mat = symmetrize<scalar_t, lno_t, size_type, device, crsMat_t>(input_mat);
  }
  lno_t nv = input_mat.numRows();
  host_scalar_view2d_t solution_x(Kokkos::ViewAllocateWithoutInitializing("X (correct)"), nv, numVecs);
  create_x_vector(solution_x);
  scalar_view2d_t x_vector(Kokkos::ViewAllocateWithoutInitializing("X"), nv, numVecs);
  Kokkos::deep_copy(x_vector, solution_x);
  scalar_view2d_t y_vector = create_y_vector_mv(input_mat, x_vector);
  auto x_host = Kokkos::create_mirror_view(x_vector);
  std::vector<mag_t> initial_norms(numVecs);
  for(lno_t i = 0; i < numVecs; i++)
  {
    scalar_t sum = 0;
    for(lno_t j = 0; j < nv; j++)
    {
      sum += solution_x(j, i) * solution_x(j, i);
    }
    initial_norms[i] = Kokkos::Details::ArithTraits<mag_t>::sqrt(
        Kokkos::Details::ArithTraits<scalar_t>::abs(sum));
  }
  int apply_count = 3;  //test symmetric, forward, backward
  //*** Point-coloring version ****
  for(int apply_type = 0; apply_type < apply_count; ++apply_type)
  {
    Kokkos::Impl::Timer timer1;
    //Zero out X before solving
    Kokkos::deep_copy(x_vector, scalar_t(0.0));
    run_gauss_seidel<crsMat_t, scalar_view2d_t, device>(
        input_mat, GS_DEFAULT, x_vector, y_vector, symmetric, apply_type);
    Kokkos::deep_copy(x_host, x_vector);
    for(lno_t i = 0; i < numVecs; i++)
    {
      scalar_t diffDot = 0;
      for(lno_t j = 0; j < numRows; j++)
      {
        scalar_t diff = x_host(j, i) - solution_x(j, i);
        diffDot += diff * diff;
      }
      mag_t res = Kokkos::Details::ArithTraits<mag_t>::sqrt(
          Kokkos::Details::ArithTraits<scalar_t>::abs(diffDot));
      EXPECT_LT(res, initial_norms[i]);
    }
  }
  //*** Cluster-coloring version ****
  int clusterSizes[3] = {2, 5, 34};
  for(int csize = 0; csize < 3; csize++)
  {
    for(int algo = 0; algo < (int) NUM_CLUSTERING_ALGORITHMS; algo++)
    {
      for(int apply_type = 0; apply_type < apply_count; ++apply_type)
      {
        Kokkos::Impl::Timer timer1;
        //Zero out X before solving
        Kokkos::deep_copy(x_vector, scalar_t(0.0));
        run_gauss_seidel<crsMat_t, scalar_view2d_t, device>(
            input_mat, GS_CLUSTER, x_vector, y_vector, symmetric, apply_type, clusterSizes[csize], (ClusteringAlgorithm) algo);
        Kokkos::deep_copy(x_host, x_vector);
        for(lno_t i = 0; i < numVecs; i++)
        {
          scalar_t diffDot = 0;
          for(lno_t j = 0; j < numRows; j++)
          {
            scalar_t diff = x_host(j, i) - solution_x(j, i);
            diffDot += diff * diff;
          }
          mag_t res = Kokkos::Details::ArithTraits<mag_t>::sqrt(
              Kokkos::Details::ArithTraits<scalar_t>::abs(diffDot));
          EXPECT_LT(res, initial_norms[i]);
        }
      }
    }
  }
}

template <typename scalar_t, typename lno_t, typename size_type, typename device>
void test_rcm(lno_t numRows, size_type nnzPerRow, lno_t bandwidth)
{
  using namespace Test;
  typedef typename KokkosSparse::CrsMatrix<scalar_t, lno_t, device, void, size_type> crsMat_t;
  typedef typename crsMat_t::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type::non_const_type lno_row_view_t;
  typedef typename graph_t::entries_type::non_const_type lno_nnz_view_t;
  typedef KokkosKernelsHandle
      <size_type, lno_t, scalar_t,
      typename device::execution_space, typename device::memory_space,typename device::memory_space> KernelHandle;
  srand(245);
  size_type nnzTotal = nnzPerRow * numRows;
  lno_t nnzVariance = nnzPerRow / 4;
  crsMat_t A = KokkosKernels::Impl::kk_generate_sparse_matrix<crsMat_t>(numRows, numRows, nnzTotal, nnzVariance, bandwidth);
  lno_row_view_t symRowmap;
  lno_nnz_view_t symEntries;
  KokkosKernels::Impl::symmetrize_graph_symbolic_hashmap
    <typename graph_t::row_map_type, typename graph_t::entries_type, lno_row_view_t, lno_nnz_view_t, typename device::execution_space>
    (numRows, A.graph.row_map, A.graph.entries, symRowmap, symEntries);
  typedef KokkosSparse::Impl::RCM<KernelHandle, typename graph_t::row_map_type::non_const_type, typename graph_t::entries_type::non_const_type> rcm_t;
  rcm_t rcm(numRows, symRowmap, symEntries);
  lno_nnz_view_t rcmOrder = rcm.rcm();
  //perm(i) = the node with timestamp i
  //make sure that perm is in fact a permutation matrix (contains each row exactly once)
  Kokkos::View<lno_t*, Kokkos::HostSpace> rcmHost("RCM row ordering", numRows);
  Kokkos::deep_copy(rcmHost, rcmOrder);
  std::set<lno_t> rowSet;
  for(lno_t i = 0; i < numRows; i++)
    rowSet.insert(rcmHost(i));
  if((lno_t) rowSet.size() != numRows)
  {
    std::cerr << "Only got back " << rowSet.size() << " unique row IDs!\n";
    return;
  }
}

template <typename scalar_t, typename lno_t, typename size_type, typename device>
void test_sequential_sor(lno_t numRows, size_type nnz, lno_t bandwidth, lno_t row_size_variance) {
  const scalar_t zero = Kokkos::Details::ArithTraits<scalar_t>::zero();
  const scalar_t one = Kokkos::Details::ArithTraits<scalar_t>::one();
  srand(245);
  typedef typename KokkosSparse::CrsMatrix<scalar_t, lno_t, device, void, size_type> crsMat_t;
  lno_t numCols = numRows;
  crsMat_t input_mat = KokkosKernels::Impl::kk_generate_diagonally_dominant_sparse_matrix<crsMat_t>(numRows,numCols,nnz,row_size_variance, bandwidth);
  auto rowmap = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), input_mat.graph.row_map);
  auto entries = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), input_mat.graph.entries);
  auto values = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), input_mat.values);
  //create raw x (unkown), y (rhs) vectors
  using vector_t = typename crsMat_t::values_type::non_const_type;
  //Create random x
  vector_t x("X", numRows);
  auto x_host = Kokkos::create_mirror_view(x);
  for(lno_t i = 0; i < numRows; i++)
  {
    x_host(i) = one * scalar_t(10.0 * rand() / RAND_MAX);
  }
  Kokkos::deep_copy(x, x_host);
  //record the correct solution, to compare against at the end
  vector_t xgold("X gold", numRows);
  Kokkos::deep_copy(xgold, x);
  vector_t y = Test::create_y_vector(input_mat, x);
  auto y_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), y);
  //initial solution is zero
  Kokkos::deep_copy(x_host, zero);
  //get the inverse diagonal (only needed on host)
  Kokkos::View<scalar_t*, Kokkos::HostSpace> invDiag("diag^-1", numRows);
  for(lno_t i = 0; i < numRows; i++)
  {
    for(size_type j = rowmap(i); j < rowmap(i + 1); j++)
    {
      if(entries(j) == i)
        invDiag(i) = one / values(j);
    }
  }
  for(int i = 0; i < 1; i++)
  {
    KokkosSparse::Impl::Sequential::gaussSeidel
      <lno_t, size_type, scalar_t, scalar_t, scalar_t>
      (numRows, 1, rowmap.data(), entries.data(), values.data(),
       y_host.data(), numRows,
       x_host.data(), numRows,
       invDiag.data(),
       one, //omega
       "F");
    KokkosSparse::Impl::Sequential::gaussSeidel
      <lno_t, size_type, scalar_t, scalar_t, scalar_t>
      (numRows, 1, rowmap.data(), entries.data(), values.data(),
       y_host.data(), numRows,
       x_host.data(), numRows,
       invDiag.data(),
       one, //omega
       "B");
  }
  //Copy solution back
  Kokkos::deep_copy(x, x_host);
  //Check against gold solution
  scalar_t xSq = KokkosBlas::dot(x, x);
  scalar_t solnDot = KokkosBlas::dot(x, xgold);
  double scaledSolutionDot = Kokkos::Details::ArithTraits<scalar_t>::abs(solnDot / xSq);
  EXPECT_TRUE(0.99 < scaledSolutionDot);
}

template <typename scalar_t, typename lno_t, typename size_type, typename device>
void test_balloon_clustering(lno_t numRows, size_type nnzPerRow, lno_t bandwidth)
{
  using namespace Test;
  typedef typename KokkosSparse::CrsMatrix<scalar_t, lno_t, device, void, size_type> crsMat_t;
  typedef typename crsMat_t::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type const_lno_row_view_t;
  typedef typename graph_t::entries_type const_lno_nnz_view_t;
  typedef typename graph_t::row_map_type::non_const_type lno_row_view_t;
  typedef typename graph_t::entries_type::non_const_type lno_nnz_view_t;
  typedef KokkosKernelsHandle
      <size_type, lno_t, scalar_t,
      typename device::execution_space, typename device::memory_space,typename device::memory_space> KernelHandle;
  srand(245);
  size_type nnzTotal = nnzPerRow * numRows;
  lno_t nnzVariance = nnzPerRow / 4;
  crsMat_t A = KokkosKernels::Impl::kk_generate_sparse_matrix<crsMat_t>(numRows, numRows, nnzTotal, nnzVariance, bandwidth);
  lno_row_view_t symRowmap;
  lno_nnz_view_t symEntries;
  KokkosKernels::Impl::symmetrize_graph_symbolic_hashmap
    <const_lno_row_view_t, const_lno_nnz_view_t, lno_row_view_t, lno_nnz_view_t, typename device::execution_space>
    (numRows, A.graph.row_map, A.graph.entries, symRowmap, symEntries);
  KokkosSparse::Impl::BalloonClustering<KernelHandle, lno_row_view_t, lno_nnz_view_t> balloon(numRows, symRowmap, symEntries);
  for(int clusterSize = 1; clusterSize <= numRows / 16; clusterSize = std::ceil(clusterSize * 1.3))
  {
    auto vertClusters = balloon.run(clusterSize);
    //validate results: make sure cluster labels are in bounds, and that the number of clusters is correct
    auto vertClustersHost = Kokkos::create_mirror_view(vertClusters);
    Kokkos::deep_copy(vertClustersHost, vertClusters);
    lno_t numClusters = (numRows + clusterSize - 1) / clusterSize;
    //check the hard constraints of the clustering
    std::set<lno_t> uniqueClusterIDs;
    for(lno_t i = 0; i < numRows; i++)
    {
      EXPECT_TRUE(vertClustersHost(i) >= 0);
      EXPECT_TRUE(vertClustersHost(i) < numClusters);
      uniqueClusterIDs.insert(vertClustersHost(i));
    }
    EXPECT_TRUE(uniqueClusterIDs.size() == static_cast<size_t>(numClusters));
  }
}

template <typename scalar_t, typename lno_t, typename size_type, typename device>
void test_sgs_zero_rows()
{
  using namespace Test;
  typedef typename KokkosSparse::CrsMatrix<scalar_t, lno_t, device, void, size_type> crsMat_t;
  typedef typename crsMat_t::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type::non_const_type row_map_type;
  typedef typename graph_t::entries_type::non_const_type entries_type;
  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;
  typedef KokkosKernelsHandle
      <size_type, lno_t, scalar_t,
      typename device::execution_space, typename device::memory_space,typename device::memory_space> KernelHandle;
  //The rowmap of a zero-row matrix can be length 0 or 1, so Gauss-Seidel should work with both
  //(the setup and apply are essentially no-ops but they shouldn't crash or throw exceptions)
  //For this test, create size-0 and size-1 rowmaps separately, and make sure each work with both point and cluster
  for(int doingCluster = 0; doingCluster < 2; doingCluster++)
  {
    for(int rowmapLen = 0; rowmapLen < 2; rowmapLen++)
    {
      KernelHandle kh;
      if(doingCluster)
        kh.create_gs_handle(CLUSTER_DEFAULT, 10);
      else
        kh.create_gs_handle(GS_DEFAULT);
      //initialized to 0
      row_map_type rowmap("Rowmap", rowmapLen);
      entries_type entries("Entries", 0);
      scalar_view_t values("Values", 0);
      //also, make sure graph symmetrization doesn't crash on zero rows
      gauss_seidel_symbolic(&kh, 0, 0, rowmap, entries, false);
      gauss_seidel_numeric(&kh, 0, 0, rowmap, entries, values, false);
      scalar_view_t x("X", 0);
      scalar_view_t y("Y", 0);
      scalar_t omega(0.9);
      symmetric_gauss_seidel_apply
        (&kh, 0, 0, rowmap, entries, values, x, y, false, true, omega, 3);
      kh.destroy_gs_handle();
    }
  }
}

#define EXECUTE_TEST(SCALAR, ORDINAL, OFFSET, DEVICE) \
TEST_F( TestCategory, sparse ## _ ## gauss_seidel_asymmetric_rank1 ## _ ## SCALAR ## _ ## ORDINAL ## _ ## OFFSET ## _ ## DEVICE ) { \
  test_gauss_seidel_rank1<SCALAR,ORDINAL,OFFSET,DEVICE>(2000, 2000 * 20, 200, 10, false); \
} \
TEST_F( TestCategory, sparse ## _ ## gauss_seidel_asymmetric_rank2 ## _ ## SCALAR ## _ ## ORDINAL ## _ ## OFFSET ## _ ## DEVICE ) { \
  test_gauss_seidel_rank2<SCALAR,ORDINAL,OFFSET,DEVICE>(2000, 2000 * 20, 200, 10, 3, false); \
} \
TEST_F( TestCategory, sparse ## _ ## gauss_seidel_symmetric_rank1 ## _ ## SCALAR ## _ ## ORDINAL ## _ ## OFFSET ## _ ## DEVICE ) { \
  test_gauss_seidel_rank1<SCALAR,ORDINAL,OFFSET,DEVICE>(2000, 2000 * 20, 200, 10, true); \
} \
TEST_F( TestCategory, sparse ## _ ## gauss_seidel_symmetric_rank2 ## _ ## SCALAR ## _ ## ORDINAL ## _ ## OFFSET ## _ ## DEVICE ) { \
  test_gauss_seidel_rank2<SCALAR,ORDINAL,OFFSET,DEVICE>(2000, 2000 * 20, 200, 10, 3, true); \
} \
TEST_F( TestCategory, sparse ## _ ## gauss_seidel_zero_rows ## _ ## SCALAR ## _ ## ORDINAL ## _ ## OFFSET ## _ ## DEVICE ) { \
  test_sgs_zero_rows<SCALAR,ORDINAL,OFFSET,DEVICE>(); \
} \
TEST_F( TestCategory, sparse ## _ ## rcm ## _ ## SCALAR ## _ ## ORDINAL ## _ ## OFFSET ## _ ## DEVICE ) { \
  test_rcm<SCALAR,ORDINAL,OFFSET,DEVICE>(10000, 50, 2000); \
} \
TEST_F( TestCategory, sparse ## _ ## balloon_clustering ## _ ## SCALAR ## _ ## ORDINAL ## _ ## OFFSET ## _ ## DEVICE ) { \
  test_balloon_clustering<SCALAR,ORDINAL,OFFSET,DEVICE>(5000, 100, 2000); \
} \
TEST_F( TestCategory, sparse ## _ ## sequential_sor ## _ ## SCALAR ## _ ## ORDINAL ## _ ## OFFSET ## _ ## DEVICE ) { \
  test_sequential_sor<SCALAR,ORDINAL,OFFSET,DEVICE>(1000, 1000 * 15, 50, 10); \
}

#if (defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(double, int, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(double, int64_t, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(double, int, size_t, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(double, int64_t, size_t, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(float, int, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(float, int64_t, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(float, int, size_t, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(float, int64_t, size_t, TestExecSpace)
#endif


#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(kokkos_complex_double, int, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(kokkos_complex_double, int64_t, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(kokkos_complex_double, int, size_t, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(kokkos_complex_double, int64_t, size_t, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(kokkos_complex_float, int, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_INT) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(kokkos_complex_float, int64_t, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(kokkos_complex_float, int, size_t, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
 && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T) ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(kokkos_complex_float, int64_t, size_t, TestExecSpace)
#endif

