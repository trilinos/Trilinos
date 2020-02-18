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

#include "KokkosKernels_Handle.hpp"
#include "KokkosKernels_IOUtils.hpp"
#include "KokkosSparse_spadd.hpp"
#include "KokkosKernels_TestParameters.hpp"

namespace KokkosKernels {
namespace Experiment {

template <typename crsMat_t>
void run_experiment(Parameters params)
{
  using namespace KokkosSparse;
  using namespace KokkosSparse::Experimental;

  using size_type = typename crsMat_t::size_type;
  using lno_t = typename crsMat_t::ordinal_type;
  using scalar_t = typename crsMat_t::value_type;
  using device_t = typename crsMat_t::device_type;
  using exec_space = typename device_t::execution_space;
  using mem_space = typename device_t::memory_space;

  using KernelHandle = KokkosKernels::Experimental::KokkosKernelsHandle
      <size_type, lno_t, scalar_t, exec_space, mem_space, mem_space>;

  std::cout << "************************************* \n";
  std::cout << "************************************* \n";
  std::cout << "Loading A from " << params.a_mtx_bin_file << '\n';
  crsMat_t A = Impl::read_kokkos_crst_matrix<crsMat_t>(params.a_mtx_bin_file);
  std::cout << "Loading B from " << params.b_mtx_bin_file << '\n';
  crsMat_t B = Impl::read_kokkos_crst_matrix<crsMat_t>(params.b_mtx_bin_file);
  //Make sure dimensions are compatible
  if(A.numRows() != B.numRows())
  {
    std::cout << "ERROR: A and B have different numbers of rows\n";
    exit(1);
  }
  if(A.numCols() != B.numCols())
  {
    std::cout << "ERROR: A and B have different numbers of columns\n";
    exit(1);
  }
  lno_t m = A.numRows();
  lno_t n = A.numCols();
  std::cout << "Read in A and B: " << m << "x" << n << '\n';

  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;
  typedef typename crsMat_t::StaticCrsGraphType::row_map_type::non_const_type lno_view_t;
  typedef typename crsMat_t::StaticCrsGraphType::entries_type::non_const_type lno_nnz_view_t;
  typedef typename crsMat_t::StaticCrsGraphType::row_map_type const_lno_view_t;
  typedef typename crsMat_t::StaticCrsGraphType::entries_type const_lno_nnz_view_t;

  lno_view_t row_mapC;
  lno_nnz_view_t entriesC;
  scalar_view_t valuesC;

  KernelHandle kh;

  if(params.assume_sorted)
    std::cout << "Assuming input matrices are sorted.\n";
  else
    std::cout << "Assuming input matrices are not sorted.\n";
  kh.create_spadd_handle(params.assume_sorted);
  auto addHandle = kh.get_spadd_handle();

  row_mapC = lno_view_t("non_const_lnow_row", m + 1);

  Kokkos::Impl::Timer timer1;

  spadd_symbolic<KernelHandle, const_lno_view_t, const_lno_nnz_view_t, const_lno_view_t, const_lno_nnz_view_t, lno_view_t, lno_nnz_view_t>
    (&kh,
      A.graph.row_map, A.graph.entries,
      B.graph.row_map, B.graph.entries,
      row_mapC);

  exec_space().fence();
  double symbolic_time = timer1.seconds();

  size_type c_nnz = addHandle->get_max_result_nnz();
  std::cout << "Result matrix will have " << c_nnz << " entries.\n";

  entriesC = lno_nnz_view_t("entriesC (empty)", c_nnz);
  valuesC = scalar_view_t("valuesC (empty)", c_nnz);

  Kokkos::Impl::Timer timer3;

  spadd_numeric(&kh,
      A.graph.row_map, A.graph.entries, A.values, 1.0, //A, alpha
      B.graph.row_map, B.graph.entries, B.values, 1.0, //B, beta
      row_mapC, entriesC, valuesC);  //C

  exec_space().fence();
  double numeric_time = timer3.seconds();

  std::cout
    << "total_time:" << symbolic_time + numeric_time
    << " symbolic_time:" << symbolic_time
    << " numeric_time:" << numeric_time << std::endl;

  if (params.verbose)
  {
    std::cout << "row_mapC:" << row_mapC.extent(0) << std::endl;
    std::cout << "entriesC:" << entriesC.extent(0) << std::endl;
    std::cout << "valuesC:" << valuesC.extent(0) << std::endl;
    KokkosKernels::Impl::print_1Dview(valuesC);
    KokkosKernels::Impl::print_1Dview(entriesC);
    KokkosKernels::Impl::print_1Dview(row_mapC);
  }
  if(params.c_mtx_bin_file)
  {
    std::cout << "Writing C (" << m << "x" << n << ") to " << params.c_mtx_bin_file << "\n";
    crsMat_t C("C", m, n, c_nnz, valuesC, row_mapC, entriesC);
    Impl::write_kokkos_crst_matrix<crsMat_t>(C, params.c_mtx_bin_file);
  }
}

}}  // namespace KokkosKernels::Experiment
