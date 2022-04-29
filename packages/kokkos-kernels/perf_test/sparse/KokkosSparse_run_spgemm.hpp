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

#include "KokkosSparse_spgemm.hpp"
#include "KokkosKernels_TestParameters.hpp"
#include "KokkosKernels_Sorting.hpp"

#define TRANPOSEFIRST false
#define TRANPOSESECOND false

namespace KokkosKernels {

namespace Experiment {
template <typename crsMat_t, typename device>
bool is_same_matrix(crsMat_t output_mat1, crsMat_t output_mat2) {
  typedef typename crsMat_t::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type::non_const_type lno_view_t;
  typedef typename graph_t::entries_type::non_const_type lno_nnz_view_t;
  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;

  size_t nrows1    = output_mat1.graph.row_map.extent(0);
  size_t nentries1 = output_mat1.graph.entries.extent(0);
  size_t nvals1    = output_mat1.values.extent(0);

  size_t nrows2    = output_mat2.graph.row_map.extent(0);
  size_t nentries2 = output_mat2.graph.entries.extent(0);
  size_t nvals2    = output_mat2.values.extent(0);

  KokkosKernels::sort_crs_matrix(output_mat1);

  if (nrows1 != nrows2) {
    std::cerr << "row count is different" << std::endl;
    return false;
  }
  if (nentries1 != nentries2) {
    std::cerr << "nentries2 is different" << std::endl;
    return false;
  }
  if (nvals1 != nvals2) {
    std::cerr << "nvals1 is different" << std::endl;
    return false;
  }

  KokkosKernels::sort_crs_matrix(output_mat2);

  bool is_identical = true;
  is_identical      = KokkosKernels::Impl::kk_is_identical_view<
      typename graph_t::row_map_type, typename graph_t::row_map_type,
      typename lno_view_t::value_type, typename device::execution_space>(
      output_mat1.graph.row_map, output_mat2.graph.row_map, 0);
  if (!is_identical) {
    std::cerr << "rowmaps differ" << std::endl;
    return false;
  }

  is_identical = KokkosKernels::Impl::kk_is_identical_view<
      lno_nnz_view_t, lno_nnz_view_t, typename lno_nnz_view_t::value_type,
      typename device::execution_space>(output_mat1.graph.entries,
                                        output_mat2.graph.entries, 0);
  if (!is_identical) {
    for (size_t i = 0; i < nrows1; ++i) {
      size_t rb      = output_mat1.graph.row_map(i);
      size_t re      = output_mat1.graph.row_map(i + 1);
      bool incorrect = false;
      for (size_t j = rb; j < re; ++j) {
        if (output_mat1.graph.entries(j) != output_mat2.graph.entries(j)) {
          incorrect = true;
          break;
        }
      }
      if (incorrect) {
        for (size_t j = rb; j < re; ++j) {
          std::cerr << "row:" << i << " j:" << j
                    << " h_ent1(j):" << output_mat1.graph.entries(j)
                    << " h_ent2(j):" << output_mat2.graph.entries(j)
                    << " rb:" << rb << " re:" << re << std::endl;
        }
      }
    }
    std::cerr << "entries differ" << std::endl;
    return false;
  }

  is_identical = KokkosKernels::Impl::kk_is_identical_view<
      scalar_view_t, scalar_view_t, typename scalar_view_t::value_type,
      typename device::execution_space>(output_mat1.values, output_mat2.values,
                                        0.000001);
  if (!is_identical) {
    std::cerr << "Incorret values" << std::endl;
  }
  return true;
}

template <typename ExecSpace, typename crsMat_t, typename crsMat_t2,
          typename crsMat_t3, typename TempMemSpace,
          typename PersistentMemSpace>
crsMat_t3 run_experiment(crsMat_t crsMat, crsMat_t2 crsMat2,
                         Parameters params) {
  using namespace KokkosSparse;
  using namespace KokkosSparse::Experimental;
  using device_t = Kokkos::Device<ExecSpace, PersistentMemSpace>;
  int algorithm  = params.algorithm;
  int repeat     = params.repeat;
  int chunk_size = params.chunk_size;

  int shmemsize                 = params.shmemsize;
  int team_size                 = params.team_size;
  int use_dynamic_scheduling    = params.use_dynamic_scheduling;
  int verbose                   = params.verbose;
  int calculate_read_write_cost = params.calculate_read_write_cost;
  // char spgemm_step = params.spgemm_step;
  int vector_size     = params.vector_size;
  int check_output    = params.check_output;
  int mkl_keep_output = params.mkl_keep_output;
  // spgemm_step++;
  typedef typename crsMat_t3::values_type::non_const_type scalar_view_t;
  typedef typename crsMat_t3::row_map_type::non_const_type lno_view_t;
  typedef typename crsMat_t3::index_type::non_const_type lno_nnz_view_t;
  typedef typename lno_nnz_view_t::value_type lno_t;
  typedef typename lno_view_t::value_type size_type;
  typedef typename scalar_view_t::value_type scalar_t;

  lno_view_t row_mapC;
  lno_nnz_view_t entriesC;
  scalar_view_t valuesC;

  typedef KokkosKernels::Experimental::KokkosKernelsHandle<
      size_type, lno_t, scalar_t, ExecSpace, TempMemSpace, PersistentMemSpace>
      KernelHandle;

  typedef typename lno_nnz_view_t::value_type idx;
  typedef typename lno_view_t::value_type size_type;

  KernelHandle kh;
  kh.set_team_work_size(chunk_size);
  kh.set_shmem_size(shmemsize);
  kh.set_suggested_team_size(team_size);
  kh.set_suggested_vector_size(vector_size);

  if (use_dynamic_scheduling) {
    kh.set_dynamic_scheduling(true);
  }
  if (verbose) {
    kh.set_verbose(true);
  }

  const idx m = crsMat.numRows();
  const idx n = crsMat2.numRows();
  const idx k = crsMat2.numCols();

  if (verbose) std::cout << "m:" << m << " n:" << n << " k:" << k << std::endl;
  if (n < crsMat.numCols()) {
    std::cerr << "left.numCols():" << crsMat.numCols()
              << " right.numRows():" << crsMat2.numRows() << std::endl;
    exit(1);
  }

  // The reference product (for verifying correctness)
  // Don't allocate them if they won't be used, but they must be declared here.
  lno_view_t row_mapC_ref;
  lno_nnz_view_t entriesC_ref;
  scalar_view_t valuesC_ref;
  // Reference output has same type as actual output
  crsMat_t3 Ccrsmat_ref;

  if (check_output) {
    if (verbose) std::cout << "Running a reference algorithm" << std::endl;
    row_mapC_ref = lno_view_t("non_const_lnow_row", m + 1);
    KernelHandle sequential_kh;
    sequential_kh.set_team_work_size(chunk_size);
    sequential_kh.set_shmem_size(shmemsize);
    sequential_kh.set_suggested_team_size(team_size);
    sequential_kh.create_spgemm_handle(KokkosSparse::SPGEMM_SERIAL);

    if (use_dynamic_scheduling) {
      sequential_kh.set_dynamic_scheduling(true);
    }

    spgemm_symbolic(&sequential_kh, m, n, k, crsMat.graph.row_map,
                    crsMat.graph.entries, TRANPOSEFIRST, crsMat2.graph.row_map,
                    crsMat2.graph.entries, TRANPOSESECOND, row_mapC_ref);

    ExecSpace().fence();

    size_type c_nnz_size = sequential_kh.get_spgemm_handle()->get_c_nnz();
    entriesC_ref         = lno_nnz_view_t(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "entriesC"),
        c_nnz_size);
    valuesC_ref = scalar_view_t(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "valuesC"), c_nnz_size);

    spgemm_numeric(&sequential_kh, m, n, k, crsMat.graph.row_map,
                   crsMat.graph.entries, crsMat.values, TRANPOSEFIRST,

                   crsMat2.graph.row_map, crsMat2.graph.entries, crsMat2.values,
                   TRANPOSESECOND, row_mapC_ref, entriesC_ref, valuesC_ref);
    ExecSpace().fence();

    Ccrsmat_ref = crsMat_t3("CorrectC", m, k, valuesC_ref.extent(0),
                            valuesC_ref, row_mapC_ref, entriesC_ref);
  }

  for (int i = 0; i < repeat; ++i) {
    kh.create_spgemm_handle(KokkosSparse::SPGEMMAlgorithm(algorithm));

    kh.get_spgemm_handle()->mkl_keep_output = mkl_keep_output;
    kh.get_spgemm_handle()->set_mkl_sort_option(params.mkl_sort_option);

    // if mkl2 input needs to be converted to 1base.
    kh.get_spgemm_handle()->mkl_convert_to_1base = true;

    // 250000 default. if cache-mode is used on KNL can increase to 1M.
    kh.get_spgemm_handle()->MaxColDenseAcc = params.MaxColDenseAcc;

    if (i == 0) {
      kh.get_spgemm_handle()->set_read_write_cost_calc(
          calculate_read_write_cost);
    }
    // do the compression whether in 2 step, or 1 step.
    kh.get_spgemm_handle()->set_compression_steps(!params.compression2step);
    // whether to scale the hash more. default is 1, so no scale.
    kh.get_spgemm_handle()->set_min_hash_size_scale(params.minhashscale);
    // max occupancy in 1-level LP hashes. LL hashes can be 100%
    kh.get_spgemm_handle()->set_first_level_hash_cut_off(
        params.first_level_hash_cut_off);
    // min reduction on FLOPs to run compression
    kh.get_spgemm_handle()->set_compression_cut_off(params.compression_cut_off);

    row_mapC = lno_view_t("non_const_lnow_row", m + 1);
    entriesC = lno_nnz_view_t("entriesC (empty)", 0);
    valuesC  = scalar_view_t("valuesC (empty)", 0);

    Kokkos::Timer timer1;
    spgemm_symbolic(&kh, m, n, k, crsMat.graph.row_map, crsMat.graph.entries,
                    TRANPOSEFIRST, crsMat2.graph.row_map, crsMat2.graph.entries,
                    TRANPOSESECOND, row_mapC);

    ExecSpace().fence();
    double symbolic_time = timer1.seconds();

    Kokkos::Timer timer3;
    size_type c_nnz_size = kh.get_spgemm_handle()->get_c_nnz();
    if (verbose) std::cout << "C SIZE:" << c_nnz_size << std::endl;
    if (c_nnz_size) {
      entriesC = lno_nnz_view_t(
          Kokkos::view_alloc(Kokkos::WithoutInitializing, "entriesC"),
          c_nnz_size);
      valuesC = scalar_view_t(
          Kokkos::view_alloc(Kokkos::WithoutInitializing, "valuesC"),
          c_nnz_size);
    }
    spgemm_numeric(&kh, m, n, k, crsMat.graph.row_map, crsMat.graph.entries,
                   crsMat.values, TRANPOSEFIRST,

                   crsMat2.graph.row_map, crsMat2.graph.entries, crsMat2.values,
                   TRANPOSESECOND, row_mapC, entriesC, valuesC);
    ExecSpace().fence();
    double numeric_time = timer3.seconds();

    std::cout << "mm_time:" << symbolic_time + numeric_time
              << " symbolic_time:" << symbolic_time
              << " numeric_time:" << numeric_time << std::endl;
  }
  if (verbose) {
    std::cout << "row_mapC:" << row_mapC.extent(0) << std::endl;
    std::cout << "entriesC:" << entriesC.extent(0) << std::endl;
    std::cout << "valuesC:" << valuesC.extent(0) << std::endl;
    KokkosKernels::Impl::print_1Dview(valuesC);
    KokkosKernels::Impl::print_1Dview(entriesC);
    KokkosKernels::Impl::print_1Dview(row_mapC);
  }
  crsMat_t3 Ccrsmat_result("CrsMatrixC", m, k, valuesC.extent(0), valuesC,
                           row_mapC, entriesC);
  if (check_output) {
    bool is_identical =
        is_same_matrix<crsMat_t3, device_t>(Ccrsmat_result, Ccrsmat_ref);
    if (!is_identical) {
      std::cerr << "Result differs. If values are differing, might be floating "
                   "point order error."
                << std::endl;
      exit(1);
    }
  }
  return Ccrsmat_result;
}

}  // namespace Experiment
}  // namespace KokkosKernels
