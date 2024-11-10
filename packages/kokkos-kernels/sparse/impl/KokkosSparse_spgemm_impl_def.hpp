//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

namespace KokkosSparse {

namespace Impl {

template <typename HandleType, typename a_row_view_t_, typename a_lno_nnz_view_t_, typename a_scalar_nnz_view_t_,
          typename b_lno_row_view_t_, typename b_lno_nnz_view_t_, typename b_scalar_nnz_view_t_>
template <typename c_row_view_t, typename c_lno_nnz_view_t, typename c_scalar_nnz_view_t>
void KokkosSPGEMM<HandleType, a_row_view_t_, a_lno_nnz_view_t_, a_scalar_nnz_view_t_, b_lno_row_view_t_,
                  b_lno_nnz_view_t_, b_scalar_nnz_view_t_>::KokkosSPGEMM_numeric(c_row_view_t &rowmapC_,
                                                                                 c_lno_nnz_view_t &entriesC_,
                                                                                 c_scalar_nnz_view_t &valuesC_) {
  // get the algorithm and execution space.
  // SPGEMMAlgorithm spgemm_algorithm =
  // this->handle->get_spgemm_handle()->get_algorithm_type();
  KokkosKernels::Impl::ExecSpaceType my_exec_space_ = KokkosKernels::Impl::get_exec_space_type<MyExecSpace>();

  if (KOKKOSKERNELS_VERBOSE) {
    std::cout << "Numeric PHASE" << std::endl;
  }

  if (spgemm_algorithm == SPGEMM_KK_SPEED || spgemm_algorithm == SPGEMM_KK_DENSE) {
    this->KokkosSPGEMM_numeric_speed(rowmapC_, entriesC_, valuesC_, my_exec_space_);
  } else {
    this->KokkosSPGEMM_numeric_hash(rowmapC_, entriesC_, valuesC_, my_exec_space_);
  }
}

template <typename HandleType, typename a_row_view_t_, typename a_lno_nnz_view_t_, typename a_scalar_nnz_view_t_,
          typename b_lno_row_view_t_, typename b_lno_nnz_view_t_, typename b_scalar_nnz_view_t_>
template <typename c_row_view_t>
void KokkosSPGEMM<HandleType, a_row_view_t_, a_lno_nnz_view_t_, a_scalar_nnz_view_t_, b_lno_row_view_t_,
                  b_lno_nnz_view_t_, b_scalar_nnz_view_t_>::KokkosSPGEMM_symbolic(c_row_view_t rowmapC_) {
  Kokkos::Profiling::pushRegion("KokkosSparse::spgemm_symbolic[NATIVE]");
  {
    if (KOKKOSKERNELS_VERBOSE) {
      std::cout << "SYMBOLIC PHASE" << std::endl;
    }
    // first calculate the number of original flops required.
    this->compute_row_flops();

    // number of rows and nnzs
    nnz_lno_t n   = this->row_mapB.extent(0) - 1;
    size_type nnz = this->entriesB.extent(0);

    bool compress_in_single_step = this->handle->get_spgemm_handle()->get_compression_step();
    // compress in single step if it is GPU.
    if (KokkosKernels::Impl::kk_is_gpu_exec_space<MyExecSpace>()) compress_in_single_step = true;

    // compressed B fields.
    row_lno_temp_work_view_t new_row_mapB(Kokkos::view_alloc(Kokkos::WithoutInitializing, "new row map"), n + 1);
    row_lno_temp_work_view_t new_row_mapB_begins;

    nnz_lno_temp_work_view_t set_index_entries;  // will be output of compress matrix.
    nnz_lno_temp_work_view_t set_entries;        // will be output of compress matrix

    // First Compress B.
    Kokkos::Timer timer1;

    if (KOKKOSKERNELS_VERBOSE) {
      std::cout << "\tCOMPRESS MATRIX-B PHASE" << std::endl;
    }

    // call compression.
    // it might not go through to the end if ratio is not high.
    bool compression_applied = this->compressMatrix(n, nnz, this->row_mapB, this->entriesB, new_row_mapB,
                                                    set_index_entries, set_entries, compress_in_single_step);

    if (KOKKOSKERNELS_VERBOSE) {
      std::cout << "\t\tCOMPRESS MATRIX-B overall time:" << timer1.seconds() << std::endl << std::endl;
    }

    timer1.reset();

    // first get the max flops for a row, which will be used for max row size.
    // If we did compression in single step, row_mapB[i] points the begining of
    // row i, and new_row_mapB[i] points to the end of row i.

    if (compression_applied) {
      nnz_lno_t maxNumRoughZeros = this->handle->get_spgemm_handle()->compressed_max_row_flops;

      if (compress_in_single_step) {
        // calling symbolic structure
        this->symbolic_c(a_row_cnt, row_mapA, entriesA, row_mapB, new_row_mapB, set_index_entries, set_entries,
                         rowmapC_, maxNumRoughZeros);

      } else {
        nnz_lno_t begin         = 0;
        auto new_row_mapB_begin = Kokkos::subview(new_row_mapB, std::make_pair(begin, n));
        auto new_row_mapB_end   = Kokkos::subview(new_row_mapB, std::make_pair(begin + 1, n + 1));

        // calling symbolic structure
        this->symbolic_c(a_row_cnt, row_mapA, entriesA, new_row_mapB_begin, new_row_mapB_end, set_index_entries,
                         set_entries, rowmapC_, maxNumRoughZeros);
      }
    } else {
      new_row_mapB               = row_lno_temp_work_view_t();
      new_row_mapB_begins        = row_lno_temp_work_view_t();
      set_index_entries          = nnz_lno_temp_work_view_t();
      set_entries                = nnz_lno_temp_work_view_t();
      nnz_lno_t maxNumRoughZeros = this->handle->get_spgemm_handle()->original_max_row_flops;
      if (KOKKOSKERNELS_VERBOSE) {
        std::cout << "SYMBOLIC PHASE -- NO COMPRESSION: maxNumRoughZeros:" << maxNumRoughZeros << std::endl;
      }

      auto new_row_mapB_begin = Kokkos::subview(this->row_mapB, std::make_pair(nnz_lno_t(0), n));
      auto new_row_mapB_end   = Kokkos::subview(this->row_mapB, std::make_pair(nnz_lno_t(1), n + 1));

      // calling symbolic structure
      this->symbolic_c_no_compression(a_row_cnt, row_mapA, entriesA, new_row_mapB_begin, new_row_mapB_end,
                                      this->entriesB, rowmapC_, maxNumRoughZeros);
    }
#ifdef KOKKOSKERNELS_ANALYZE_MEMORYACCESS
    double read_write_cost = this->handle->get_spgemm_handle()->get_read_write_cost_calc();
    if (read_write_cost) {
      this->print_read_write_cost(rowmapC_);
    }
#endif
  }
  Kokkos::Profiling::popRegion();
}

template <typename HandleType, typename a_row_view_t_, typename a_lno_nnz_view_t_, typename a_scalar_nnz_view_t_,
          typename b_lno_row_view_t_, typename b_lno_nnz_view_t_, typename b_scalar_nnz_view_t_>
template <typename c_row_view_t, typename c_nnz_view_t>
void KokkosSPGEMM<HandleType, a_row_view_t_, a_lno_nnz_view_t_, a_scalar_nnz_view_t_, b_lno_row_view_t_,
                  b_lno_nnz_view_t_,
                  b_scalar_nnz_view_t_>::write_matrix_to_plot(nnz_lno_t &num_colors,
                                                              nnz_lno_persistent_work_host_view_t &h_color_xadj,
                                                              nnz_lno_persistent_work_view_t &color_adj,
                                                              c_row_view_t &rowmapC, c_nnz_view_t &entryIndicesC_) {
  std::cout << "writing to plot" << std::endl;

  nnz_lno_persistent_work_host_view_t h_color_adj = Kokkos::create_mirror_view(color_adj);
  Kokkos::deep_copy(h_color_adj, color_adj);
  auto h_rowmapC = Kokkos::create_mirror_view(rowmapC);
  Kokkos::deep_copy(h_rowmapC, rowmapC);
  auto h_entryIndicesC = Kokkos::create_mirror_view(entryIndicesC_);
  Kokkos::deep_copy(h_entryIndicesC, entryIndicesC_);

  for (nnz_lno_t i = 0; i < num_colors; ++i) {
    nnz_lno_t color_begin = h_color_xadj(i);
    nnz_lno_t color_end   = h_color_xadj(i + 1);

    std::string colorind = "";
    std::stringstream ss;
    ss << i;

    ss >> colorind;
    colorind += ".coords";
    std::fstream fs;
    fs.open(colorind.c_str(), std::fstream::out);

    std::cout << "COLOR:" << i << " colorbegin:" << color_begin << " colorend:" << color_end
              << " size:" << color_end - color_begin << std::endl;
    for (nnz_lno_t j = color_begin; j < color_end; ++j) {
      nnz_lno_t row = h_color_adj(j);
      for (size_type k = h_rowmapC(row); k < h_rowmapC(row + 1); ++k) {
        nnz_lno_t column = h_entryIndicesC(k);
        // std::cout << row << " " << column << std::endl;
        fs << row << " " << column << std::endl;
      }
    }
    fs.close();
  }

  std::fstream fs;
  fs.open("plot1.gnuplot", std::fstream::out);
  for (nnz_lno_t i = 0; i < num_colors; ++i) {
    std::string colorind = "\"";
    std::stringstream ss;
    ss << i;

    ss >> colorind;
    colorind += ".coords\"";
    if (i > 0) fs << "re";
    fs << "plot " << colorind << std::endl;
  }
  fs << "pause -1" << std::endl;
  fs.close();
}

}  // namespace Impl
}  // namespace KokkosSparse
