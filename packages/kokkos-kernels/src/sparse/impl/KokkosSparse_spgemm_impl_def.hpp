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


namespace KokkosSparse{


namespace Impl{


template <typename HandleType,
typename a_row_view_t_, typename a_lno_nnz_view_t_, typename a_scalar_nnz_view_t_,
typename b_lno_row_view_t_, typename b_lno_nnz_view_t_, typename b_scalar_nnz_view_t_  >
template <typename c_row_view_t, typename c_lno_nnz_view_t, typename c_scalar_nnz_view_t>
void KokkosSPGEMM
  <HandleType, a_row_view_t_, a_lno_nnz_view_t_, a_scalar_nnz_view_t_,
    b_lno_row_view_t_, b_lno_nnz_view_t_, b_scalar_nnz_view_t_>::
    KokkosSPGEMM_numeric(
      c_row_view_t &rowmapC_,
      c_lno_nnz_view_t &entriesC_,
      c_scalar_nnz_view_t &valuesC_){

    //get the algorithm and execution space.
    //SPGEMMAlgorithm spgemm_algorithm = this->handle->get_spgemm_handle()->get_algorithm_type();
    KokkosKernels::Impl::ExecSpaceType my_exec_space = KokkosKernels::Impl::get_exec_space_type<MyExecSpace>();

    if (KOKKOSKERNELS_VERBOSE){
      std::cout << "Numeric PHASE" << std::endl;
    }

    if (spgemm_algorithm == SPGEMM_KK_SPEED)
    {
      this->KokkosSPGEMM_numeric_speed(rowmapC_, entriesC_, valuesC_, my_exec_space);
    }
    else if ( spgemm_algorithm == SPGEMM_KK_COLOR ||
              spgemm_algorithm == SPGEMM_KK_MULTICOLOR ||
              spgemm_algorithm == SPGEMM_KK_MULTICOLOR2){
      this->KokkosSPGEMM_numeric_color(rowmapC_, entriesC_, valuesC_, spgemm_algorithm);
    }
    else if (spgemm_algorithm == SPGEMM_KK_MEMORY2){
        this->KokkosSPGEMM_numeric_hash2(rowmapC_, entriesC_, valuesC_, my_exec_space);
    }
    else if (spgemm_algorithm == SPGEMM_KK_OUTERMULTIMEM ){
      this->KokkosSPGEMM_numeric_outer(rowmapC_, entriesC_, valuesC_, my_exec_space);
    }
    else {
      this->KokkosSPGEMM_numeric_hash(rowmapC_, entriesC_, valuesC_, my_exec_space);
    }

  }

template <typename HandleType,
typename a_row_view_t_, typename a_lno_nnz_view_t_, typename a_scalar_nnz_view_t_,
typename b_lno_row_view_t_, typename b_lno_nnz_view_t_, typename b_scalar_nnz_view_t_  >
template <typename c_row_view_t>
void KokkosSPGEMM
  <HandleType, a_row_view_t_, a_lno_nnz_view_t_, a_scalar_nnz_view_t_,
    b_lno_row_view_t_, b_lno_nnz_view_t_, b_scalar_nnz_view_t_>::
    KokkosSPGEMM_symbolic(c_row_view_t rowmapC_){

  //SPGEMMAlgorithm spgemm_algorithm = this->handle->get_spgemm_handle()->get_algorithm_type();
  {

    //number of rows and nnzs
    nnz_lno_t n = this->row_mapB.dimension_0() - 1;
    size_type nnz = this->entriesB.dimension_0();
    KokkosKernels::Impl::ExecSpaceType my_exec_space = KokkosKernels::Impl::get_exec_space_type<MyExecSpace>();
    bool compress_in_single_step = this->handle->get_spgemm_handle()->get_compression_step();
    if (my_exec_space == KokkosKernels::Impl::Exec_CUDA) {
	compress_in_single_step = true;
    }
    //compressed b
    row_lno_temp_work_view_t new_row_mapB(Kokkos::ViewAllocateWithoutInitializing("new row map"), n+1);
    row_lno_temp_work_view_t new_row_mapB_begins;

    nnz_lno_temp_work_view_t set_index_entries; //will be output of compress matrix.
    nnz_lno_temp_work_view_t set_entries; //will be output of compress matrix


    if (KOKKOSKERNELS_VERBOSE){
      std::cout << "SYMBOLIC PHASE" << std::endl;
    }
    //First Compress B.
    Kokkos::Impl::Timer timer1;

    if (KOKKOSKERNELS_VERBOSE){
      std::cout << "\tCOMPRESS MATRIX-B PHASE" << std::endl;
    }
    //get the compressed matrix.
    this->compressMatrix(n, nnz, this->row_mapB, this->entriesB, new_row_mapB, set_index_entries, set_entries, compress_in_single_step);

    if (KOKKOSKERNELS_VERBOSE){
      std::cout << "\t\tCOMPRESS MATRIX-B overall time:" << timer1.seconds()
                                  << std::endl << std::endl;
    }

    timer1.reset();

    //first get the max flops for a row, which will be used for max row size.
    //If we did compression in single step, row_mapB[i] points the begining of row i,
    //and new_row_mapB[i] points to the end of row i.
    nnz_lno_t maxNumRoughZeros = 0;
    if (compress_in_single_step){
      maxNumRoughZeros = this->getMaxRoughRowNNZ(a_row_cnt, row_mapA, entriesA, row_mapB, new_row_mapB);
      if (KOKKOSKERNELS_VERBOSE){
        std::cout << "\tMax Row Flops:" << maxNumRoughZeros  << std::endl;
        std::cout << "\tMax Row Flop Calc Time:" << timer1.seconds()  << std::endl;
      }

      //calling symbolic structure
      this->symbolic_c(a_row_cnt, row_mapA, entriesA,
          row_mapB, new_row_mapB, set_index_entries, set_entries,
          rowmapC_, maxNumRoughZeros);

    }
    else {
      nnz_lno_t begin = 0;
      auto new_row_mapB_begin = Kokkos::subview (new_row_mapB, std::make_pair (begin, n - 1));
      auto new_row_mapB_end = Kokkos::subview (new_row_mapB, std::make_pair (begin + 1, n));
      //KokkosKernels::Impl::print_1Dview(new_row_mapB);
      //KokkosKernels::Impl::print_1Dview(new_row_mapB_begin);
      //KokkosKernels::Impl::print_1Dview(new_row_mapB_end);
      //But for 2 step it is a bit different.
      //new_row_mapB is complete and holds content of row i is in between new_row_mapB[i] - new_row_mapB[i+1]
      maxNumRoughZeros = this->getMaxRoughRowNNZ(a_row_cnt, row_mapA, entriesA, new_row_mapB_begin, new_row_mapB_end);
      if (KOKKOSKERNELS_VERBOSE){
        std::cout << "\tMax Row Flops:" << maxNumRoughZeros  << std::endl;
        std::cout << "\tMax Row Flop Calc Time:" << timer1.seconds()  << std::endl;
        std::cout << "\t Compression Ratio: " << set_index_entries.dimension_0() << " / " << nnz
            << " = " << set_index_entries.dimension_0() / double (nnz) << std::endl;
      }

      //calling symbolic structure
      this->symbolic_c(a_row_cnt, row_mapA, entriesA,
          new_row_mapB_begin, new_row_mapB_end, set_index_entries, set_entries,
          rowmapC_, maxNumRoughZeros);
    }

#ifdef KOKKOSKERNELS_ANALYZE_MEMORYACCESS
    double read_write_cost = this->handle->get_spgemm_handle()->get_read_write_cost_calc();
    if (read_write_cost){
      this->print_read_write_cost(rowmapC_);
    }
#endif
  }
}


template <typename HandleType,
typename a_row_view_t_, typename a_lno_nnz_view_t_, typename a_scalar_nnz_view_t_,
typename b_lno_row_view_t_, typename b_lno_nnz_view_t_, typename b_scalar_nnz_view_t_  >
template <typename c_row_view_t, typename c_nnz_view_t>
void KokkosSPGEMM
  <HandleType, a_row_view_t_, a_lno_nnz_view_t_, a_scalar_nnz_view_t_,
    b_lno_row_view_t_, b_lno_nnz_view_t_, b_scalar_nnz_view_t_>::
    write_matrix_to_plot(
    nnz_lno_t &num_colors,
    nnz_lno_persistent_work_host_view_t &h_color_xadj,
    nnz_lno_persistent_work_view_t &color_adj,
    c_row_view_t &rowmapC, c_nnz_view_t &entryIndicesC_){
  std::cout << "writing to plot" << std::endl;

  nnz_lno_persistent_work_host_view_t h_color_adj = Kokkos::create_mirror_view (color_adj);
  Kokkos::deep_copy (h_color_adj, color_adj);
  auto h_rowmapC = Kokkos::create_mirror_view (rowmapC);
  Kokkos::deep_copy (h_rowmapC, rowmapC);
  auto h_entryIndicesC = Kokkos::create_mirror_view (entryIndicesC_);
  Kokkos::deep_copy (h_entryIndicesC, entryIndicesC_);

  for (nnz_lno_t i = 0; i < num_colors; ++i){
    nnz_lno_t color_begin = h_color_xadj(i);
    nnz_lno_t color_end = h_color_xadj(i + 1);

    std::string colorind = "";
    std::stringstream ss;
    ss << i;


    ss >> colorind;
    colorind += ".coords";
    std::fstream fs;
    fs.open(colorind.c_str(), std::fstream::out);

    std::cout << "COLOR:" << i << " colorbegin:" << color_begin << " colorend:" << color_end << " size:" << color_end - color_begin << std::endl;
    for (nnz_lno_t j = color_begin; j < color_end; ++j){
      nnz_lno_t row = h_color_adj(j);
      for (size_type k = h_rowmapC(row); k < h_rowmapC(row + 1); ++k){
        nnz_lno_t column = h_entryIndicesC(k);
        //std::cout << row << " " << column << std::endl;
        fs << row << " " << column << std::endl;
      }
    }
    fs.close();
  }



  std::fstream fs;
  fs.open("plot1.gnuplot", std::fstream::out);
  for (nnz_lno_t i = 0; i < num_colors; ++i){
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


}
}

