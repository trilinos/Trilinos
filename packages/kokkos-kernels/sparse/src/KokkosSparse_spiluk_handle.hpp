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

#include <Kokkos_Core.hpp>
#include <iostream>
#include <string>
#include <KokkosKernels_HashmapAccumulator.hpp>

#ifndef _SPILUKHANDLE_HPP
#define _SPILUKHANDLE_HPP

// #define EXPAND_FACT 3

namespace KokkosSparse {
namespace Experimental {

// TP2 algorithm has issues with some offset-ordinal combo to be addressed
enum class SPILUKAlgorithm {
  SEQLVLSCHD_TP1 /*, SEQLVLSCHED_TP2*/
};

template <class size_type_, class lno_t_, class scalar_t_, class ExecutionSpace, class TemporaryMemorySpace,
          class PersistentMemorySpace>
class SPILUKHandle {
 public:
  using HandleExecSpace             = ExecutionSpace;
  using HandleTempMemorySpace       = TemporaryMemorySpace;
  using HandlePersistentMemorySpace = PersistentMemorySpace;

  using execution_space = ExecutionSpace;
  using memory_space    = HandlePersistentMemorySpace;

  using TeamPolicy  = Kokkos::TeamPolicy<execution_space>;
  using RangePolicy = Kokkos::RangePolicy<execution_space>;

  using size_type       = typename std::remove_const<size_type_>::type;
  using const_size_type = const size_type;

  using nnz_lno_t       = typename std::remove_const<lno_t_>::type;
  using const_nnz_lno_t = const nnz_lno_t;

  using nnz_scalar_t       = typename std::remove_const<scalar_t_>::type;
  using const_nnz_scalar_t = const nnz_scalar_t;

  using nnz_row_view_t = Kokkos::View<size_type *, HandlePersistentMemorySpace>;

  using nnz_lno_view_t = Kokkos::View<nnz_lno_t *, HandlePersistentMemorySpace>;

  using nnz_value_view_t = typename Kokkos::View<nnz_scalar_t *, HandlePersistentMemorySpace>;

  using nnz_row_view_host_t = typename Kokkos::View<size_type *, Kokkos::HostSpace>;

  using nnz_lno_view_host_t = typename Kokkos::View<nnz_lno_t *, Kokkos::HostSpace>;

  using signed_integral_t = typename std::make_signed<typename nnz_row_view_t::non_const_value_type>::type;
  using signed_nnz_lno_view_t =
      Kokkos::View<signed_integral_t *, typename nnz_row_view_t::array_layout, typename nnz_row_view_t::device_type,
                   typename nnz_row_view_t::memory_traits>;

  using work_view_t = Kokkos::View<nnz_lno_t **, Kokkos::LayoutRight, HandlePersistentMemorySpace>;

 private:
  nnz_row_view_t level_list;  // level IDs which the rows belong to
  nnz_lno_view_t level_idx;   // the list of rows in each level
  nnz_lno_view_t level_ptr;   // the starting index (into the view level_idx) of each level
  // Make hlevel_ptr a separate allocation, since it will be accessed on host
  // between kernel launches.
  nnz_lno_view_host_t hlevel_ptr;
  nnz_lno_view_host_t level_nchunks;        // number of chunks of rows at each level
  nnz_lno_view_host_t level_nrowsperchunk;  // maximum number of rows among chunks at each level
  work_view_t iw;                           // working view for mapping dense indices to sparse indices

  size_type nrows;
  size_type nlevels;
  size_type nnzL;
  size_type nnzU;
  size_type block_size;
  size_type level_maxrows;          // max. number of rows among levels
  size_type level_maxrowsperchunk;  // max.number of rows among chunks among levels

  bool symbolic_complete;

  SPILUKAlgorithm algm;

  int team_size;
  int vector_size;

 public:
  SPILUKHandle(SPILUKAlgorithm choice, const size_type nrows_, const size_type nnzL_, const size_type nnzU_,
               const size_type block_size_ = 0, bool symbolic_complete_ = false)
      : level_list(),
        level_idx(),
        level_ptr(),
        hlevel_ptr(),
        level_nchunks(),
        level_nrowsperchunk(),
        iw(),
        nrows(nrows_),
        nlevels(0),
        nnzL(nnzL_),
        nnzU(nnzU_),
        block_size(block_size_),
        level_maxrows(0),
        level_maxrowsperchunk(0),
        symbolic_complete(symbolic_complete_),
        algm(choice),
        team_size(-1),
        vector_size(-1) {}

  void reset_handle(const size_type nrows_, const size_type nnzL_, const size_type nnzU_,
                    const size_type block_size_ = Kokkos::ArithTraits<size_type>::max()) {
    set_nrows(nrows_);
    set_num_levels(0);
    set_nnzL(nnzL_);
    set_nnzU(nnzU_);
    // user likely does not want to reset block size to 0, so set default
    // to size_type::max
    if (block_size_ != Kokkos::ArithTraits<size_type>::max()) {
      set_block_size(block_size_);
    }
    set_level_maxrows(0);
    set_level_maxrowsperchunk(0);
    level_list          = nnz_row_view_t("level_list", nrows_);
    level_idx           = nnz_lno_view_t("level_idx", nrows_);
    level_ptr           = nnz_lno_view_t("level_ptr", nrows_ + 1);
    hlevel_ptr          = nnz_lno_view_host_t("hlevel_ptr", nrows_ + 1);
    level_nchunks       = nnz_lno_view_host_t();
    level_nrowsperchunk = nnz_lno_view_host_t();
    iw                  = work_view_t();
    reset_symbolic_complete();
  }

  virtual ~SPILUKHandle(){};

  void set_algorithm(SPILUKAlgorithm choice) { algm = choice; }

  SPILUKAlgorithm get_algorithm() { return algm; }

  KOKKOS_INLINE_FUNCTION
  nnz_row_view_t get_level_list() const { return level_list; }

  KOKKOS_INLINE_FUNCTION
  nnz_lno_view_t get_level_idx() const { return level_idx; }

  KOKKOS_INLINE_FUNCTION
  nnz_lno_view_t get_level_ptr() const { return level_ptr; }

  inline nnz_lno_view_host_t get_host_level_ptr() const { return hlevel_ptr; }

  KOKKOS_INLINE_FUNCTION
  nnz_lno_view_host_t get_level_nchunks() const { return level_nchunks; }

  void alloc_level_nchunks(const size_type nlevels_) { level_nchunks = nnz_lno_view_host_t("level_nchunks", nlevels_); }

  KOKKOS_INLINE_FUNCTION
  nnz_lno_view_host_t get_level_nrowsperchunk() const { return level_nrowsperchunk; }

  void alloc_level_nrowsperchunk(const size_type nlevels_) {
    level_nrowsperchunk = nnz_lno_view_host_t("level_nrowsperchunk", nlevels_);
  }

  KOKKOS_INLINE_FUNCTION
  work_view_t get_iw() const { return iw; }

  void alloc_iw(const size_type nrows_, const size_type ncols_) {
    iw = work_view_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "iw"), nrows_, ncols_);
    Kokkos::deep_copy(iw, nnz_lno_t(-1));
  }

  KOKKOS_INLINE_FUNCTION
  size_type get_nrows() const { return nrows; }

  KOKKOS_INLINE_FUNCTION
  void set_nrows(const size_type nrows_) { this->nrows = nrows_; }

  KOKKOS_INLINE_FUNCTION
  size_type get_nnzL() const { return nnzL; }

  KOKKOS_INLINE_FUNCTION
  void set_nnzL(const size_type nnzL_) { this->nnzL = nnzL_; }

  KOKKOS_INLINE_FUNCTION
  size_type get_nnzU() const { return nnzU; }

  KOKKOS_INLINE_FUNCTION
  void set_nnzU(const size_type nnzU_) { this->nnzU = nnzU_; }

  KOKKOS_INLINE_FUNCTION
  size_type get_block_size() const { return block_size; }

  KOKKOS_INLINE_FUNCTION
  void set_block_size(const size_type block_size_) { this->block_size = block_size_; }

  KOKKOS_INLINE_FUNCTION
  size_type get_level_maxrows() const { return level_maxrows; }

  KOKKOS_INLINE_FUNCTION
  void set_level_maxrows(const size_type level_maxrows_) { this->level_maxrows = level_maxrows_; }

  KOKKOS_INLINE_FUNCTION
  size_type get_level_maxrowsperchunk() const { return level_maxrowsperchunk; }

  KOKKOS_INLINE_FUNCTION
  void set_level_maxrowsperchunk(const size_type level_maxrowsperchunk_) {
    this->level_maxrowsperchunk = level_maxrowsperchunk_;
  }

  bool is_symbolic_complete() const { return symbolic_complete; }

  bool is_block_enabled() const { return block_size > 0; }

  size_type get_num_levels() const { return nlevels; }
  void set_num_levels(size_type nlevels_) { this->nlevels = nlevels_; }

  void set_symbolic_complete() { this->symbolic_complete = true; }
  void reset_symbolic_complete() { this->symbolic_complete = false; }

  void set_team_size(const int ts) { this->team_size = ts; }
  int get_team_size() const { return this->team_size; }

  void set_vector_size(const int vs) { this->vector_size = vs; }
  int get_vector_size() const { return this->vector_size; }

  void print_algorithm() {
    if (algm == SPILUKAlgorithm::SEQLVLSCHD_TP1) std::cout << "SEQLVLSCHD_TP1" << std::endl;

    /*if ( algm == SPILUKAlgorithm::SEQLVLSCHED_TP2 ) {
      std::cout << "SEQLVLSCHED_TP2" << std::endl;;
      std::cout << "WARNING: With CUDA this is currently only reliable with
    int-int ordinal-offset pair" << std::endl;
    }
    */
  }
};

}  // namespace Experimental
}  // namespace KokkosSparse

#endif
