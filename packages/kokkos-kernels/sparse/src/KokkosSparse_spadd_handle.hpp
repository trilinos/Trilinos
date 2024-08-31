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

#ifndef _SPADDHANDLE_HPP
#define _SPADDHANDLE_HPP

namespace KokkosSparse {

template <class lno_row_view_t_, class lno_nnz_view_t_, class scalar_nnz_view_t_, class ExecutionSpace,
          class MemorySpace>
class SPADDHandle {
 public:
  typedef typename lno_nnz_view_t_::non_const_type nnz_lno_view_t;
  typedef typename lno_row_view_t_::non_const_type nnz_row_view_t;
  typedef typename lno_row_view_t_::non_const_value_type size_type;
  typedef ExecutionSpace execution_space;

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
  struct SpaddCusparseData {
    size_t nbytes;
    void* workspace;
    cusparseMatDescr_t descrA, descrB, descrC;

    SpaddCusparseData() : nbytes(0), workspace(nullptr), descrA(nullptr), descrB(nullptr), descrC(nullptr) {}

    ~SpaddCusparseData() {
      Kokkos::kokkos_free<MemorySpace>(workspace);
      cusparseDestroyMatDescr(descrA);
      cusparseDestroyMatDescr(descrB);
      cusparseDestroyMatDescr(descrC);
    }
  };
#endif

#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE
  struct SpaddRocsparseData {
    rocsparse_mat_descr descrA, descrB, descrC;

    SpaddRocsparseData() : descrA(nullptr), descrB(nullptr), descrC(nullptr) {}

    ~SpaddRocsparseData() {
      rocsparse_destroy_mat_descr(descrA);
      rocsparse_destroy_mat_descr(descrB);
      rocsparse_destroy_mat_descr(descrC);
    }
  };
#endif

 private:
  // if both are true, the input matrices are strict CRS
  bool input_sorted;  // column indices in a row are sorted
  bool input_merged;  // column indices in a row are unique (i.e., merged)

  size_type result_nnz_size;

  bool called_symbolic;
  bool called_numeric;

  // a_pos and b_pos are used by the unsorted version of the kernel
  // both have same length as a_entries and b_entries
  // each entry provides the index in C row where the corresponding entry is
  // added
  nnz_lno_view_t a_pos;
  nnz_lno_view_t b_pos;

 public:
  /// \brief sets the result nnz size.
  /// \param a_pos_in The offset into a.
  /// \param b_pos_in The offset into b.
  void set_a_b_pos(const nnz_lno_view_t& a_pos_in, const nnz_lno_view_t& b_pos_in) {
    a_pos = a_pos_in;
    b_pos = b_pos_in;
  }

  nnz_lno_view_t get_a_pos() { return a_pos; }

  nnz_lno_view_t get_b_pos() { return b_pos; }

  /// \brief sets the result nnz size.
  /// \param result_nnz_size_ size of the output matrix.
  void set_c_nnz(size_type result_nnz_size_) { this->result_nnz_size = result_nnz_size_; }

  /**
   * \brief returns the result nnz size.
   */
  size_type get_c_nnz() { return this->result_nnz_size; }

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
  SpaddCusparseData cusparseData;
#endif

#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE
  SpaddRocsparseData rocsparseData;
#endif

  /**
   * \brief Default constructor.
   */
  SPADDHandle(bool input_is_sorted, bool input_is_merged = false)
      : input_sorted(input_is_sorted),
        input_merged(input_is_merged),
        result_nnz_size(0),
        called_symbolic(false),
        called_numeric(false) {}

  virtual ~SPADDHandle(){};

  bool is_symbolic_called() { return this->called_symbolic; }
  bool is_numeric_called() { return this->called_numeric; }

  // setters
  void set_call_symbolic(bool call = true) { this->called_symbolic = call; }
  void set_call_numeric(bool call = true) { this->called_numeric = call; }

  bool is_input_sorted() { return input_sorted; }
  bool is_input_merged() { return input_merged; }
  bool is_input_strict_crs() { return input_sorted && input_merged; }
};

}  // namespace KokkosSparse

#endif
