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

#ifndef _KOKKOSSPGEMMMKL_HPP
#define _KOKKOSSPGEMMMKL_HPP

#include "KokkosKernels_config.h"
#include "KokkosSparse_Utils_mkl.hpp"

#ifdef KOKKOSKERNELS_ENABLE_TPL_MKL
#include "mkl_spblas.h"

namespace KokkosSparse {
namespace Impl {

// multiplies two sparse MKL matrices and returns sparse MKL matrix
template <typename value_type>
inline static MKLSparseMatrix<value_type> mkl_spmm(
    sparse_operation_t operation, const MKLSparseMatrix<value_type> &A,
    const MKLSparseMatrix<value_type> &B) {
  sparse_matrix_t C;
  KOKKOSKERNELS_MKL_SAFE_CALL(mkl_sparse_spmm(operation, A, B, &C));
  return MKLSparseMatrix<value_type>(C);
}

template <typename KernelHandle, typename a_rowmap_view_type,
          typename a_index_view_type, typename a_values_view_type,
          typename b_rowmap_view_type, typename b_index_view_type,
          typename b_values_view_type, typename c_rowmap_view_type,
          typename c_index_view_type, typename c_values_view_type>
class MKL_SPGEMM {
 public:
  typedef typename KernelHandle::nnz_lno_t nnz_lno_t;
  typedef typename KernelHandle::size_type size_type;
  typedef typename KernelHandle::nnz_scalar_t value_type;
  typedef typename KernelHandle::HandleExecSpace MyExecSpace;
  typedef typename Kokkos::View<int *, Kokkos::HostSpace> int_tmp_view_t;

 public:
  static void mkl_symbolic(KernelHandle *handle, nnz_lno_t m, nnz_lno_t n,
                           nnz_lno_t k, a_rowmap_view_type row_mapA,
                           a_index_view_type entriesA, bool transposeA,
                           b_rowmap_view_type row_mapB,
                           b_index_view_type entriesB, bool transposeB,
                           c_rowmap_view_type row_mapC, bool verbose = false) {
    if (m < 1 || n < 1 || k < 1 || entriesA.extent(0) < 1 ||
        entriesB.extent(0) < 1) {
      // set correct values in non-empty 0-nnz corner case
      handle->set_c_nnz(0);
      Kokkos::deep_copy(row_mapC, 0);
      return;
    }

    Kokkos::Timer timer;
    using scalar_t = typename KernelHandle::nnz_scalar_t;

    const auto export_rowmap = [&](MKL_INT num_rows, MKL_INT *rows_start,
                                   MKL_INT * /*columns*/,
                                   scalar_t * /*values*/) {
      if (handle->mkl_keep_output) {
        Kokkos::Timer copy_time;
        const nnz_lno_t nnz = rows_start[num_rows];
        handle->set_c_nnz(nnz);
        copy(make_host_view(rows_start, num_rows + 1), row_mapC);
        if (verbose)
          std::cout << "\tMKL rowmap export time:" << copy_time.seconds()
                    << std::endl;
      }
    };

    // use dummy values for A and B inputs
    a_values_view_type tmp_valsA(
        Kokkos::ViewAllocateWithoutInitializing("tmp_valuesA"),
        entriesA.extent(0));
    b_values_view_type tmp_valsB(
        Kokkos::ViewAllocateWithoutInitializing("tmp_valuesB"),
        entriesB.extent(0));

    spmm(handle, m, n, k, row_mapA, entriesA, tmp_valsA, transposeA, row_mapB,
         entriesB, tmp_valsB, transposeB, verbose, export_rowmap);

    if (verbose)
      std::cout << "MKL symbolic time:" << timer.seconds() << std::endl;
  }

  static void mkl_numeric(
      KernelHandle *handle, nnz_lno_t m, nnz_lno_t n, nnz_lno_t k,
      a_rowmap_view_type row_mapA, a_index_view_type entriesA,
      a_values_view_type valuesA, bool transposeA, b_rowmap_view_type row_mapB,
      b_index_view_type entriesB, b_values_view_type valuesB, bool transposeB,
      c_rowmap_view_type /* row_mapC */, c_index_view_type entriesC,
      c_values_view_type valuesC, bool verbose = false) {
    Kokkos::Timer timer;

    const auto export_values =
        [&](MKL_INT num_rows, MKL_INT *rows_start, MKL_INT *columns,
            typename KernelHandle::nnz_scalar_t *values) {
          if (handle->mkl_keep_output) {
            Kokkos::Timer copy_time;
            const nnz_lno_t nnz = rows_start[num_rows];
            copy(make_host_view(columns, nnz), entriesC);
            copy(make_host_view(values, nnz), valuesC);
            if (verbose)
              std::cout << "\tMKL values export time:" << copy_time.seconds()
                        << std::endl;
          }
        };

    spmm(handle, m, n, k, row_mapA, entriesA, valuesA, transposeA, row_mapB,
         entriesB, valuesB, transposeB, verbose, export_values);

    if (verbose)
      std::cout << "MKL numeric time:" << timer.seconds() << std::endl;
  }

 private:
  template <typename CB>
  static void spmm(KernelHandle * /* handle */, nnz_lno_t m, nnz_lno_t n,
                   nnz_lno_t k, a_rowmap_view_type row_mapA,
                   a_index_view_type entriesA, a_values_view_type valuesA,

                   bool transposeA, b_rowmap_view_type row_mapB,
                   b_index_view_type entriesB, b_values_view_type valuesB,
                   bool transposeB, bool verbose, const CB &callback) {
    if (!std::is_same<nnz_lno_t, int>::value) {
      throw std::runtime_error("MKL requires local ordinals to be integer.\n");
    }

    if (m < 1 || n < 1 || k < 1 || entriesA.extent(0) < 1 ||
        entriesB.extent(0) < 1) {
      return;
    }

    const auto create_mirror = [](auto view) {
      return Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), view);
    };

    auto h_rowsA      = create_mirror(row_mapA);
    auto h_rowsB      = create_mirror(row_mapB);
    const int *a_xadj = reinterpret_cast<const int *>(h_rowsA.data());
    const int *b_xadj = reinterpret_cast<const int *>(h_rowsB.data());
    int_tmp_view_t a_xadj_v, b_xadj_v;

    if (!std::is_same<size_type, int>::value) {
      if (entriesA.extent(0) > INT_MAX || entriesB.extent(0) > INT_MAX) {
        throw std::runtime_error(
            "MKL requires integer values for size type for SPGEMM. Copying "
            "to "
            "integer will cause overflow.\n");
      }
      static_assert(
          std::is_same<typename int_tmp_view_t::value_type,
                       typename int_tmp_view_t::non_const_value_type>::value,
          "deep_copy requires non-const destination type");

      Kokkos::Timer copy_time;
      a_xadj_v = int_tmp_view_t("tmpa", m + 1);
      b_xadj_v = int_tmp_view_t("tmpb", n + 1);
      Kokkos::deep_copy(a_xadj_v, h_rowsA);
      Kokkos::deep_copy(b_xadj_v, h_rowsB);
      a_xadj = (int *)a_xadj_v.data();
      b_xadj = (int *)b_xadj_v.data();
      if (verbose)
        std::cout << "\tMKL int-type temp rowmap copy time:"
                  << copy_time.seconds() << std::endl;
    }

    auto h_valsA           = create_mirror(valuesA);
    auto h_valsB           = create_mirror(valuesB);
    auto h_entriesA        = create_mirror(entriesA);
    auto h_entriesB        = create_mirror(entriesB);
    const int *a_adj       = reinterpret_cast<const int *>(h_entriesA.data());
    const int *b_adj       = reinterpret_cast<const int *>(h_entriesB.data());
    const value_type *a_ew = h_valsA.data();
    const value_type *b_ew = h_valsB.data();

    // Hack: we discard const with pointer casts here to work around MKL
    // requiring mutable input and our symbolic interface not providing it
    using Matrix = MKLSparseMatrix<value_type>;
    Matrix A(m, n, (int *)a_xadj, (int *)a_adj, (value_type *)a_ew);
    Matrix B(n, k, (int *)b_xadj, (int *)b_adj, (value_type *)b_ew);

    sparse_operation_t operation;
    if (transposeA && transposeB) {
      operation = SPARSE_OPERATION_TRANSPOSE;
    } else if (!(transposeA || transposeB)) {
      operation = SPARSE_OPERATION_NON_TRANSPOSE;
    } else {
      throw std::runtime_error(
          "MKL either transpose both matrices, or none for SPGEMM\n");
    }

    Kokkos::Timer timer1;
    Matrix C = mkl_spmm(operation, A, B);
    if (verbose) {
      std::cout << "\tMKL spmm (";
      if (std::is_same<float, value_type>::value)
        std::cout << "FLOAT";
      else if (std::is_same<double, value_type>::value)
        std::cout << "DOUBLE";
      else
        std::cout << "?";
      std::cout << ") time:" << timer1.seconds() << std::endl;
    }

    MKL_INT num_rows, num_cols, *rows_start, *columns;
    value_type *values;
    C.export_data(num_rows, num_cols, rows_start, columns, values);
    callback(m, rows_start, columns, values);

    A.destroy();
    B.destroy();
    C.destroy();
  }

  template <typename from_view_type, typename dst_view_type>
  inline static void copy(from_view_type from, dst_view_type to) {
    auto h_from =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), from);
    auto h_to = Kokkos::create_mirror_view(Kokkos::HostSpace(), to);
    Kokkos::deep_copy(h_to, h_from);  // view copy (for different element types)
    Kokkos::deep_copy(to, h_to);
  }

  template <typename T,
            typename view_type = Kokkos::View<const T *, Kokkos::HostSpace>>
  inline static view_type make_host_view(const T *data, size_t num_elems) {
    return view_type(data, num_elems);
  }
};

template <typename KernelHandle, typename a_rowmap_type, typename a_index_type,
          typename b_rowmap_type, typename b_index_type, typename c_rowmap_type,
          typename nnz_lno_t = typename KernelHandle::nnz_lno_t>
void mkl_symbolic(KernelHandle *handle, nnz_lno_t m, nnz_lno_t n, nnz_lno_t k,
                  a_rowmap_type row_mapA, a_index_type entriesA,
                  bool transposeA, b_rowmap_type row_mapB,
                  b_index_type entriesB, bool transposeB,
                  c_rowmap_type row_mapC, bool verbose = false) {
  using values_type  = typename KernelHandle::scalar_temp_work_view_t;
  using c_index_type = b_index_type;
  using mkl = MKL_SPGEMM<KernelHandle, a_rowmap_type, a_index_type, values_type,
                         b_rowmap_type, b_index_type, values_type,
                         c_rowmap_type, c_index_type, values_type>;
  mkl::mkl_symbolic(handle, m, n, k, row_mapA, entriesA, transposeA, row_mapB,
                    entriesB, transposeB, row_mapC, verbose);
}

template <typename KernelHandle, typename a_rowmap_type, typename a_index_type,
          typename a_values_type, typename b_rowmap_type, typename b_index_type,
          typename b_values_type, typename c_rowmap_type, typename c_index_type,
          typename c_values_type,
          typename nnz_lno_t = typename KernelHandle::nnz_lno_t>
void mkl_numeric(KernelHandle *handle, nnz_lno_t m, nnz_lno_t n, nnz_lno_t k,
                 a_rowmap_type row_mapA, a_index_type entriesA,
                 a_values_type valuesA, bool transposeA, b_rowmap_type row_mapB,
                 b_index_type entriesB, b_values_type valuesB, bool transposeB,
                 c_rowmap_type row_mapC, c_index_type entriesC,
                 c_values_type valuesC, bool verbose = false) {
  using mkl =
      MKL_SPGEMM<KernelHandle, a_rowmap_type, a_index_type, a_values_type,
                 b_rowmap_type, b_index_type, b_values_type, c_rowmap_type,
                 c_index_type, c_values_type>;
  mkl::mkl_numeric(handle, m, n, k, row_mapA, entriesA, valuesA, transposeA,
                   row_mapB, entriesB, valuesB, transposeB, row_mapC, entriesC,
                   valuesC, verbose);
}

}  // namespace Impl
}  // namespace KokkosSparse

#endif  // KOKKOSKERNELS_ENABLE_TPL_MKL
#endif  // _KOKKOSSPGEMMMKL_HPP
