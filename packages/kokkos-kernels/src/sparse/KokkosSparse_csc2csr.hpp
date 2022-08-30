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

#include "KokkosKernels_Utils.hpp"
#include <Kokkos_StdAlgorithms.hpp>

#ifndef _KOKKOSSPARSE_CSC2CSR_HPP
#define _KOKKOSSPARSE_CSC2CSR_HPP
namespace KokkosSparse {
namespace Impl {
template <class OrdinalType, class SizeType, class ValViewType,
          class RowIdViewType, class ColMapViewType>
class Csc2Csr {
 private:
  using CrsST             = typename ValViewType::value_type;
  using CrsOT             = OrdinalType;
  using CrsET             = typename ValViewType::execution_space;
  using CrsMT             = void;
  using CrsSzT            = SizeType;
  using CrsType           = CrsMatrix<CrsST, CrsOT, CrsET, CrsMT, CrsSzT>;
  using CrsValsViewType   = typename CrsType::values_type;
  using CrsRowMapViewType = typename CrsType::row_map_type::non_const_type;
  using CrsColIdViewType  = typename CrsType::index_type;

  OrdinalType __nrows;
  OrdinalType __ncols;
  SizeType __nnz;
  ValViewType __vals;
  RowIdViewType __row_ids;
  ColMapViewType __col_map;

  RowIdViewType __crs_row_cnt;

  CrsValsViewType __crs_vals;
  CrsRowMapViewType __crs_row_map;
  CrsRowMapViewType __crs_row_map_scratch;
  CrsColIdViewType __crs_col_ids;

 public:
  struct AlgoTags {
    struct s1RowCnt {};
    struct s2RowMap {};
    struct s3Copy {};
  };

  using s1RowCntTag = typename AlgoTags::s1RowCnt;
  using s3CopyTag   = typename AlgoTags::s3Copy;

 private:
  using TeamPolicyType = Kokkos::TeamPolicy<s3CopyTag, CrsET>;

  int __suggested_team_size, __suggested_vec_size, __league_size;

  template <class FunctorType>
  void __run(FunctorType &functor) {
    // s1RowCntTag
    {
      Kokkos::parallel_for("Csc2Csr",
                           Kokkos::RangePolicy<s1RowCntTag, CrsET>(0, __nnz),
                           functor);
      CrsET().fence();
    }
    // s2RowMapTag
    {
      namespace KE = Kokkos::Experimental;
      CrsET crsET;
      // Use exclusive scan so we can allocate the row map uninitialized and
      // avoid accessing device views on the host.
      KE::exclusive_scan(crsET, KE::cbegin(__crs_row_cnt),
                         KE::cend(__crs_row_cnt), KE::begin(__crs_row_map), 0);
      CrsET().fence();
      Kokkos::deep_copy(__crs_row_map_scratch, __crs_row_map);
      CrsET().fence();
    }
    // s3CopyTag
    {
      TeamPolicyType teamPolicy(__ncols, __suggested_team_size,
                                __suggested_vec_size);
      Kokkos::parallel_for("Csc2Csr", teamPolicy, functor);
      CrsET().fence();
    }
    // TODO: s3CopySortCompressTag
  }

 public:
  template <class MemberType>
  class __Functor {
   private:
    OrdinalType __nrows;
    OrdinalType __ncols;
    SizeType __nnz;
    ValViewType __vals;
    CrsValsViewType __crs_vals;
    RowIdViewType __row_ids;
    CrsRowMapViewType __crs_row_map;
    CrsRowMapViewType __crs_row_map_scratch;
    ColMapViewType __col_map;
    CrsColIdViewType __crs_col_ids;
    RowIdViewType __crs_row_cnt;

   public:
    __Functor(OrdinalType nrows, OrdinalType ncols, SizeType nnz,
              ValViewType vals, CrsValsViewType crs_vals, RowIdViewType row_ids,
              CrsRowMapViewType crs_row_map,
              CrsRowMapViewType crs_row_map_scratch, ColMapViewType col_map,
              CrsColIdViewType crs_col_ids, RowIdViewType crs_row_cnt)
        : __nrows(nrows),
          __ncols(ncols),
          __nnz(nnz),
          __vals(vals),
          __crs_vals(crs_vals),
          __row_ids(row_ids),
          __crs_row_map(crs_row_map),
          __crs_row_map_scratch(crs_row_map_scratch),
          __col_map(col_map),
          __crs_col_ids(crs_col_ids),
          __crs_row_cnt(crs_row_cnt){};

    KOKKOS_INLINE_FUNCTION
    void operator()(const s3CopyTag &, const MemberType &member) const {
      auto j         = member.league_rank();
      auto col_start = __col_map(j);
      auto col_len   = __col_map(j + 1) - col_start;

      Kokkos::parallel_for(
          Kokkos::TeamVectorRange(member, 0, col_len), [&](const int &k) {
            auto idx = col_start + k;
            auto i   = __row_ids(idx);
            auto crs_idx =
                Kokkos::atomic_fetch_inc(&__crs_row_map_scratch.data()[i]);
            __crs_col_ids(crs_idx) = j;
            __crs_vals(crs_idx)    = __vals(idx);
          });
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const s1RowCntTag &, const int &thread_id) const {
      Kokkos::atomic_inc(&__crs_row_cnt.data()[__row_ids(thread_id)]);
    }
  };

  Csc2Csr(OrdinalType nrows, OrdinalType ncols, SizeType nnz, ValViewType vals,
          RowIdViewType row_ids, ColMapViewType col_map, int league_size = 2)
      : __nrows(nrows),
        __ncols(ncols),
        __nnz(nnz),
        __vals(vals),
        __row_ids(row_ids),
        __col_map(col_map),
        __league_size(league_size) {
    __crs_vals = CrsValsViewType(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "__crs_vals"), nnz);
    __crs_row_map = CrsRowMapViewType(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "__crs_row_map"),
        nrows + 1);
    __crs_row_map_scratch =
        CrsRowMapViewType(Kokkos::view_alloc(Kokkos::WithoutInitializing,
                                             "__crs_row_map_scratch"),
                          nrows + 1);
    __crs_col_ids = CrsColIdViewType(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "__crs_col_ids"), nnz);

    __crs_row_cnt = RowIdViewType("__crs_row_cnt", __nrows + 1);

    __Functor<typename TeamPolicyType::member_type> functor(
        __nrows, __ncols, __nnz, __vals, __crs_vals, __row_ids, __crs_row_map,
        __crs_row_map_scratch, __col_map, __crs_col_ids, __crs_row_cnt);

    KokkosKernels::Impl::get_suggested_vector_size<int64_t, CrsET>(
        __suggested_vec_size, __nrows, __nnz);
    __suggested_team_size =
        KokkosKernels::Impl::get_suggested_team_size<TeamPolicyType>(
            functor, __suggested_vec_size);

    __run(functor);
  }

  CrsType get_csrMat() {
    return CrsType("csc2csr", __nrows, __ncols, __nnz, __crs_vals,
                   __crs_row_map, __crs_col_ids);
  }
};
}  // namespace Impl
///
/// \brief Converts a csc matrix to a CrsMatrix.
/// \tparam OrdinalType The view value type associated with the RowIdViewType
/// \tparam SizeType The type of nnz
/// \tparam ValViewType The values view type
/// \tparam RowIdViewType The row ids view type
/// \tparam ColMapViewType The column map view type
/// \param nrows The number of rows in the csc matrix
/// \param ncols The number of columns in the csc matrix
/// \param nnz The number of non-zeros in the csc matrix
/// \param vals The values view of the csc matrix
/// \param row_ids The row ids view of the csc matrix
/// \param col_map The column map view of the csc matrix
/// \return A KokkosSparse::CrsMatrix.
template <class OrdinalType, class SizeType, class ValViewType,
          class RowIdViewType, class ColMapViewType>
auto csc2csr(OrdinalType nrows, OrdinalType ncols, SizeType nnz,
             ValViewType vals, RowIdViewType row_ids, ColMapViewType col_map,
             int league_size) {
  using Csc2csrType = Impl::Csc2Csr<OrdinalType, SizeType, ValViewType,
                                    RowIdViewType, ColMapViewType>;
  Csc2csrType csc2Csr(nrows, ncols, nnz, vals, row_ids, col_map, league_size);
  return csc2Csr.get_csrMat();
}
}  // namespace KokkosSparse
#endif  //  _KOKKOSSPARSE_CSC2CSR_HPP
