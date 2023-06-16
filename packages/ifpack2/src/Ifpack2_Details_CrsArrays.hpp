/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

/// \file Ifpack2_Details_CrsArrays.hpp
/// \brief Provides functions for retrieving local CRS arrays
///   (row pointers, column indices, and values) from Tpetra::RowMatrix.
///   This is used by Ifpack2's FastILU wrapper.

#ifndef __IFPACK2_CRSARRAYS_DECL_HPP__
#define __IFPACK2_CRSARRAYS_DECL_HPP__

#include <Tpetra_RowMatrix.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_KokkosCompat_DefaultNode.hpp>
#include <Tpetra_BlockCrsMatrix_Helpers_decl.hpp>
#include <KokkosSparse_CrsMatrix.hpp>
#include <Ifpack2_LocalFilter.hpp>
#include <Ifpack2_ReorderFilter.hpp>

namespace Ifpack2
{
namespace Details
{

// template <typename View1, typename View2, typename View3>
// std::vector<std::vector<typename View3::non_const_value_type>> decompress_matrix(
//   const View1& row_map,
//   const View2& entries,
//   const View3& values,
//   const int block_size)
// {
//   using size_type = typename View1::non_const_value_type;
//   using lno_t     = typename View2::non_const_value_type;
//   using scalar_t  = typename View3::non_const_value_type;

//   const size_type nrows = row_map.extent(0) - 1;
//   std::vector<std::vector<scalar_t>> result;
//   result.resize(nrows);
//   for (auto& row : result) {
//     row.resize(nrows, 0.0);
//   }

//   std::cout << "cols: " << entries.extent(0) << std::endl;

//   for (size_type row_idx = 0; row_idx < nrows; ++row_idx) {
//     const size_type row_nnz_begin = row_map(row_idx);
//     const size_type row_nnz_end   = row_map(row_idx + 1);
//     for (size_type row_nnz = row_nnz_begin; row_nnz < row_nnz_end; ++row_nnz) {
//       const lno_t col_idx      = entries(row_nnz);
//       const scalar_t value     = values.extent(0) > 0 ? values(row_nnz) : 1;
//       result[row_idx][col_idx] = value;
//     }
//   }

//   return result;
// }

// template <typename scalar_t>
// void print_matrix(const std::vector<std::vector<scalar_t>>& matrix) {
//   for (const auto& row : matrix) {
//     for (const auto& item : row) {
//       std::printf("%.2f ", item);
//     }
//     std::cout << std::endl;
//   }
// }

//Utility for getting the local values, rowptrs and colinds (in Kokkos::Views) for any RowMatrix
//Used by Fic, Filu and Fildl but may also be useful in other classes
template<typename Scalar, typename ImplScalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
struct CrsArrayReader
{
  typedef typename Node::device_type device_type;
  typedef typename device_type::execution_space execution_space;
  typedef Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> TRowMatrix;
  typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> TCrsMatrix;
  typedef Tpetra::BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> TBrsMatrix;
  typedef Ifpack2::LocalFilter<TRowMatrix> Filter;
  typedef Ifpack2::ReorderFilter<TRowMatrix> ReordFilter;
  typedef KokkosSparse::CrsMatrix<ImplScalar, LocalOrdinal, execution_space> KCrsMatrix;
  typedef Kokkos::View<LocalOrdinal*, execution_space> OrdinalArray;
  typedef Kokkos::View<ImplScalar*, execution_space> ScalarArray;
  typedef typename OrdinalArray::HostMirror  OrdinalArrayHost;
  //! The execution space to used to run the row access functors.
  typedef Kokkos::Serial functor_space;
  typedef Kokkos::RangePolicy<functor_space, int> RangePol;

  //! Get the values (in device view) of the local rows of A.
  // \param[in] A The matrix
  // \param[out] vals The values array on device (allocated inside function)
  // \param[in] rowptrs The rowptrs host view provided by getStructure()
  static void getValues(const TRowMatrix* A, ScalarArray& vals, OrdinalArrayHost& rowptrs)
  {
    auto Acrs = dynamic_cast<const TCrsMatrix*>(A);
    auto Abrs = dynamic_cast<const TBrsMatrix*>(A);
    if(Acrs)
    {
      getValuesCrs(Acrs, vals);
      return;
    }
    if (Abrs) {
      getValuesBrs(Abrs, vals);
      return;
    }
    using range_type = Kokkos::pair<int, int>;
    using local_inds_host_view_type = typename TRowMatrix::nonconst_local_inds_host_view_type;
    using values_host_view_type     = typename TRowMatrix::nonconst_values_host_view_type;
    using scalar_type               = typename values_host_view_type::value_type;

    LocalOrdinal nrows = A->getLocalNumRows();
    size_t nnz = A->getLocalNumEntries();
    size_t maxNnz = A->getLocalMaxNumRowEntries();

    vals = ScalarArray("Values", nnz);
    auto valsHost = Kokkos::create_mirror(vals);
    local_inds_host_view_type lclColInds ("lclColinds", maxNnz);

    nnz = 0;
    for(LocalOrdinal i = 0; i < nrows; i++) {
      size_t NumEntries = A->getNumEntriesInLocalRow(i);
      auto constLclValues = Kokkos::subview (valsHost, range_type (nnz, nnz+NumEntries));
      values_host_view_type lclValues (const_cast<scalar_type*>(constLclValues.data()), NumEntries);

      A->getLocalRowCopy (i, lclColInds, lclValues, NumEntries);
      nnz += NumEntries;
    }
    Kokkos::deep_copy(vals, valsHost);
  }

  //! Get the structure (rowptrs and colinds) of the local rows of A.
  // \param[in] A The matrix
  // \param[out] rowptrsHost The rowptrs array, in host space (allocated inside function)
  // \param[out] rowptrs The rowptrs host array, in device space (allocated inside function). Will have exactly the same values as rowptrsHost.
  // \param[out] colinds The colinds array, in device space (allocated inside function)
  static void getStructure(const TRowMatrix* A, OrdinalArrayHost& rowptrsHost, OrdinalArray& rowptrs, OrdinalArray& colinds)
  {
    auto Acrs = dynamic_cast<const TCrsMatrix*>(A);
    auto Abrs = dynamic_cast<const TBrsMatrix*>(A);
    if(Acrs)
    {
      getStructureCrs(Acrs, rowptrsHost, rowptrs, colinds);
      return;
    }
    if (Abrs) {
      getStructureBrs(Abrs, rowptrsHost, rowptrs, colinds);
      return;
    }
    //Need to allocate new array, then copy in one row at a time
    //note: actual rowptrs in the CrsMatrix implementation is size_t, but
    //FastILU needs them as LocalOrdinal so this function provides an OrdinalArray
    LocalOrdinal nrows = A->getLocalNumRows();
    rowptrsHost = OrdinalArrayHost("RowPtrs (host)", nrows + 1);

    using range_type = Kokkos::pair<int, int>;
    using values_host_view_type     = typename TRowMatrix::nonconst_values_host_view_type;
    using local_inds_host_view_type = typename TRowMatrix::nonconst_local_inds_host_view_type;
    using local_ind_type            = typename local_inds_host_view_type::value_type;
    size_t nnz = A->getLocalNumEntries();
    size_t maxNnz = A->getLocalMaxNumRowEntries();

    colinds = OrdinalArray("ColInds", nnz);
    auto colindsHost = Kokkos::create_mirror(colinds);
    values_host_view_type lclValues ("lclValues", maxNnz);

    nnz = 0;
    rowptrsHost[0] = nnz;
    for(LocalOrdinal i = 0; i < nrows; i++) {
      size_t NumEntries = A->getNumEntriesInLocalRow(i);
      auto constLclValues = Kokkos::subview (colindsHost, range_type (nnz, nnz+NumEntries));
      local_inds_host_view_type lclColInds (const_cast<local_ind_type*>(constLclValues.data()), NumEntries);
      A->getLocalRowCopy (i, lclColInds, lclValues, NumEntries);

      nnz += NumEntries;
      rowptrsHost[i+1] = nnz;
    }

    rowptrs = OrdinalArray("RowPtrs", nrows + 1);
    Kokkos::deep_copy(rowptrs, rowptrsHost);
    Kokkos::deep_copy(colinds, colindsHost);
  }

  private:

  //! Faster specialization of getValues() for when A is a Tpetra::CrsMatrix.
  static void getValuesCrs(const TCrsMatrix* A, ScalarArray& values_)
  {
    std::cout << "getValuesCrs" << std::endl;
    auto localA = A->getLocalMatrixDevice();
    auto values = localA.values;
    auto nnz = values.extent(0);
    values_ = ScalarArray("Values", nnz );
    Kokkos::deep_copy(values_, values);
  }

  //! Faster specialization of getStructure() for when A is a Tpetra::CrsMatrix.
  static void getStructureCrs(const TCrsMatrix* A, OrdinalArrayHost& rowptrsHost_, OrdinalArray& rowptrs_, OrdinalArray& colinds_)
  {
    std::cout << "getStructureCrs" << std::endl;
    //rowptrs really have data type size_t, but need them as LocalOrdinal, so must convert manually
    auto localA = A->getLocalMatrixDevice();
    auto rowptrs = localA.graph.row_map;
    auto colinds = localA.graph.entries;
    auto numRows = A->getLocalNumRows();
    auto nnz = colinds.extent(0);
    //allocate rowptrs, it's a deep copy (colinds is a shallow copy so not necessary for it)
    rowptrs_ = OrdinalArray("RowPtrs", numRows + 1);
    colinds_ = OrdinalArray("ColInds", nnz );
    Kokkos::deep_copy(rowptrs_, rowptrs);
    Kokkos::deep_copy(colinds_, colinds);
    // deep-copy to host
    rowptrsHost_ = Kokkos::create_mirror(rowptrs_);
    Kokkos::deep_copy(rowptrsHost_, rowptrs_);
  }

  //! Faster specialization of getValues() for when A is a Tpetra::BlockCrsMatrix.
  static void getValuesBrs(const TBrsMatrix* A, ScalarArray& values_)
  {
    std::cout << "getValuesBrs" << std::endl;

    auto localA = A->getLocalMatrixDevice();
    auto values = localA.values;
    auto nnz = values.extent(0);
    values_ = ScalarArray("Values", nnz );
    Kokkos::deep_copy(values_, values);

    // Tpetra::blockCrsMatrixWriter(*A, std::cout);
    // auto crs_matrix = Tpetra::convertToCrsMatrix(*A);
    // OrdinalArrayHost rowptrsHost;
    // OrdinalArray rowptrs;
    // OrdinalArray colinds;
    // ScalarArray values_tmp;
    // getStructureCrs(&(*crs_matrix), rowptrsHost, rowptrs, colinds);
    // getValuesCrs(&(*crs_matrix), values_tmp);

    // const auto decomp = decompress_matrix(rowptrs, colinds, values_tmp, 1);
    // print_matrix(decomp);
  }

  //! Faster specialization of getStructure() for when A is a Tpetra::BlockCrsMatrix.
  static void getStructureBrs(const TBrsMatrix* A, OrdinalArrayHost& rowptrsHost_, OrdinalArray& rowptrs_, OrdinalArray& colinds_)
  {
    std::cout << "getStructureBrs" << std::endl;
    //rowptrs really have data type size_t, but need them as LocalOrdinal, so must convert manually
    auto localA = A->getLocalMatrixDevice();
    auto rowptrs = localA.graph.row_map;
    auto colinds = localA.graph.entries;
    auto numRows = A->getLocalNumRows();
    auto nnz = colinds.extent(0);
    //allocate rowptrs, it's a deep copy (colinds is a shallow copy so not necessary for it)
    rowptrs_ = OrdinalArray("RowPtrs", numRows + 1);
    colinds_ = OrdinalArray("ColInds", nnz );
    Kokkos::deep_copy(rowptrs_, rowptrs);
    Kokkos::deep_copy(colinds_, colinds);
    // deep-copy to host
    rowptrsHost_ = Kokkos::create_mirror(rowptrs_);
    Kokkos::deep_copy(rowptrsHost_, rowptrs_);
  }

};

} //Details
} //Ifpack2

#endif

