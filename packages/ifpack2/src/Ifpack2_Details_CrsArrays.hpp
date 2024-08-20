// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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

//Utility for getting the local values, rowptrs and colinds (in Kokkos::Views) for any RowMatrix
//Used by Fic, Filu and Fildl but may also be useful in other classes
template<typename Scalar, typename ImplScalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
struct CrsArrayReader
{
  typedef typename Node::device_type device_type;
  typedef typename device_type::execution_space execution_space;
  typedef Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> TRowMatrix;
  typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> TCrsMatrix;
  typedef Tpetra::BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> TBcrsMatrix;
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
    auto Abcrs = dynamic_cast<const TBcrsMatrix*>(A);
    if(Acrs)
    {
      getValuesCrs(Acrs, vals);
      return;
    }
    if (Abcrs) {
      getValuesBcrs(Abcrs, vals);
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
    auto Abcrs = dynamic_cast<const TBcrsMatrix*>(A);
    if(Acrs)
    {
      getStructureCrs(Acrs, rowptrsHost, rowptrs, colinds);
      return;
    }
    if (Abcrs) {
      getStructureBcrs(Abcrs, rowptrsHost, rowptrs, colinds);
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
    auto localA = A->getLocalMatrixDevice();
    auto values = localA.values;
    auto nnz = values.extent(0);
    values_ = ScalarArray("Values", nnz );
    Kokkos::deep_copy(values_, values);
  }

  //! Faster specialization of getStructure() for when A is a Tpetra::CrsMatrix.
  static void getStructureCrs(const TCrsMatrix* A, OrdinalArrayHost& rowptrsHost_, OrdinalArray& rowptrs_, OrdinalArray& colinds_)
  {
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
  static void getValuesBcrs(const TBcrsMatrix* A, ScalarArray& values_)
  {
    auto localA = A->getLocalMatrixDevice();
    auto values = localA.values;
    auto nnz = values.extent(0);
    values_ = ScalarArray("Values", nnz );
    Kokkos::deep_copy(values_, values);
  }

  //! Faster specialization of getStructure() for when A is a Tpetra::BlockCrsMatrix.
  static void getStructureBcrs(const TBcrsMatrix* A, OrdinalArrayHost& rowptrsHost_, OrdinalArray& rowptrs_, OrdinalArray& colinds_)
  {
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

