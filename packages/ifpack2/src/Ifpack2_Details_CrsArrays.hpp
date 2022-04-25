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
#include <Kokkos_DefaultNode.hpp>
#include <KokkosSparse_CrsMatrix.hpp>
#include <Ifpack2_LocalFilter.hpp>
#include <Ifpack2_ReorderFilter.hpp>

namespace Ifpack2
{
namespace Details
{

//Utility for getting the local values, rowptrs and colinds (in Kokkos::Views) for any RowMatrix
//Used by Fic, Filu and Fildl but may also be useful in other classes
template<typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
struct CrsArrayReader
{
  typedef typename Node::device_type device_type;
  typedef typename device_type::execution_space execution_space;
  typedef Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> TRowMatrix;
  typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> TCrsMatrix;
  typedef Ifpack2::LocalFilter<TRowMatrix> Filter;
  typedef Ifpack2::ReorderFilter<TRowMatrix> ReordFilter;
  typedef KokkosSparse::CrsMatrix<Scalar, LocalOrdinal, execution_space> KCrsMatrix;
  typedef Kokkos::View<LocalOrdinal*, execution_space> OrdinalArray;
  typedef Kokkos::View<Scalar*, execution_space> ScalarArray;
  typedef Kokkos::View<LocalOrdinal*, Kokkos::HostSpace> OrdinalArrayHost;
  typedef Kokkos::View<Scalar*, Kokkos::HostSpace> ScalarArrayHost;
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
    if(Acrs)
    {
      getValuesCrs(Acrs, vals);
      return;
    }
    using range_type = Kokkos::pair<int, int>;
    using local_inds_host_view_type = typename TRowMatrix::local_inds_host_view_type;
    using values_host_view_type = typename TRowMatrix::values_host_view_type;

    LocalOrdinal nrows = A->getLocalNumRows();
    size_t nnz = A->getLocalNumEntries();
    vals = ScalarArray("Values", nnz);

    nnz = 0;
    for(LocalOrdinal i = 0; i < nrows; i++) {
      local_inds_host_view_type lclColInds;
      values_host_view_type lclValues;
      A->getLocalRowView (i, lclColInds, lclValues);

      auto lclValues_ = Kokkos::subview (vals, range_type (nnz, nnz+lclValues.extent(0)));
      Kokkos::deep_copy(lclValues_, lclValues);
      nnz += lclValues.extent(0);
    }
  }

  //! Get the structure (rowptrs and colinds) of the local rows of A.
  // \param[in] A The matrix
  // \param[out] rowptrsHost The rowptrs array, in host space (allocated inside function)
  // \param[out] rowptrs The rowptrs host array, in device space (allocated inside function). Will have exactly the same values as rowptrsHost.
  // \param[out] colinds The colinds array, in device space (allocated inside function)
  static void getStructure(const TRowMatrix* A, OrdinalArrayHost& rowptrsHost, OrdinalArray& rowptrs, OrdinalArray& colinds)
  {
    auto Acrs = dynamic_cast<const TCrsMatrix*>(A);
    if(Acrs)
    {
      getStructureCrs(Acrs, rowptrs, colinds);
      return;
    }
    //Need to allocate new array, then copy in one row at a time
    //note: actual rowptrs in the CrsMatrix implementation is size_t, but
    //FastILU needs them as LocalOrdinal so this function provides an OrdinalArray
    LocalOrdinal nrows = A->getLocalNumRows();
    rowptrsHost = OrdinalArrayHost("RowPtrs (host)", nrows + 1);

    using range_type = Kokkos::pair<int, int>;
    using local_inds_host_view_type = typename TRowMatrix::local_inds_host_view_type;
    auto graph = A->getGraph();
    size_t nnz = A->getLocalNumEntries();
    colinds = OrdinalArray("ColInds", nnz);

    nnz = 0;
    rowptrsHost[0] = nnz;
    for(LocalOrdinal i = 0; i < nrows; i++) {
      local_inds_host_view_type lclColInds;
      graph->getLocalRowView (i, lclColInds);

      auto lclColInds_ = Kokkos::subview (colinds, range_type (nnz, nnz+lclColInds.extent(0)));
      Kokkos::deep_copy(lclColInds_, lclColInds);
      nnz += lclColInds.extent(0);
      rowptrsHost[i+1] = nnz;
    }
    rowptrs = OrdinalArray("RowPtrs", nrows + 1);
    Kokkos::deep_copy(rowptrs, rowptrsHost);
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
  static void getStructureCrs(const TCrsMatrix* A, OrdinalArray& rowptrs_, OrdinalArray& colinds_)
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
  }
};

} //Details
} //Ifpack2

#endif

