/*
//@HEADER
// ************************************************************************
//
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER
*/
#ifndef KOKKOS_CRSMATRIX_H_
#define KOKKOS_CRSMATRIX_H_

/// \file Kokkos_CrsMatrix.hpp
/// \brief Kokkos' sparse matrix interface

// FIXME (mfh 29 Sep 2013) There should never be a reason for using
// assert() in these routines.  Routines that _could_ fail should
// either return an error code (if they are device functions) or throw
// an exception (if they are not device functions).
#include <assert.h>
#include <algorithm>

#include <Kokkos_View.hpp>
#include <Kokkos_Atomic.hpp>
#ifdef KOKKOS_HAVE_CUDA
#  include <Kokkos_Cuda.hpp>
#endif
#include <Kokkos_Macros.hpp>
#include <Kokkos_StaticCrsGraph.hpp>
#include <Kokkos_MV.hpp>

#ifndef _OPENMP
#include <Kokkos_Threads.hpp>
#endif // ! _OPENMP

#ifdef KOKKOS_USE_CUSPARSE
#  include <cusparse_v2.h>
#  include <Kokkos_CrsMatrix_CuSparse.hpp>
#endif // KOKKOS_USE_CUSPARSE

#ifdef KOKKOS_USE_MKL
#  include <mkl.h>
#  include <mkl_spblas.h>
#  include <Kokkos_CrsMatrix_MKL.hpp>
#endif // KOKKOS_USE_MKL

namespace Kokkos {

/// \class SparseRowView
/// \brief View of a row of a sparse matrix.
/// \tparam MatrixType Sparse matrix type, such as (but not limited to) CrsMatrix.
///
/// This class provides a generic view of a row of a sparse matrix.
/// We intended this class to view a row of a CrsMatrix, but
/// MatrixType need not necessarily be CrsMatrix.
///
/// The row view is suited for computational kernels like sparse
/// matrix-vector multiply, as well as for modifying entries in the
/// sparse matrix.  Whether the view is const or not, depends on
/// whether MatrixType is a const or nonconst view of the matrix.  If
/// you always want a const view, use SparseRowViewConst (see below).
///
/// Here is an example loop over the entries in the row:
/// \code
/// typedef typename SparseRowView<MatrixType>::scalar_type scalar_type;
/// typedef typename SparseRowView<MatrixType>::ordinal_type ordinal_type;
///
/// SparseRowView<MatrixType> A_i = ...;
/// const int numEntries = A_i.length;
/// for (int k = 0; k < numEntries; ++k) {
///   scalar_type A_ij = A_i.value (k);
///   ordinal_type j = A_i.colidx (k);
///   // ... do something with A_ij and j ...
/// }
/// \endcode
///
/// MatrixType must provide the \c scalar_type and \c ordinal_type
/// typedefs.  In addition, it must make sense to use SparseRowView to
/// view a row of MatrixType.  In particular, the values and column
/// indices of a row must be accessible using the <tt>values</tt>
/// resp. <tt>colidx</tt> arrays given to the constructor of this
/// class, with a constant <tt>stride</tt> between successive entries.
/// The stride is one for the compressed sparse row storage format (as
/// is used by CrsMatrix), but may be greater than one for other
/// sparse matrix storage formats (e.g., ELLPACK or jagged diagonal).
template<class MatrixType>
struct SparseRowView {
  //! The type of the values in the row.
  typedef typename MatrixType::scalar_type scalar_type;
  //! The type of the column indices in the row.
  typedef typename MatrixType::ordinal_type ordinal_type;

private:
  //! Array of values in the row.
  scalar_type* values_;
  //! Array of (local) column indices in the row.
  ordinal_type* colidx_;
  //! Stride between successive entries in the row.
  const int stride_;

public:
  /// \brief Constructor
  ///
  /// \param values [in] Array of the row's values.
  /// \param colidx [in] Array of the row's column indices.
  /// \param stride [in] (Constant) stride between matrix entries in
  ///   each of the above arrays.
  /// \param count [in] Number of entries in the row.
  KOKKOS_INLINE_FUNCTION
  SparseRowView (scalar_type* const values,
                 ordinal_type* const colidx,
                 const int stride,
                 const int count) :
    values_ (values), colidx_ (colidx), stride_ (stride), length (count)
  {}

  /// \brief Number of entries in the row.
  ///
  /// This is a public const field rather than a public const method,
  /// in order to avoid possible overhead of a method call if the
  /// compiler is unable to inline that method call.
  const int length;

  /// \brief Reference to the value of entry i in this row of the sparse matrix.
  ///
  /// "Entry i" is not necessarily the entry with column index i, nor
  /// does i necessarily correspond to the (local) row index.
  KOKKOS_INLINE_FUNCTION
  scalar_type& value (const int& i) const {
    return values_[i*stride_];
  }

  /// \brief Reference to the column index of entry i in this row of the sparse matrix.
  ///
  /// "Entry i" is not necessarily the entry with column index i, nor
  /// does i necessarily correspond to the (local) row index.
  KOKKOS_INLINE_FUNCTION
  ordinal_type& colidx (const int& i) const {
    return colidx_[i*stride_];
  }
};


/// \class SparseRowViewConst
/// \brief Const view of a row of a sparse matrix.
/// \tparam MatrixType Sparse matrix type, such as (but not limited to) CrsMatrix.
///
/// This class is like SparseRowView, except that it provides a const
/// view.  This class exists in order to let users get a const view of
/// a row of a nonconst matrix.
template<class MatrixType>
struct SparseRowViewConst {
  //! The type of the values in the row.
  typedef const typename MatrixType::nonconst_scalar_type scalar_type;
  //! The type of the column indices in the row.
  typedef const typename MatrixType::nonconst_ordinal_type ordinal_type;

private:
  //! Array of values in the row.
  scalar_type* values_;
  //! Array of (local) column indices in the row.
  ordinal_type* colidx_;
  //! Stride between successive entries in the row.
  const int stride_;

public:
  /// \brief Constructor
  ///
  /// \param values [in] Array of the row's values.
  /// \param colidx [in] Array of the row's column indices.
  /// \param stride [in] (Constant) stride between matrix entries in
  ///   each of the above arrays.
  /// \param count [in] Number of entries in the row.
  KOKKOS_INLINE_FUNCTION
  SparseRowViewConst (scalar_type* const values,
                      ordinal_type* const colidx,
                      const int stride,
                      const int count) :
    values_ (values), colidx_ (colidx), stride_ (stride), length (count)
  {}

  /// \brief Number of entries in the row.
  ///
  /// This is a public const field rather than a public const method,
  /// in order to avoid possible overhead of a method call if the
  /// compiler is unable to inline that method call.
  const int length;

  /// \brief (Const) reference to the value of entry i in this row of the sparse matrix.
  ///
  /// "Entry i" is not necessarily the entry with column index i, nor
  /// does i necessarily correspond to the (local) row index.
  KOKKOS_INLINE_FUNCTION
  scalar_type& value (const int& i) const {
    return values_[i*stride_];
  }

  /// \brief (Const) reference to the column index of entry i in this row of the sparse matrix.
  ///
  /// "Entry i" is not necessarily the entry with column index i, nor
  /// does i necessarily correspond to the (local) row index.
  KOKKOS_INLINE_FUNCTION
  ordinal_type& colidx (const int& i) const {
    return colidx_[i*stride_];
  }
};


/// \class CrsMatrix
/// \brief Compressed sparse row implementation of a sparse matrix.
/// \tparam ScalarType The type of entries in the sparse matrix.
/// \tparam OrdinalType The type of column indices in the sparse matrix.
/// \tparam Device The Kokkos Device type.
/// \tparam MemoryTraits Traits describing how Kokkos manages and
///   accesses data.  The default parameter suffices for most users.
/// 
/// "Crs" stands for "compressed row sparse."  This is the phrase
/// Trilinos traditionally uses to describe compressed sparse row
/// storage for sparse matrices, as described, for example, in Saad
/// (2nd ed.).
template<typename ScalarType,
         typename OrdinalType,
         class Device,
         class MemoryTraits = void,
         typename SizeType = size_t>

class CrsMatrix {
public:
  typedef Device      device_type;
  typedef ScalarType  scalar_type;
  typedef OrdinalType ordinal_type;
  typedef MemoryTraits memory_traits;
  typedef SizeType size_type;

  // FIXME (mfh 28 Sep 2013) Cuda::host_mirror_device_type is Threads.
  // Shouldn't CrsMatrix::host_device_type always be the same as its
  // Device's host_mirror_device_type?
  //
  // OpenMP is the default host type if you turned on OpenMP when
  // building.  OpenMP is not on by default, so if you specified in
  // the build that you wanted OpenMP, then we say that the default
  // host type is OpenMP instead of Threads.
#ifdef _OPENMP
  typedef Kokkos::OpenMP host_device_type;
#else
  typedef Kokkos::Threads host_device_type;
#endif

  //! Type of a host-memory mirror of the sparse matrix.
  typedef CrsMatrix<ScalarType, OrdinalType, host_device_type, MemoryTraits> HostMirror;

  /// \brief Type of the graph structure of the sparse matrix.
  ///
  /// FIXME (mfh 29 Sep 2013) It doesn't make much sense to use int
  /// (the fourth template parameter of CrsArray below) as SizeType,
  /// if OrdinalType is bigger than int.  We should use
  /// Kokkos::Impl::if_c to pick the default SizeType, possibly as
  /// follows:
  ///
  /// \code
  /// typedef Impl::if_c< (sizeof(OrdinalType) >= sizeof(typename ViewTraits<OrdinalType*, Kokkos::LayoutLeft, Device, void>::size_type)),
  ///   OrdinalType,
  ///   typename ViewTraits<OrdinalType*, Kokkos::LayoutLeft, Device, void>::size_type >
  /// size_type;
  /// \endcode
  ///
  /// The first argument of if_c is a bool condition.  If true,
  /// OrdinalType is size_type, else CrsArray's default size_type is
  /// size_type.  We took the ViewTraits expression from the default
  /// value of the fourth template parameter of CrsArray.  I have
  /// tested that the above expression compiles.
  ///
  /// There is also some argument to be made that size_type should be
  /// chosen dynamically, as a function of the number of entries.
  /// It's entirely possible that a (very large) local sparse matrix
  /// could have dimensions (and therefore column indices) that fit in
  /// int32_t, but more entries than can be addressed by int32_t or
  /// even uint32_t.
  ///
  /// (CRT 1 Oct 2013) approached the issue above by giving the matrix a
  /// fifth template parameter which propagates to the CrsArray (now
  /// StaticCrsGraph). This defaults to size_t. We still might look for a
  /// better solution though.

  //! Type of the graph structure of the sparse matrix.
  typedef Kokkos::StaticCrsGraph<OrdinalType, Kokkos::LayoutLeft, Device,SizeType> StaticCrsGraphType;

  //! Type of column indices in the sparse matrix.
  typedef typename StaticCrsGraphType::entries_type index_type;
  //! Type of the "row map" (which contains the offset for each row's data).
  typedef typename StaticCrsGraphType::row_map_type row_map_type;
  //! Kokkos Array type of the entries (values) in the sparse matrix.
  typedef Kokkos::View<scalar_type*, Kokkos::LayoutLeft, device_type, MemoryTraits> values_type;
  //! Const version of the type of the entries in the sparse matrix.
  typedef typename values_type::const_scalar_type  const_scalar_type;
  //! Nonconst version of the type of the entries in the sparse matrix.
  typedef typename values_type::non_const_scalar_type  non_const_scalar_type;

#ifdef KOKKOS_USE_CUSPARSE
  cusparseHandle_t cusparse_handle;
  cusparseMatDescr_t cusparse_descr;
#endif // KOKKOS_USE_CUSPARSE
  StaticCrsGraphType graph;
  values_type values;

  // FIXME (mfh 28 Sep 2013) std::vector should never appear in this
  // class, except perhaps as an input format for compatibility.

  std::vector<OrdinalType> h_entries_;
  std::vector<OrdinalType> rows_;


  /// \brief Default constructor; constructs an empty sparse matrix.
  ///
  /// FIXME (mfh 09 Aug 2013) numRows, numCols, and nnz should be
  /// properties of the graph, not the matrix.  Then CrsMatrix needs
  /// methods to get these from the graph.
  CrsMatrix() : _numRows (0), _numCols (0), _nnz (0) {}

  /// \brief Constructor that copies raw arrays of host data in
  ///   coordinate format.
  ///
  /// On input, each entry of the sparse matrix is stored in val[k],
  /// with row index rows[k] and column index cols[k].  We assume that
  /// the entries are sorted in increasing order by row index.
  ///
  /// This constructor is mainly useful for benchmarking or for
  /// reading the sparse matrix's data from a file.
  ///
  /// \param label [in] The sparse matrix's label.
  /// \param nrows [in] The number of rows.
  /// \param ncols [in] The number of columns.
  /// \param annz [in] The number of entries.
  /// \param val [in] The entries.
  /// \param rows [in] The row indices.  rows[k] is the row index of
  ///   val[k].
  /// \param cols [in] The column indices.  cols[k] is the column
  ///   index of val[k].
  /// \param pad [in] If true, pad the sparse matrix's storage with
  ///   zeros in order to improve cache alignment and / or
  ///   vectorization.
  ///
  /// FIXME (mfh 21 Jun 2013) The \c pad argument is currently not used.
  CrsMatrix (const std::string &label,
             OrdinalType nrows,
             OrdinalType ncols,
             OrdinalType annz,
             ScalarType* val,
             OrdinalType* rows,
             OrdinalType* cols,
             bool pad = false)
  {
    import (label, nrows, ncols, annz, val, rows, cols);

    // FIXME (mfh 09 Aug 2013) Specialize this on the Device type.
    // Only use cuSPARSE for the Cuda Device.
#ifdef KOKKOS_USE_CUSPARSE
    // FIXME (mfh 09 Aug 2013) This is actually static initialization
    // of the library; you should do it once for the whole program,
    // not once per matrix.  We need to protect this somehow.
    cusparseCreate (&cusparse_handle);

    // This is a per-matrix attribute.  It encapsulates things like
    // whether the matrix is lower or upper triangular, etc.  Ditto
    // for other TPLs like MKL.
    cusparseCreateMatDescr (&cusparse_descr);
#endif // KOKKOS_USE_CUSPARSE
  }

  /// \brief Constructor that accepts a row map, column indices, and values.
  ///
  /// The matrix will store and use the row map, indices, and values
  /// directly (by view, not by deep copy).
  ///
  /// \param label [in] The sparse matrix's label.
  /// \param nrows [in] The number of rows.
  /// \param ncols [in] The number of columns.
  /// \param annz [in] The number of entries.
  /// \param vals [in/out] The entries.
  /// \param rows [in/out] The row map (containing the offsets to the
  ///   data in each row).
  /// \param cols [in/out] The column indices.
  CrsMatrix (const std::string &label,
             OrdinalType nrows,
             OrdinalType ncols,
             OrdinalType annz,
             values_type vals,
             row_map_type rows,
             index_type cols) :
    _numRows (nrows),
    _numCols (ncols),
    _nnz (annz)
  {
    graph.row_map = rows;
    graph.entries = cols;
    values = vals;
#ifdef KOKKOS_USE_CUSPARSE
    cusparseCreate(&cusparse_handle);
    cusparseCreateMatDescr(&cusparse_descr);
#endif // KOKKOS_USE_CUSPARSE
  }

  /// \brief Constructor that accepts a a static graph, and values.
  ///
  /// The matrix will store and use the row map, indices, and values
  /// directly (by view, not by deep copy).
  ///
  /// \param label [in] The sparse matrix's label.
  /// \param nrows [in] The number of rows.
  /// \param ncols [in] The number of columns.
  /// \param annz [in] The number of entries.
  /// \param vals [in/out] The entries.
  /// \param rows [in/out] The row map (containing the offsets to the
  ///   data in each row).
  /// \param cols [in/out] The column indices.
  CrsMatrix (const std::string &label,
             OrdinalType ncols,
             values_type vals,
             StaticCrsGraphType graph_) :
    _numRows (graph_.row_map.dimension_0()-1),
    _numCols (ncols),
    _nnz (graph_.entries.dimension_0()),
    values(vals),
    graph(graph_)
  {
#ifdef KOKKOS_USE_CUSPARSE
    cusparseCreate(&cusparse_handle);
    cusparseCreateMatDescr(&cusparse_descr);
#endif // KOKKOS_USE_CUSPARSE
  }
  
  void
  import (const std::string &label,
          OrdinalType nrows,
          OrdinalType ncols,
          OrdinalType annz,
          ScalarType* val,
          OrdinalType* rows,
          OrdinalType* cols);

  //! This is a method only for testing that creates a random sparse matrix.
  void
  generate (const std::string &label,
            OrdinalType nrows,
            OrdinalType ncols,
            OrdinalType target_nnz,
            OrdinalType varianz_nel_row,
            OrdinalType width_row);

  void
  generate (const std::string &label,
      OrdinalType nrows,
      OrdinalType ncols,
      OrdinalType cols_per_row);

  void
  generate (const std::string &label);

  void
  generateHostGraph (OrdinalType nrows,
      OrdinalType ncols,
      OrdinalType cols_per_row);

  // FIXME (mfh 29 Sep 2013) See notes on the three-argument version
  // of this method below.
  void
  insertInGraph(OrdinalType row, OrdinalType col)
  {
    insertInGraph(row, &col, 1);
  }

  // FIXME (mfh 29 Sep 2013) There should not be an "insertInGraph"
  // method.  If you want to insert into the graph, you should get the
  // graph and insert into it.  If you want to change the structure of
  // the matrix, you should be required to specify a value to put in
  // the new spot.
  //
  // Furthermore, this should be a device function, by analogy with
  // UnorderedMap.
  void
  insertInGraph (const OrdinalType row, OrdinalType *cols, const size_t ncol)
  {
    OrdinalType* const start = &h_entries_[rows_[row]];
    OrdinalType* const end   = &h_entries_[rows_[row+1]];
    for (size_t i = 0; i < ncol; ++i) {
      OrdinalType *iter = start;
      while (iter < end && *iter != -1 && *iter != cols[i]) {
        ++iter;
      }

      // FIXME (mfh 29 Sep 2013) Use of assert() statements is only
      // acceptable for debugging.  It's legitimate for insertions to
      // fail.  We should use the techniques that Dan Sunderland uses
      // in UnorderedMap; for example:
      //
      // 1. Insertion should return an indication of success or failure
      // 2. The graph should keep track of the number of failed insertions

      assert (iter != end );
      *iter = cols[i];
    }
  }

  // FIXME (mfh 29 Sep 2013) We need a way to disable atomic updates
  // for ScalarType types that do not support them.  We're pretty much
  // limited to ScalarType = float, double, and {u}int{32,64}_t.  It
  // could make sense to do atomic add updates elementwise for complex
  // numbers, but that's about it unless we have transactional memory
  // extensions.  Dan Sunderland explained to me that the "array of
  // atomic int 'locks'" approach (for ScalarType that don't directly
  // support atomic updates) won't work on GPUs.
  KOKKOS_INLINE_FUNCTION
  void
  sumIntoValues (const OrdinalType row, 
		 const OrdinalType cols[],
		 const size_t ncol, 
		 ScalarType vals[], 
		 const bool force_atomic = false) const 
  {
    SparseRowView<CrsMatrix> row_view = this->row (row);
    const int length = row_view.length;
    for (size_t i = 0; i < ncol; ++i) {
      for (int j = 0; j < length; ++j) {
        if (row_view.colidx(j) == cols[i]) {
          if (force_atomic) {
            atomic_fetch_add(&row_view.value(j), vals[i]);
	  } else {
            row_view.value(j) += vals[i];
	  }
        }
      }
    }
  }

  // FIXME (mfh 29 Sep 2013) See above notes on sumIntoValues.
  KOKKOS_INLINE_FUNCTION
  void
  replaceValues (const OrdinalType row, 
		 const OrdinalType cols[], 
		 const size_t ncol, 
		 ScalarType vals[], 
		 const bool force_atomic = false) const 
  {
    SparseRowView<CrsMatrix> row_view = this->row (row);
    const int length = row_view.length;
    for (size_t i = 0; i < ncol; ++i) {
      for (int j = 0; j < length; ++j) {
        if (row_view.colidx(j) == cols[i]) {
          if (force_atomic) {
            atomic_exchange(&row_view.value(j), vals[i]);
	  } else {
            row_view.value(j) = vals[i];
	  }
        }
      }
    }
  }

  // FIXME (mfh 29 Sep 2013) It doesn't really make sense to template
  // on the scalar or ordinal types of the input, since direct
  // assignment of the underlying Views (which is how this operator
  // works) won't work if the types aren't compatible.  It would make
  // more sense to template on things like the Device and
  // MemoryTraits.
  //
  // COMMENT: (CRT 1 Oct 2013) the alternative it to template on the incoming type,
  // But that way it still matches this function even if you try to assign a vector
  // to a matrix. Using the same scalar type and ordinal type as the 'this' matrix does
  // not necessaryily work because of things like const / non-const

  template<typename aScalarType, typename aOrdinalType, class aDevice, class aMemoryTraits,typename aSizeType>
  CrsMatrix&
  operator= (const CrsMatrix<aScalarType,aOrdinalType,aDevice,aMemoryTraits, aSizeType>& mtx)
  {
    _numRows = mtx.numRows();
    _numCols = mtx.numCols();
    _nnz = mtx.nnz();
    graph = mtx.graph;
    values = mtx.values;
    return *this;
  }

  //! The number of rows in the sparse matrix.
  KOKKOS_INLINE_FUNCTION
  ordinal_type numRows() const {
    return _numRows;
  }
  //! The number of columns in the sparse matrix.
  KOKKOS_INLINE_FUNCTION
  ordinal_type numCols() const {
    return _numCols;
  }
  //! The number of stored entries in the sparse matrix.
  KOKKOS_INLINE_FUNCTION
  ordinal_type nnz() const {
    return _nnz;
  }

  friend struct SparseRowView<CrsMatrix>;

  //! Return a view of row i of the matrix.
  KOKKOS_INLINE_FUNCTION
  SparseRowView<CrsMatrix> row (int i) const {
    const int start = graph.row_map(i);
    const int end = graph.row_map(i+1);
    return SparseRowView<CrsMatrix> (&values(start), &graph.entries(start), 1, end - start);
  }

  //! Return a const view of row i of the matrix.
  KOKKOS_INLINE_FUNCTION
  SparseRowViewConst<CrsMatrix> rowConst (int i) const {
    const int start = graph.row_map(i);
    const int end = graph.row_map(i+1);
    return SparseRowViewConst<CrsMatrix> (&values(start), &graph.entries(start), 1, end - start);
  }

private:
  ordinal_type _numRows;
  ordinal_type _numCols;
  ordinal_type _nnz;
};

template< typename ScalarType , typename OrdinalType, class Device, class MemoryTraits, typename SizeType >
void
CrsMatrix<ScalarType , OrdinalType, Device, MemoryTraits, SizeType >::
import (const std::string &label,
        OrdinalType nrows,
        OrdinalType ncols,
        OrdinalType annz,
        ScalarType* val,
        OrdinalType* rows,
        OrdinalType* cols)
{
  std::string str = label;
  values = values_type (str.append (".values"), annz);

  _numRows = nrows;
  _numCols = ncols;
  _nnz = annz;

  // FIXME (09 Aug 2013) CrsArray only takes std::vector for now.
  // We'll need to fix that.
  std::vector<int> row_lengths (_numRows, 0);

  // FIXME (mfh 21 Jun 2013) This calls for a parallel_for kernel.
  for (int i = 0; i < _numRows; ++i) {
    row_lengths[i] = rows[i + 1] - rows[i];
  }

  str = label;
  graph = Kokkos::create_staticcrsgraph<StaticCrsGraphType> (str.append (".graph"), row_lengths);
  typename values_type::HostMirror h_values = Kokkos::create_mirror_view (values);
  typename index_type::HostMirror h_entries = Kokkos::create_mirror_view (graph.entries);

  // FIXME (mfh 21 Jun 2013) This needs to be a parallel copy.
  // Furthermore, why are the arrays copied twice? -- once here, to a
  // host view, and once below, in the deep copy?
  for (OrdinalType i = 0; i < _nnz; ++i) {
    if (val) {
      h_values(i) = val[i];
    }
    h_entries(i) = cols[i];
  }

  Kokkos::deep_copy (values, h_values);
  Kokkos::deep_copy (graph.entries, h_entries);
}

template<typename ScalarType, typename OrdinalType, class Device, class MemoryTraits, typename SizeType>
void
CrsMatrix<ScalarType, OrdinalType, Device, MemoryTraits, SizeType >::
generate (const std::string &label,
          OrdinalType nrows,
          OrdinalType ncols,
          OrdinalType target_nnz,
          OrdinalType varianz_nel_row,
          OrdinalType width_row)
{
  _numRows = nrows;
  _numCols = ncols;

  graph.row_map = row_map_type ("CrsMatrix::rowPtr", nrows + 1);
  typename row_map_type::HostMirror h_row_map = Kokkos::create_mirror_view (graph.row_map);

  // FIXME (mfh 21 Jun 2013) What is this method actualy doing?  It
  // looks like it's not setting the structure or values of the matrix
  // at all.

  OrdinalType elements_per_row = target_nnz / nrows;
  srand(13721);
  h_row_map(0) = 0;

  for (int row = 0; row < nrows; ++row) {
   // int varianz = (1.0 * rand() / INT_MAX - 0.5) * varianz_nel_row;
   // h_row_map(row + 1) = h_row_map(row) + elements_per_row + varianz;
  }

  _nnz = h_row_map(nrows);
  values = values_type("CrsMatrix::values", _nnz);
  graph.entries = index_type("CrsMatrix::colInd", _nnz);
  typename values_type::HostMirror h_values = Kokkos::create_mirror_view(values);
  typename index_type::HostMirror h_entries = Kokkos::create_mirror_view(graph.entries);

  for(int row = 0; row < nrows; row++) {
    for(int k = h_row_map(row); k < h_row_map(row + 1); row++) {
      //int pos = (1.0 * rand() / INT_MAX - 0.5) * width_row;

      //if(pos < 0) pos += ncols;

     // if(pos >= ncols) pos -= ncols;

     // h_entries(k) = pos;
     // h_values(k) = 100.0 * rand() / INT_MAX - 50.0;
    }
  }

  Kokkos::deep_copy(values, h_values);
  Kokkos::deep_copy(graph.entries, h_entries);
  Kokkos::deep_copy(graph.row_map, h_row_map);

}

template<typename ScalarType, typename OrdinalType, class Device, class MemoryTraits, typename SizeType>
void
CrsMatrix<ScalarType, OrdinalType, Device, MemoryTraits, SizeType >::
generate (const std::string &label,
    OrdinalType nrows,
    OrdinalType ncols,
    OrdinalType cols_per_row)
{
  _numRows = nrows;
  _numCols = ncols;
  _nnz = nrows*cols_per_row;

  std::string str = label;
  values = values_type (str.append (".values"), _nnz);


  std::vector<int> row_lengths (_numRows, 0);

  // FIXME (mfh 21 Jun 2013) This calls for a parallel_for kernel.
  for (int i = 0; i < _numRows; ++i) {
    row_lengths[i] = cols_per_row;
  }

  str = label;
  graph = Kokkos::create_staticcrsgraph<StaticCrsGraphType> (str.append (".graph"), row_lengths);
  typename values_type::HostMirror h_values = Kokkos::create_mirror_view (values);
  typename index_type::HostMirror h_entries = Kokkos::create_mirror_view (graph.entries);

  // FIXME (mfh 21 Jun 2013) Why is this copy not a parallel copy?
  // Furthermore, why are the arrays copied twice? -- once here, to a
  // host view, and once below, in the deep copy?
  for (OrdinalType i = 0; i < _nnz; ++i) {
    h_values(i) = ScalarType();
    h_entries(i) = OrdinalType();
  }

  Kokkos::deep_copy (values, h_values);
  Kokkos::deep_copy (graph.entries, h_entries);
}
template<typename ScalarType, typename OrdinalType, class Device, class MemoryTraits, typename SizeType>
void
CrsMatrix<ScalarType, OrdinalType, Device, MemoryTraits, SizeType >::
generate (const std::string &label)
{
  // Compress the entries
  size_t ptr_from= 0, ptr_to = 0;
  int cur_row = 0;
  while ( ptr_from < h_entries_.size()) {
    size_t  row_stop = rows_[cur_row+1];
    while (ptr_from < row_stop) {
      if ( h_entries_[ptr_from] == OrdinalType(-1)) {
        ptr_from = row_stop;
      } else {
        h_entries_[ptr_to++] = h_entries_[ptr_from++];
      }
    }
    rows_[++cur_row] = ptr_to;
  }
  OrdinalType nrows = rows_.size()-1;
  OrdinalType nnz = ptr_to;

  h_entries_.resize(nnz);

  //sort the rows
  for (OrdinalType i=0; i<nrows; ++i )
    std::sort(&h_entries_[i], &h_entries_[i+1]);

  // generate the matrix
  import(label, nrows, nrows, nnz, NULL, &rows_[0], &h_entries_[0]);

}


template<typename ScalarType, typename OrdinalType, class Device, class MemoryTraits, typename SizeType>
void
CrsMatrix<ScalarType, OrdinalType, Device, MemoryTraits, SizeType >::
generateHostGraph ( OrdinalType nrows,
    OrdinalType ncols,
    OrdinalType cols_per_row)
{
  _numRows = nrows;
  _numCols = ncols;
  _nnz = nrows*cols_per_row;

  h_entries_.resize(_nnz, OrdinalType(-1));
  rows_.resize(_numRows+1);
  rows_[0] = 0;
  for (OrdinalType i = 0; i < _numRows; ++i)
    rows_[i+1] = rows_[i]+cols_per_row;

}

// FIXME (mfh 09 Aug 2013) These "shuffle" operations need to move
// into kokkos/core, because they are fundamental to Kokkos and not
// specific to sparse matrices.
//
// Shuffle only makes sense on >= Kepler GPUs; it doesn't work on CPUs
// or other GPUs.  We provide a generic definition (which is trivial
// and doesn't do what it claims to do) because we don't actually use
// this function unless we are on a suitable GPU, with a suitable
// Scalar type.  (For example, in the mat-vec, the "ThreadsPerRow"
// internal parameter depends both on the Device and the Scalar type,
// and it controls whether shfl_down() gets called.)
template<typename Scalar>
KOKKOS_INLINE_FUNCTION
Scalar shfl_down(const Scalar &val, const int& delta, const int& width){
  return val;
}

template<>
KOKKOS_INLINE_FUNCTION
unsigned int shfl_down<unsigned int>(const unsigned int &val, const int& delta, const int& width){
#ifdef __CUDA_ARCH__
  #if (__CUDA_ARCH__ >= 300)
    unsigned int tmp1 = val;
    int tmp = *reinterpret_cast<int*>(&tmp1);
    tmp = __shfl_down(tmp,delta,width);
    return *reinterpret_cast<unsigned int*>(&tmp);
  #else
    return val;
  #endif
#else
  return val;
#endif
}

template<>
KOKKOS_INLINE_FUNCTION
int shfl_down<int>(const int &val, const int& delta, const int& width){
#ifdef __CUDA_ARCH__
  #if (__CUDA_ARCH__ >= 300)
    return __shfl_down(val,delta,width);
  #else
    return val;
  #endif
#else
  return val;
#endif
}

template<>
KOKKOS_INLINE_FUNCTION
float shfl_down<float>(const float &val, const int& delta, const int& width){
#ifdef __CUDA_ARCH__
  #if (__CUDA_ARCH__ >= 300)
    return __shfl_down(val,delta,width);
  #else
    return val;
  #endif
#else
  return val;
#endif
}

template<>
KOKKOS_INLINE_FUNCTION
double shfl_down<double>(const double &val, const int& delta, const int& width){
#ifdef __CUDA_ARCH__
  #if (__CUDA_ARCH__ >= 300)
    int lo = __double2loint(val);
    int hi = __double2hiint(val);
    lo = __shfl_down(lo,delta,width);
    hi = __shfl_down(hi,delta,width);
    return __hiloint2double(hi,lo);
  #else
    return val;
  #endif
#else
  return val;
#endif
}


template<class RangeVector,
         class CrsMatrix,
         class DomainVector,
         class CoeffVector1,
         class CoeffVector2,
         int doalpha,
         int dobeta,
         int ThreadsPerRow>
struct MV_MultiplyFunctor {
  typedef typename CrsMatrix::device_type                   device_type ;
  typedef typename CrsMatrix::ordinal_type                    size_type ;
  typedef typename CrsMatrix::non_const_scalar_type         scalar_type ;
  typedef typename Kokkos::View<scalar_type*, typename CrsMatrix::device_type> range_values;

  CoeffVector1 beta;
  CoeffVector2 alpha;
  CrsMatrix  m_A ;
  DomainVector  m_x ;
  RangeVector  m_y ;
  size_type n;

  //--------------------------------------------------------------------------

  template<int UNROLL>
  KOKKOS_INLINE_FUNCTION
  void strip_mine (const size_type i, const size_type kk) const {
    const size_type iRow = i / ThreadsPerRow;
    const int lane = i % ThreadsPerRow;
    scalar_type sum[UNROLL];
    // FIXME (mfh 29 Sep 2013) These pragmas ("ivdep", "unroll", and
    // "loop count") should be protected by macros that identify the
    // compilers which support them.
    //
    // FIXME (mfh 29 Sep 2013) Is a "loop count" of "24" a good idea
    // when the matrix has many fewer entries per row than that?
#pragma ivdep
#pragma unroll
    for (size_type k = 0 ; k < UNROLL ; ++k) {
      // NOTE (mfh 09 Aug 2013) This requires that assignment from int
      // (in this case, 0) to scalar_type be defined.  It's not for
      // types like arprec and dd_real.
      //
      // mfh 29 Sep 2013: On the other hand, arprec and dd_real won't
      // work on CUDA devices anyway, since their methods aren't
      // device functions.  arprec has other issues (e.g., dynamic
      // memory allocation, and the array-of-structs memory layout
      // which is unfavorable to GPUs), but could be dealt with in the
      // same way as Sacado's AD types.
      sum[k] = 0;
    }

    if (doalpha != -1) {
      const SparseRowView<CrsMatrix> row = m_A.row(iRow);

#pragma ivdep
#pragma loop count (15)
#pragma unroll
      for (size_type iEntry = lane; iEntry < row.length; iEntry += ThreadsPerRow) {
        const scalar_type val = row.value(iEntry);
        const size_type ind = row.colidx(iEntry);

#pragma unroll
        for (size_type k = 0; k < UNROLL; ++k) {
          sum[k] +=  val * m_x(ind, kk + k);
        }
      }
    } else {
      const SparseRowView<CrsMatrix> row = m_A.row(iRow);

#pragma ivdep
#pragma loop count (15)
#pragma unroll
      for(size_type iEntry = lane ; iEntry < row.length ; iEntry+=ThreadsPerRow) {
        const scalar_type val = row.value(iEntry);
        const size_type ind = row.colidx(iEntry);

#pragma unroll
        for (size_type k = 0; k < UNROLL; ++k) {
          sum[k] -= val * m_x(ind, kk + k);
        }
      }
    }
    for (int i=0; i < UNROLL; ++i) {
      if (ThreadsPerRow > 1)
        sum[i] += shfl_down(sum[i], 1,ThreadsPerRow);
      if (ThreadsPerRow > 2)
        sum[i] += shfl_down(sum[i], 2,ThreadsPerRow);
      if (ThreadsPerRow > 4)
        sum[i] += shfl_down(sum[i], 4,ThreadsPerRow);
      if (ThreadsPerRow > 8)
        sum[i] += shfl_down(sum[i], 8,ThreadsPerRow);
      if (ThreadsPerRow > 16)
        sum[i] += shfl_down(sum[i], 16,ThreadsPerRow);
    }
    if (lane==0) {
      if(doalpha * doalpha != 1) {
#pragma ivdep
#pragma unroll
        for (size_type k = 0; k < UNROLL; ++k) {
          sum[k] *= alpha(kk + k);
        }
      }

      if (dobeta == 0) {
#pragma ivdep
#pragma unroll
        for (size_type k = 0; k < UNROLL; ++k) {
          m_y(iRow, kk + k) = sum[k];
        }
      } else if(dobeta == 1) {
#pragma ivdep
#pragma unroll
        for (size_type k = 0; k < UNROLL; ++k) {
          m_y(iRow, kk + k) += sum[k];
        }
      } else if (dobeta == -1) {
#pragma ivdep
#pragma unroll
        for (size_type k = 0; k < UNROLL; ++k) {
          m_y(iRow, kk + k) = -m_y(iRow, kk + k) +  sum[k];
        }
      } else {
#pragma ivdep
#pragma unroll
        for (size_type k = 0; k < UNROLL; ++k) {
          m_y(iRow, kk + k) = beta(kk + k) * m_y(iRow, kk + k) + sum[k] ;
        }
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void strip_mine_1 (const size_type i) const {
    const size_type iRow = i/ThreadsPerRow;
    const int lane = i%ThreadsPerRow;
    scalar_type sum = 0;

    if(doalpha != -1) {
      const SparseRowView<CrsMatrix> row = m_A.row(iRow);

#pragma loop count (15)
      for(size_type iEntry = lane; iEntry < row.length; iEntry += ThreadsPerRow) {
        sum += row.value(iEntry) * m_x(row.colidx(iEntry),0);
      }
    } else {
      const SparseRowView<CrsMatrix> row = m_A.row(iRow);

#pragma loop count (15)
      for(size_type iEntry = lane; iEntry < row.length; iEntry += ThreadsPerRow) {
        sum -= row.value(iEntry) * m_x(row.colidx(iEntry),0);
      }
    }

    if(ThreadsPerRow > 1)
      sum += shfl_down(sum, 1,ThreadsPerRow);
    if(ThreadsPerRow > 2)
      sum += shfl_down(sum, 2,ThreadsPerRow);
    if(ThreadsPerRow > 4)
      sum += shfl_down(sum, 4,ThreadsPerRow);
    if(ThreadsPerRow > 8)
      sum += shfl_down(sum, 8,ThreadsPerRow);
    if(ThreadsPerRow > 16)
      sum += shfl_down(sum, 16,ThreadsPerRow);

    if (lane == 0) {
      if(doalpha * doalpha != 1) {
        sum *= alpha(0);
      }

      if (dobeta == 0) {
        m_y(iRow, 0) = sum ;
      } else if (dobeta == 1) {
        m_y(iRow, 0) += sum ;
      } else if (dobeta == -1) {
        m_y(iRow, 0) = -m_y(iRow, 0) +  sum;
      } else {
        m_y(iRow, 0) = beta(0) * m_y(iRow, 0) + sum;
      }
    }
  }


  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type iRow) const {
    size_type kk = 0;

#ifdef KOKKOS_FAST_COMPILE
    for (; kk + 4 <= n; kk += 4) {
      strip_mine<4>(iRow, kk);
    }
    for( ; kk < n; ++kk) {
      strip_mine<1>(iRow, kk);
    }
#else
#  ifdef __CUDA_ARCH__
    if ((n > 8) && (n % 8 == 1)) {
      strip_mine<9>(iRow, kk);
      kk += 9;
    }
    for(; kk + 8 <= n; kk += 8)
      strip_mine<8>(iRow, kk);
    if(kk < n)
      switch(n - kk) {
#  else // NOT a CUDA device
        if ((n > 16) && (n % 16 == 1)) {
          strip_mine<17>(iRow, kk);
          kk += 17;
        }

        for (; kk + 16 <= n; kk += 16) {
          strip_mine<16>(iRow, kk);
        }

        if(kk < n)
          switch(n - kk) {
          case 15:
            strip_mine<15>(iRow, kk);
            break;

          case 14:
            strip_mine<14>(iRow, kk);
            break;

          case 13:
            strip_mine<13>(iRow, kk);
            break;

          case 12:
            strip_mine<12>(iRow, kk);
            break;

          case 11:
            strip_mine<11>(iRow, kk);
            break;

          case 10:
            strip_mine<10>(iRow, kk);
            break;

          case 9:
            strip_mine<9>(iRow, kk);
            break;

          case 8:
            strip_mine<8>(iRow, kk);
            break;
#  endif // __CUDA_ARCH__
          case 7:
            strip_mine<7>(iRow, kk);
            break;

          case 6:
            strip_mine<6>(iRow, kk);
            break;

          case 5:
            strip_mine<5>(iRow, kk);
            break;

          case 4:
            strip_mine<4>(iRow, kk);
            break;

          case 3:
            strip_mine<3>(iRow, kk);
            break;

          case 2:
            strip_mine<2>(iRow, kk);
            break;

          case 1:
            strip_mine_1(iRow);
            break;
          }
#endif // KOKKOS_FAST_COMPILE
      }
  };

  template<class RangeVector,
	   class CrsMatrix,
	   class DomainVector,
	   class CoeffVector1,
	   class CoeffVector2,
	   int doalpha,
	   int dobeta,
	   int ThreadsPerRow>
  struct MV_MultiplySingleFunctor {
    typedef typename CrsMatrix::device_type                   device_type ;
    typedef typename CrsMatrix::ordinal_type                    size_type ;
    typedef typename CrsMatrix::non_const_scalar_type         scalar_type ;
    typedef typename Kokkos::View<scalar_type*, typename CrsMatrix::device_type> range_values;

    CoeffVector1 beta;
    CoeffVector2 alpha;
    CrsMatrix  m_A ;
    DomainVector  m_x ;
    RangeVector  m_y ;
    size_type n;

    KOKKOS_INLINE_FUNCTION
    void operator()(const size_type i) const {
      const size_type iRow = i/ThreadsPerRow;
      const int lane = i%ThreadsPerRow;
      scalar_type sum = 0;

      if (doalpha != -1) {
	const SparseRowView<CrsMatrix> row = m_A.row(iRow);

#pragma loop count (15)
#pragma unroll
	for (size_type iEntry = lane; iEntry < row.length; iEntry += ThreadsPerRow) {
	  sum += row.value(iEntry) * m_x(row.colidx(iEntry));
	}
      } else {
	const SparseRowView<CrsMatrix> row = m_A.row(iRow);

#pragma loop count (15)
#pragma unroll
	for (size_type iEntry = lane; iEntry < row.length; iEntry += ThreadsPerRow) {
	  sum -= row.value(iEntry) * m_x(row.colidx(iEntry));
	}
      }
      if (ThreadsPerRow > 1)
	sum += shfl_down(sum, 1,ThreadsPerRow);
      if (ThreadsPerRow > 2)
	sum += shfl_down(sum, 2,ThreadsPerRow);
      if (ThreadsPerRow > 4)
	sum += shfl_down(sum, 4,ThreadsPerRow);
      if (ThreadsPerRow > 8)
	sum += shfl_down(sum, 8,ThreadsPerRow);
      if (ThreadsPerRow > 16)
	sum += shfl_down(sum, 16,ThreadsPerRow);

      if (lane == 0) {
	if (doalpha * doalpha != 1) {
	  sum *= alpha(0);
	}

	if (dobeta == 0) {
	  m_y(iRow) = sum ;
	} else if (dobeta == 1) {
	  m_y(iRow) += sum ;
	} else if (dobeta == -1) {
	  m_y(iRow) = -m_y(iRow) +  sum;
	} else {
	  m_y(iRow) = beta(0) * m_y(iRow) + sum;
	}
      }
    }
  };

  namespace Impl {

    template <class RangeVector,
	      class CrsMatrix,
	      class DomainVector,
	      class CoeffVector1,
	      class CoeffVector2>
    void
    MV_Multiply_Check_Compatibility (const CoeffVector1 &betav,
				     const RangeVector &y,
				     const CoeffVector2 &alphav,
				     const CrsMatrix &A,
				     const DomainVector &x,
				     const int& doalpha,
				     const int& dobeta)
    {
      typename DomainVector::size_type numVecs = x.dimension_1();
      typename DomainVector::size_type numRows = A.numRows();
      typename DomainVector::size_type numCols = A.numCols();

      if (y.dimension_1() != numVecs) {
	std::ostringstream msg;
	msg << "Error in CRSMatrix - Vector Multiply (y = by + aAx): 2nd dimensions of y and x do not match\n";
	msg << "\t Labels are: y(" << RangeVector::memory_space::query_label(y.ptr_on_device()) << ") b("
	    << CoeffVector1::memory_space::query_label(betav.ptr_on_device()) << ") a("
	    << CoeffVector2::memory_space::query_label(alphav.ptr_on_device()) << ") x("
	    << CrsMatrix::values_type::memory_space::query_label(A.values.ptr_on_device()) << ") x("
	    << DomainVector::memory_space::query_label(x.ptr_on_device()) << ")\n";
	msg << "\t Dimensions are: y(" << y.dimension_0() << "," << y.dimension_1() << ") x(" << x.dimension_0() << "," << x.dimension_1() << ")\n";
	Impl::throw_runtime_exception( msg.str() );
      }
      if (numRows > y.dimension_0()) {
	std::ostringstream msg;
	msg << "Error in CRSMatrix - Vector Multiply (y = by + aAx): dimensions of y and A do not match\n";
	msg << "\t Labels are: y(" << RangeVector::memory_space::query_label(y.ptr_on_device()) << ") b("
	    << CoeffVector1::memory_space::query_label(betav.ptr_on_device()) << ") a("
	    << CoeffVector2::memory_space::query_label(alphav.ptr_on_device()) << ") x("
	    << CrsMatrix::values_type::memory_space::query_label(A.values.ptr_on_device()) << ") x("
	    << DomainVector::memory_space::query_label(x.ptr_on_device()) << ")\n";
	msg << "\t Dimensions are: y(" << y.dimension_0() << "," << y.dimension_1() << ") A(" << A.numCols() << "," << A.numRows() << ")\n";
	Impl::throw_runtime_exception( msg.str() );
      }
      if (numCols > x.dimension_0()) {
	std::ostringstream msg;
	msg << "Error in CRSMatrix - Vector Multiply (y = by + aAx): dimensions of x and A do not match\n";
	msg << "\t Labels are: y(" << RangeVector::memory_space::query_label(y.ptr_on_device()) << ") b("
	    << CoeffVector1::memory_space::query_label(betav.ptr_on_device()) << ") a("
	    << CoeffVector2::memory_space::query_label(alphav.ptr_on_device()) << ") x("
	    << CrsMatrix::values_type::memory_space::query_label(A.values.ptr_on_device()) << ") x("
	    << DomainVector::memory_space::query_label(x.ptr_on_device()) << ")\n";
	msg << "\t Dimensions are: x(" << x.dimension_0() << "," << x.dimension_1() << ") A(" << A.numCols() << "," << A.numRows() << ")\n";
	Impl::throw_runtime_exception( msg.str() );
      }
      if (dobeta==2) {
	if (betav.dimension_0()!=numVecs) {
	  std::ostringstream msg;
	  msg << "Error in CRSMatrix - Vector Multiply (y = by + aAx): 2nd dimensions of y and b do not match\n";
	  msg << "\t Labels are: y(" << RangeVector::memory_space::query_label(y.ptr_on_device()) << ") b("
	      << CoeffVector1::memory_space::query_label(betav.ptr_on_device()) << ") a("
	      << CoeffVector2::memory_space::query_label(alphav.ptr_on_device()) << ") x("
	      << CrsMatrix::values_type::memory_space::query_label(A.values.ptr_on_device()) << ") x("
	      << DomainVector::memory_space::query_label(x.ptr_on_device()) << ")\n";
	  msg << "\t Dimensions are: y(" << y.dimension_0() << "," << y.dimension_1() << ") b(" << betav.dimension_0() << ")\n";
	  Impl::throw_runtime_exception( msg.str() );
	}
      }
      if(doalpha==2) {
	if(alphav.dimension_0()!=numVecs) {
	  std::ostringstream msg;
	  msg << "Error in CRSMatrix - Vector Multiply (y = by + aAx): 2nd dimensions of x and b do not match\n";
	  msg << "\t Labels are: y(" << RangeVector::memory_space::query_label(y.ptr_on_device()) << ") b("
	      << CoeffVector1::memory_space::query_label(betav.ptr_on_device()) << ") a("
	      << CoeffVector2::memory_space::query_label(alphav.ptr_on_device()) << ") x("
	      << CrsMatrix::values_type::memory_space::query_label(A.values.ptr_on_device()) << ") x("
	      << DomainVector::memory_space::query_label(x.ptr_on_device()) << ")\n";
	  msg << "\t Dimensions are: x(" << x.dimension_0() << "," << x.dimension_1() << ") b(" << betav.dimension_0() << ")\n";
	  Impl::throw_runtime_exception( msg.str() );
	}
      }
    }
  } // namespace Impl

  template<class Device, class ScalarType>
  struct ThreadsPerRow {
    static const int value = 1;
  };

#ifdef KOKKOS_HAVE_CUDA

  template<>
  struct ThreadsPerRow<Kokkos::Cuda,unsigned int> {
    static const int value = 8;
  };

  template<>
  struct ThreadsPerRow<Kokkos::Cuda,int> {
    static const int value = 8;
  };

  template<>
  struct ThreadsPerRow<Kokkos::Cuda,float> {
    static const int value = 8;
  };

  template<>
  struct ThreadsPerRow<Kokkos::Cuda,double> {
    static const int value = 8;
  };

#endif // KOKKOS_HAVE_CUDA

  template <class RangeVector,
	    class TCrsMatrix,
	    class DomainVector,
	    class CoeffVector1,
	    class CoeffVector2,
	    int doalpha,
	    int dobeta>
  void
  MV_Multiply (typename Kokkos::Impl::enable_if<DomainVector::Rank == 2, const CoeffVector1>::type& betav,
	       const RangeVector &y,
	       const CoeffVector2 &alphav,
	       const TCrsMatrix &A,
	       const DomainVector &x)
  {
    if (doalpha == 0) {
      if (dobeta==2) {
	      MV_MulScalar(y,betav,y);
      } else {
	      MV_MulScalar(y,dobeta,y);
      }
      return;
    } else {
      typedef View< typename RangeVector::non_const_data_type ,
                    typename RangeVector::array_layout ,
                    typename RangeVector::device_type ,
                    typename RangeVector::memory_traits >
      RangeVectorType;

      typedef View< typename DomainVector::const_data_type ,
                    typename DomainVector::array_layout ,
                    typename DomainVector::device_type ,
                    Kokkos::MemoryRandomRead >
      DomainVectorType;

      typedef View< typename CoeffVector1::const_data_type ,
                    typename CoeffVector1::array_layout ,
                    typename CoeffVector1::device_type ,
                    Kokkos::MemoryRandomRead >
      CoeffVector1Type;

      typedef View< typename CoeffVector2::const_data_type ,
                    typename CoeffVector2::array_layout ,
                    typename CoeffVector2::device_type ,
                    Kokkos::MemoryRandomRead >
      CoeffVector2Type;

      typedef CrsMatrix<typename TCrsMatrix::const_scalar_type,
                        typename TCrsMatrix::ordinal_type,
                        typename TCrsMatrix::device_type,
                        typename TCrsMatrix::memory_traits,
                        typename TCrsMatrix::size_type> CrsMatrixType;

      Impl::MV_Multiply_Check_Compatibility(betav,y,alphav,A,x,doalpha,dobeta);

      #ifndef KOKKOS_FAST_COMPILE
      MV_MultiplyFunctor<RangeVectorType, CrsMatrixType, DomainVectorType,
                         CoeffVector1Type, CoeffVector2Type, doalpha, dobeta,
	          ThreadsPerRow<typename TCrsMatrix::device_type,typename TCrsMatrix::non_const_scalar_type>::value> op ;
      const typename CrsMatrixType::ordinal_type nrow = A.numRows();
      op.m_A = A ;
      op.m_x = x ;
      op.m_y = y ;
      op.beta = betav;
      op.alpha = alphav;
      op.n = x.dimension(1);
      Kokkos::parallel_for(nrow*ThreadsPerRow<typename TCrsMatrix::device_type,typename TCrsMatrix::non_const_scalar_type>::value , op);

#else // NOT KOKKOS_FAST_COMPILE

      MV_MultiplyFunctor<RangeVectorType, CrsMatrixType, DomainVectorType, CoeffVector1Type, CoeffVector2Type, 2, 2,
	          ThreadsPerRow<typename TCrsMatrix::device_type, typename TCrsMatrix::non_const_scalar_type>::value> op ;

      int numVecs = x.dimension_1();
      CoeffVector1 beta = betav;
      CoeffVector2 alpha = alphav;

      if (doalpha != 2) {
	      alpha = CoeffVector2("CrsMatrix::auto_a", numVecs);
	      typename CoeffVector2::HostMirror h_a = Kokkos::create_mirror_view(alpha);
	      typename CoeffVector2::scalar_type s_a = (typename CoeffVector2::scalar_type) doalpha;

	      for (int i = 0; i < numVecs; ++i)
	        h_a(i) = s_a;

	      Kokkos::deep_copy(alpha, h_a);
      }

      if (dobeta != 2) {
	      beta = CoeffVector1("CrsMatrix::auto_b", numVecs);
	      typename CoeffVector1::HostMirror h_b = Kokkos::create_mirror_view(beta);
	      typename CoeffVector1::scalar_type s_b = (typename CoeffVector1::scalar_type) dobeta;

	      for(int i = 0; i < numVecs; i++)
	        h_b(i) = s_b;

	      Kokkos::deep_copy(beta, h_b);
      }

      const typename CrsMatrixType::ordinal_type nrow = A.numRows();
      op.m_A = A;
      op.m_x = x;
      op.m_y = y;
      op.beta = beta;
      op.alpha = alpha;
      op.n = x.dimension_1();
      Kokkos::parallel_for (nrow * ThreadsPerRow<typename TCrsMatrix::device_type,
                            typename TCrsMatrix::non_const_scalar_type>::value, op);
#endif // KOKKOS_FAST_COMPILE
    }
  }

  template<class RangeVector,
	   class TCrsMatrix,
	   class DomainVector,
	   class CoeffVector1,
	   class CoeffVector2,
	   int doalpha,
	   int dobeta>
  void
  MV_Multiply (typename Kokkos::Impl::enable_if<DomainVector::Rank == 1, const CoeffVector1>::type& betav,
	       const RangeVector &y,
	       const CoeffVector2 &alphav,
	       const TCrsMatrix& A,
	       const DomainVector& x)
  {
    if (doalpha == 0) {
      if (dobeta==2) {
	      V_MulScalar(y,betav,y);
      }
      else {
	      V_MulScalar(y,dobeta,y);
      }
      return;
    } else {
      typedef View< typename RangeVector::non_const_data_type ,
                    typename RangeVector::array_layout ,
                    typename RangeVector::device_type ,
                    typename RangeVector::memory_traits >
      RangeVectorType;

      typedef View< typename DomainVector::const_data_type ,
                    typename DomainVector::array_layout ,
                    typename DomainVector::device_type ,
                    Kokkos::MemoryRandomRead >
      DomainVectorType;

      typedef View< typename CoeffVector1::const_data_type ,
                    typename CoeffVector1::array_layout ,
                    typename CoeffVector1::device_type ,
                    Kokkos::MemoryRandomRead >
      CoeffVector1Type;

      typedef View< typename CoeffVector2::const_data_type ,
                    typename CoeffVector2::array_layout ,
                    typename CoeffVector2::device_type ,
                    Kokkos::MemoryRandomRead >
      CoeffVector2Type;

      typedef CrsMatrix<typename TCrsMatrix::const_scalar_type,
                        typename TCrsMatrix::ordinal_type,
                        typename TCrsMatrix::device_type,
                        typename TCrsMatrix::memory_traits,
                        typename TCrsMatrix::size_type>
      CrsMatrixType;

      Impl::MV_Multiply_Check_Compatibility(betav,y,alphav,A,x,doalpha,dobeta);

#ifndef KOKKOS_FAST_COMPILE
      MV_MultiplySingleFunctor<RangeVectorType, CrsMatrixType, DomainVectorType,
                               CoeffVector1Type, CoeffVector2Type, doalpha, dobeta
	,ThreadsPerRow<typename CrsMatrixType::device_type,typename CrsMatrixType::non_const_scalar_type>::value> op ;
      const typename CrsMatrixType::ordinal_type nrow = A.numRows();
      op.m_A = A ;
      op.m_x = x ;
      op.m_y = y ;
      op.beta = betav;
      op.alpha = alphav;
      op.n = x.dimension(1);
      Kokkos::parallel_for (nrow * ThreadsPerRow<typename TCrsMatrix::device_type,
                            typename TCrsMatrix::non_const_scalar_type>::value, op);

#else // NOT KOKKOS_FAST_COMPILE

      MV_MultiplySingleFunctor<RangeVectorType, CrsMatrixType, DomainVectorType,
                               CoeffVector1Type, CoeffVector2Type, 2, 2,
	ThreadsPerRow<typename CrsMatrixType::device_type, typename CrsMatrixType::non_const_scalar_type>::value> op;

      int numVecs = x.dimension_1();
      CoeffVector1 beta = betav;
      CoeffVector2 alpha = alphav;

      if(doalpha!=2) {
	      alpha = CoeffVector2("CrsMatrix::auto_a", numVecs);
	      typename CoeffVector2::HostMirror h_a = Kokkos::create_mirror_view(alpha);
	      typename CoeffVector2::scalar_type s_a = (typename CoeffVector2::scalar_type) doalpha;

	      for(int i = 0; i < numVecs; i++)
	        h_a(i) = s_a;

	      Kokkos::deep_copy(alpha, h_a);
      }
      if(dobeta!=2) {
	      beta = CoeffVector1("CrsMatrix::auto_b", numVecs);
	      typename CoeffVector1::HostMirror h_b = Kokkos::create_mirror_view(beta);
	      typename CoeffVector1::scalar_type s_b = (typename CoeffVector1::scalar_type) dobeta;

	      for(int i = 0; i < numVecs; i++)
	        h_b(i) = s_b;

	      Kokkos::deep_copy(beta, h_b);
      }
      const typename CrsMatrixType::ordinal_type nrow = A.numRows();
      op.m_A = A ;
      op.m_x = x ;
      op.m_y = y ;
      op.beta = beta;
      op.alpha = alpha;
      op.n = x.dimension_1();
      Kokkos::parallel_for (nrow * ThreadsPerRow<typename TCrsMatrix::device_type,
                            typename TCrsMatrix::non_const_scalar_type>::value, op);
#endif // KOKKOS_FAST_COMPILE
    }
  }

  template<class RangeVector, class CrsMatrix, class DomainVector, class CoeffVector1, class CoeffVector2>
  void
  MV_Multiply (const CoeffVector1& betav,
	       const RangeVector& y,
	       const CoeffVector2& alphav,
	       const CrsMatrix& A,
	       const DomainVector& x,
	       int beta,
	       int alpha)
  {
    if (beta == 0) {
      if(alpha == 0)
	MV_Multiply<RangeVector, CrsMatrix, DomainVector, CoeffVector1, CoeffVector2, 0, 0>(betav, y, alphav, A ,  x);
      else if(alpha == 1)
	MV_Multiply<RangeVector, CrsMatrix, DomainVector, CoeffVector1, CoeffVector2, 1, 0>(betav, y, alphav, A ,  x);
      else if(alpha == -1)
	MV_Multiply < RangeVector, CrsMatrix, DomainVector, CoeffVector1, CoeffVector2, -1, 0 > (betav, y, alphav, A ,  x);
      else
	MV_Multiply<RangeVector, CrsMatrix, DomainVector, CoeffVector1, CoeffVector2, 2, 0>(betav, y, alphav, A ,  x);
    } else if(beta == 1) {
      if(alpha == 0)
	return;
      else if(alpha == 1)
	MV_Multiply<RangeVector, CrsMatrix, DomainVector, CoeffVector1, CoeffVector2, 1, 1>(betav, y, alphav, A ,  x);
      else if(alpha == -1)
	MV_Multiply < RangeVector, CrsMatrix, DomainVector, CoeffVector1, CoeffVector2, -1, 1 > (betav, y, alphav, A ,  x);
      else
	MV_Multiply<RangeVector, CrsMatrix, DomainVector, CoeffVector1, CoeffVector2, 2, 1>(betav, y, alphav, A ,  x);
    } else if(beta == -1) {
      if(alpha == 0)
	MV_Multiply<RangeVector, CrsMatrix, DomainVector, CoeffVector1, CoeffVector2, 0, -1>(betav, y, alphav, A ,  x);
      else if(alpha == 1)
	MV_Multiply < RangeVector, CrsMatrix, DomainVector, CoeffVector1, CoeffVector2, 1, -1 > (betav, y, alphav, A ,  x);
      else if(alpha == -1)
	MV_Multiply < RangeVector, CrsMatrix, DomainVector, CoeffVector1, CoeffVector2, -1, -1 > (betav, y, alphav, A ,  x);
      else
	MV_Multiply < RangeVector, CrsMatrix, DomainVector, CoeffVector1, CoeffVector2, 2, -1 > (betav, y, alphav, A ,  x);
    } else {
      if(alpha == 0)
	MV_Multiply<RangeVector, CrsMatrix, DomainVector, CoeffVector1, CoeffVector2, 0, 2>(betav, y, alphav, A ,  x);
      else if(alpha == 1)
	MV_Multiply<RangeVector, CrsMatrix, DomainVector, CoeffVector1, CoeffVector2, 1, 2>(betav, y, alphav, A ,  x);
      else if(alpha == -1)
	MV_Multiply < RangeVector, CrsMatrix, DomainVector, CoeffVector1, CoeffVector2, -1, 2 > (betav, y, alphav, A ,  x);
      else
	MV_Multiply<RangeVector, CrsMatrix, DomainVector, CoeffVector1, CoeffVector2, 2, 2>(betav, y, alphav, A ,  x);
    }
  }

  template <class RangeVector, class CrsMatrix, class DomainVector,
	    class Value1, class Layout1, class Device1, class MemoryManagement1,
	    class Value2, class Layout2, class Device2, class MemoryManagement2>
  void
  MV_Multiply (const Kokkos::View<Value1, Layout1, Device1, MemoryManagement1>& betav,
	       const RangeVector& y,
	       const Kokkos::View<Value2, Layout2, Device2, MemoryManagement2>& alphav,
	       const CrsMatrix& A,
	       const DomainVector& x)
  {
    return MV_Multiply (betav, y, alphav, A, x, 2, 2);
  }

  template <class RangeVector, class CrsMatrix, class DomainVector,
	    class Value1, class Layout1, class Device1, class MemoryManagement1>
  void
  MV_Multiply (const RangeVector& y,
	       const Kokkos::View<Value1, Layout1, Device1, MemoryManagement1>& alphav,
	       const CrsMatrix& A,
	       const DomainVector& x)
  {
    return MV_Multiply (alphav, y, alphav, A, x, 0, 2);
  }

  template<class RangeVector, class CrsMatrix, class DomainVector>
  void
  MV_Multiply (const RangeVector& y,
	       const CrsMatrix& A,
	       const DomainVector& x)
  {
    // FIXME (mfh 21 Jun 2013) The way this code is supposed to work, is
    // that it tests at run time for each TPL in turn.  Shouldn't it
    // rather dispatch on the Device type?  But I suppose the "try"
    // functions do that.
    //
    // We want to condense this a bit: "Try TPLs" function that tests
    // all the suitable TPLs at run time.  This would be a run-time test
    // that compares the Scalar and Device types to those accepted by
    // the TPL(s).
#ifdef KOKKOS_USE_CUSPARSE
    if (MV_Multiply_Try_CuSparse (0.0, y, 1.0, A, x)) {
      return;
    }
#endif // KOKKOS_USE_CUSPARSE
#ifdef KOKKOS_USE_MKL
    if (MV_Multiply_Try_MKL (0.0, y, 1.0, A, x)) {
      return;
    }
#endif // KOKKOS_USE_MKL
    typedef Kokkos::View<typename DomainVector::scalar_type*, typename DomainVector::device_type> aVector;
    aVector a;

    return MV_Multiply (a, y, a, A, x, 0, 1);
  }

  template<class RangeVector, class CrsMatrix, class DomainVector>
  void
  MV_Multiply (const RangeVector& y,
	       typename DomainVector::const_scalar_type s_a,
	       const CrsMatrix& A,
	       const DomainVector& x)
  {
#ifdef KOKKOS_USE_CUSPARSE
    if (MV_Multiply_Try_CuSparse (0.0, y, s_a, A, x)) {
      return;
    }
#endif // KOKKOS_USE_CUSPARSE
#ifdef KOKKOS_USE_MKL
    if (MV_Multiply_Try_MKL (0.0, y, s_a, A, x)) {
      return;
    }
#endif // KOKKOS_USE_MKL
    typedef Kokkos::View<typename RangeVector::scalar_type*, typename RangeVector::device_type> aVector;
    aVector a;
    const int numVecs = x.dimension_1();

    if (s_a == -1) {
      return MV_Multiply (a, y, a, A, x, 0, -1);
    } else if (s_a == 1) {
      return MV_Multiply (a, y, a, A, x, 0, 1);
    }

    if (s_a != 0) {
      a = aVector("a", numVecs);
      typename aVector::HostMirror h_a = Kokkos::create_mirror_view (a);
      for (int i = 0; i < numVecs; ++i) {
	h_a(i) = s_a;
      }
      Kokkos::deep_copy(a, h_a);
      return MV_Multiply (a, y, a, A, x, 0, 2);
    }
  }

  template<class RangeVector, class CrsMatrix, class DomainVector>
  void
  MV_Multiply (typename RangeVector::const_scalar_type s_b,
	       const RangeVector& y,
	       typename DomainVector::const_scalar_type s_a,
	       const CrsMatrix& A,
	       const DomainVector& x)
  {
#ifdef KOKKOS_USE_CUSPARSE
    if (MV_Multiply_Try_CuSparse (s_b, y, s_a, A, x)) {
      return;
    }
#endif // KOKKOSE_USE_CUSPARSE
#ifdef KOKKOS_USE_MKL
    if (MV_Multiply_Try_MKL (s_b, y, s_a, A, x)) {
      return;
    }
#endif // KOKKOS_USE_MKL
    typedef Kokkos::View<typename RangeVector::scalar_type*, typename RangeVector::device_type> aVector;
    aVector a;
    aVector b;
    int numVecs = x.dimension_1();

    if(numVecs==1)
    if (s_b == 0) {
      if (s_a == 0)
        return MV_Multiply (a, y, a, A, x, 0, 0);
      else if (s_a == 1)
        return MV_Multiply (a, y, a, A, x, 0, 1);
      else if (s_a == -1)
        return MV_Multiply (a, y, a, A, x, 0, -1);
      else {
        a = aVector("a", numVecs);
        typename aVector::HostMirror h_a = Kokkos::create_mirror_view(a);
        for (int i = 0; i < numVecs; ++i) {
          h_a(i) = s_a;
        }
        Kokkos::deep_copy (a, h_a);
        return MV_Multiply (a, y, a, A, x, 0, 2);
      }
    } else if (s_b == 1) {
      if (s_a == 0)
        return MV_Multiply (a, y, a, A, x, 1, 0);
      else if (s_a == 1)
        return MV_Multiply (a, y, a, A, x, 1, 1);
      else if (s_a == -1)
        return MV_Multiply (a, y, a, A, x, 1, -1);
      else {
        a = aVector("a", numVecs);
        typename aVector::HostMirror h_a = Kokkos::create_mirror_view(a);
        for (int i = 0; i < numVecs; ++i) {
          h_a(i) = s_a;
        }
        Kokkos::deep_copy (a, h_a);
        return MV_Multiply (a, y, a, A, x, 1, 2);
      }
    } else if (s_b == -1) {
      if (s_a == 0)
        return MV_Multiply (a, y, a, A, x, -1, 0);
      else if (s_a == 1)
        return MV_Multiply (a, y, a, A, x, -1, 1);
      else if (s_a == -1)
        return MV_Multiply (a, y, a, A, x, -1, -1);
      else {
        a = aVector("a", numVecs);
        typename aVector::HostMirror h_a = Kokkos::create_mirror_view(a);
        for (int i = 0; i < numVecs; ++i) {
          h_a(i) = s_a;
        }
        Kokkos::deep_copy (a, h_a);
        return MV_Multiply (a, y, a, A, x, -1, 2);
      }
    } else {
      b = aVector("b", numVecs);
      typename aVector::HostMirror h_b = Kokkos::create_mirror_view(b);
      for (int i = 0; i < numVecs; ++i) {
        h_b(i) = s_b;
      }
      Kokkos::deep_copy(b, h_b);

      if (s_a == 0)
        return MV_Multiply (b, y, a, A, x, 2, 0);
      else if (s_a == 1)
        return MV_Multiply (b, y, a, A, x, 2, 1);
      else if (s_a == -1)
        return MV_Multiply (b, y, a, A, x, 2, -1);
      else {
        a = aVector("a", numVecs);
        typename aVector::HostMirror h_a = Kokkos::create_mirror_view(a);
        for (int i = 0; i < numVecs; ++i) {
          h_a(i) = s_a;
        }
        Kokkos::deep_copy (a, h_a);
        return MV_Multiply (b, y, a, A, x, 2, 2);
      }
    }
  }

  namespace KokkosCrsMatrix {
    /// \brief Copy the CrsMatrix B into the CrsMatrix A.
    /// \tparam CrsMatrixDst CrsMatrix specialization of the destination
    ///   (Dst) matrix A.
    /// \tparam CrsMatrixSrc CrsMatrix specialization of the source
    ///   (Src) matrix B.
    ///
    /// The two CrsMatrix specializations CrsMatrixDst and CrsMatrixSrc
    /// need not be the same.  However, it must be possible to deep_copy
    /// their column indices and their values.
    ///
    /// The target matrix must already be allocated, and must have the
    /// same number of rows and number of entries as the source matrix.
    /// It need not have the same row map as the source matrix.
    template <class CrsMatrixDst, class CrsMatrixSrc>
    void deep_copy (CrsMatrixDst A, CrsMatrixSrc B) {
      Kokkos::deep_copy(A.graph.entries, B.graph.entries);
      // FIXME (mfh 09 Aug 2013) This _should_ copy the row map.  We
      // couldn't do it before because the row map was const, forbidding
      // deep_copy.
      //
      //Kokkos::deep_copy(A.graph.row_map,B.graph.row_map);
      Kokkos::deep_copy(A.values, B.values);

      // FIXME (mfh 09 Aug 2013) Be sure to copy numRows, numCols, and
      // nnz as well.
      // (CRT 25 Sep 2013) don't copy rather check that they match.
      // Deep_copy in Kokkos is intended for copy between compatible objects.
    }
  } // namespace KokkosCrsMatrix

} // namespace Kokkos

#endif /* KOKKOS_CRSMATRIX_H_ */
