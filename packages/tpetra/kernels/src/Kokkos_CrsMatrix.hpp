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

#include <Kokkos_Core.hpp>
#include <Kokkos_StaticCrsGraph.hpp>
#include <Kokkos_MV.hpp>

#ifdef KOKKOS_USE_CUSPARSE
#  include <cusparse.h>
#  include <Kokkos_CrsMatrix_CuSparse.hpp>
#endif // KOKKOS_USE_CUSPARSE

#ifdef KOKKOS_USE_MKL
#  include <mkl.h>
#  include <mkl_spblas.h>
#  include <Kokkos_CrsMatrix_MKL.hpp>
#endif // KOKKOS_USE_MKL

//#include <Kokkos_Vectorization.hpp>
#include <impl/Kokkos_Error.hpp>

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
/// typedef typename SparseRowView<MatrixType>::value_type value_type;
/// typedef typename SparseRowView<MatrixType>::ordinal_type ordinal_type;
///
/// SparseRowView<MatrixType> A_i = ...;
/// const int numEntries = A_i.length;
/// for (int k = 0; k < numEntries; ++k) {
///   value_type A_ij = A_i.value (k);
///   ordinal_type j = A_i.colidx (k);
///   // ... do something with A_ij and j ...
/// }
/// \endcode
///
/// MatrixType must provide the \c value_type and \c ordinal_type
/// typedefs.  In addition, it must make sense to use SparseRowView to
/// view a row of MatrixType.  In particular, the values and column
/// indices of a row must be accessible using the <tt>values</tt>
/// resp. <tt>colidx</tt> arrays given to the constructor of this
/// class, with a constant <tt>stride</tt> between successive entries.
/// The stride is one for the compressed sparse row storage format (as
/// is used by CrsMatrix), but may be greater than one for other
/// sparse matrix storage formats (e.g., ELLPACK or jagged diagonal).
template<class MatrixType, class SizeType = typename MatrixType::size_type>
struct SparseRowView {
  //! The type of the values in the row.
  typedef typename MatrixType::value_type value_type;
  //! The type of the column indices in the row.
  typedef typename MatrixType::ordinal_type ordinal_type;
  //! The type of array offsets and strides.
  typedef SizeType size_type;

private:
  //! Array of values in the row.
  value_type* values_;
  //! Array of (local) column indices in the row.
  ordinal_type* colidx_;
  /// \brief Stride between successive entries in the row.
  ///
  /// For compressed sparse row (CSR) storage, this is always one.
  /// This might be greater than one for storage formats like ELLPACK.
  const size_type stride_;

public:
  /// \brief Constructor
  ///
  /// \param values [in] Array of the row's values.
  /// \param colidx [in] Array of the row's column indices.
  /// \param stride [in] (Constant) stride between matrix entries in
  ///   each of the above arrays.
  /// \param count [in] Number of entries in the row.
  KOKKOS_INLINE_FUNCTION
  SparseRowView (value_type* const values,
                 ordinal_type* const colidx__,
                 const size_type& stride,
                 const size_type& count) :
    values_ (values), colidx_ (colidx__), stride_ (stride), length (count)
  {}

  /// \brief Constructor
  ///
  /// \param values [in] Array of the row's values.
  /// \param colidx [in] Array of the row's column indices.
  /// \param stride [in] (Constant) stride between matrix entries in
  ///   each of the above arrays.
  /// \param count [in] Number of entries in the row.
  KOKKOS_INLINE_FUNCTION
  SparseRowView (const typename MatrixType::values_type& values,
                 const typename MatrixType::index_type& colidx__,
                 const size_type& stride,
                 const size_type& count,
                 const size_type& idx) :
    values_ (&values(idx)), colidx_ (&colidx__(idx)), stride_ (stride), length (count)
  {}

  /// \brief Number of entries in the row.
  ///
  /// This is a public const field rather than a public const method,
  /// in order to avoid possible overhead of a method call if the
  /// compiler is unable to inline that method call.
  ///
  /// We assume that rows contain no duplicate entries (i.e., entries
  /// with the same column index).  Thus, a row may have up to
  /// A.numCols() entries.  This means that the correct type of
  /// 'length' is ordinal_type.
  const ordinal_type length;

  /// \brief Reference to the value of entry i in this row of the sparse matrix.
  ///
  /// "Entry i" is not necessarily the entry with column index i, nor
  /// does i necessarily correspond to the (local) row index.
  KOKKOS_INLINE_FUNCTION
  value_type& value (const ordinal_type& i) const {
    return values_[i*stride_];
  }

  /// \brief Reference to the column index of entry i in this row of the sparse matrix.
  ///
  /// "Entry i" is not necessarily the entry with column index i, nor
  /// does i necessarily correspond to the (local) row index.
  KOKKOS_INLINE_FUNCTION
  ordinal_type& colidx (const ordinal_type& i) const {
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
template<class MatrixType, class SizeType = typename MatrixType::size_type>
struct SparseRowViewConst {
  //! The type of the values in the row.
  typedef const typename MatrixType::non_const_value_type value_type;
  //! The type of the column indices in the row.
  typedef const typename MatrixType::non_const_ordinal_type ordinal_type;
  //! The type of array offsets and strides.
  typedef SizeType size_type;

private:
  //! Array of values in the row.
  value_type* values_;
  //! Array of (local) column indices in the row.
  ordinal_type* colidx_;
  /// \brief Stride between successive entries in the row.
  ///
  /// For compressed sparse row (CSR) storage, this is always one.
  /// This might be greater than one for storage formats like ELLPACK.
  const size_type stride_;

public:
  /// \brief Constructor
  ///
  /// \param values [in] Array of the row's values.
  /// \param colidx [in] Array of the row's column indices.
  /// \param stride [in] (Constant) stride between matrix entries in
  ///   each of the above arrays.
  /// \param count [in] Number of entries in the row.
  KOKKOS_INLINE_FUNCTION
  SparseRowViewConst (value_type* const values,
                      ordinal_type* const colidx__,
                      const size_type& stride,
                      const size_type& count) :
    values_ (values), colidx_ (colidx__), stride_ (stride), length (count)
  {}

  /// \brief Constructor
  ///
  /// \param values [in] Array of the row's values.
  /// \param colidx [in] Array of the row's column indices.
  /// \param stride [in] (Constant) stride between matrix entries in
  ///   each of the above arrays.
  /// \param count [in] Number of entries in the row.
  KOKKOS_INLINE_FUNCTION
  SparseRowViewConst (const typename MatrixType::values_type& values,
                      const typename MatrixType::index_type& colidx__,
                      const size_type& stride,
                      const size_type& count,
                      const size_type& idx) :
    values_ (&values(idx)), colidx_ (&colidx__(idx)), stride_ (stride), length (count)
  {}

  /// \brief Number of entries in the row.
  ///
  /// This is a public const field rather than a public const method,
  /// in order to avoid possible overhead of a method call if the
  /// compiler is unable to inline that method call.
  ///
  /// We assume that rows contain no duplicate entries (i.e., entries
  /// with the same column index).  Thus, a row may have up to
  /// A.numCols() entries.  This means that the correct type of
  /// 'length' is ordinal_type.
  const ordinal_type length;

  /// \brief (Const) reference to the value of entry i in this row of
  ///   the sparse matrix.
  ///
  /// "Entry i" is not necessarily the entry with column index i, nor
  /// does i necessarily correspond to the (local) row index.
  KOKKOS_INLINE_FUNCTION
  value_type& value (const ordinal_type& i) const {
    return values_[i*stride_];
  }

  /// \brief (Const) reference to the column index of entry i in this
  ///   row of the sparse matrix.
  ///
  /// "Entry i" is not necessarily the entry with column index i, nor
  /// does i necessarily correspond to the (local) row index.
  KOKKOS_INLINE_FUNCTION
  ordinal_type& colidx (const ordinal_type& i) const {
    return colidx_[i*stride_];
  }
};

// A simple struct for storing a kernel launch configuration.
// This is currently used by CrsMatrix to allow the user to have some control
// over how kernels are launched, however it is currently only exercised by
// Stokhos.  This is a simpler case of "state" needed by TPLs, and at this point
// is just a hack until we figure out how to support state in a general,
// extensible way.
struct DeviceConfig {
  struct Dim3 {
    size_t x, y, z;
    Dim3(const size_t x_, const size_t y_ = 1, const size_t z_ = 1) :
      x(x_), y(y_), z(z_) {}
  };

  Dim3 block_dim;
  size_t num_blocks;
  size_t num_threads_per_block;

  DeviceConfig(const size_t num_blocks_ = 0,
               const size_t threads_per_block_x_ = 0,
               const size_t threads_per_block_y_ = 0,
               const size_t threads_per_block_z_ = 1) :
    block_dim(threads_per_block_x_,threads_per_block_y_,threads_per_block_z_),
    num_blocks(num_blocks_),
    num_threads_per_block(block_dim.x * block_dim.y * block_dim.z)
    {}
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
template<class ScalarType,
         class OrdinalType,
         class Device,
         class MemoryTraits = void,
         class SizeType = typename Kokkos::ViewTraits<OrdinalType*, Device, void, void>::size_type>
class CrsMatrix {
private:
  typedef typename Kokkos::ViewTraits<ScalarType*,Device,void,void>::host_mirror_space host_mirror_space ;
public:
  //! Type of the matrix's execution space.
  typedef typename Device::execution_space execution_space;
  //! Type of each value in the matrix.
  typedef ScalarType value_type;
  //! Type of each (column) index in the matrix.
  typedef OrdinalType ordinal_type;
  typedef MemoryTraits memory_traits;
  /// \brief Type of each entry of the "row map."
  ///
  /// The "row map" corresponds to the \c ptr array of row offsets in
  /// compressed sparse row (CSR) storage.
  typedef SizeType size_type;

  //! Type of a host-memory mirror of the sparse matrix.
  typedef CrsMatrix<ScalarType, OrdinalType, host_mirror_space, MemoryTraits> HostMirror;
  //! Type of the graph structure of the sparse matrix.
  typedef Kokkos::StaticCrsGraph<OrdinalType, Kokkos::LayoutLeft, Device, SizeType> StaticCrsGraphType;
  //! Type of column indices in the sparse matrix.
  typedef typename StaticCrsGraphType::entries_type index_type;
  //! Nonconst version of the type of column indices in the sparse matrix.
  typedef typename index_type::non_const_value_type non_const_ordinal_type;
  //! Type of the "row map" (which contains the offset for each row's data).
  typedef typename StaticCrsGraphType::row_map_type row_map_type;
  //! Kokkos Array type of the entries (values) in the sparse matrix.
  typedef Kokkos::View<value_type*, Kokkos::LayoutRight, execution_space, MemoryTraits> values_type;
  //! Const version of the type of the entries in the sparse matrix.
  typedef typename values_type::const_value_type const_value_type;
  //! Nonconst version of the type of the entries in the sparse matrix.
  typedef typename values_type::non_const_value_type non_const_value_type;

#ifdef KOKKOS_USE_CUSPARSE
  cusparseHandle_t cusparse_handle;
  cusparseMatDescr_t cusparse_descr;
#endif // KOKKOS_USE_CUSPARSE

  /// \name Storage of the actual sparsity structure and values.
  ///
  /// CrsMatrix uses the compressed sparse row (CSR) storage format to
  /// store the sparse matrix.  CSR is also called "compressed row
  /// storage"; hence the name, which it inherits from Tpetra and from
  /// Epetra before it.
  //@{
  //! The graph (sparsity structure) of the sparse matrix.
  StaticCrsGraphType graph;
  //! The 1-D array of values of the sparse matrix.
  values_type values;
  //@}

  /// \brief Launch configuration that can be used by
  ///   overloads/specializations of MV_multiply().
  ///
  /// This is a hack and needs to be replaced by a general
  /// state mechanism.
  DeviceConfig dev_config;

  /// \brief Default constructor; constructs an empty sparse matrix.
  ///
  /// FIXME (mfh 09 Aug 2013) numCols and nnz should be properties of
  /// the graph, not the matrix.  Then CrsMatrix needs methods to get
  /// these from the graph.
  CrsMatrix () :
    numCols_ (0)
  {}

  //! Copy constructor (shallow copy).
  template<typename SType,
           typename OType,
           class DType,
           class MTType,
           typename IType>
  CrsMatrix (const CrsMatrix<SType,OType,DType,MTType,IType> & B) :
    graph (B.graph),
    values (B.values),
    numCols_ (B.numCols ())
  {}

  /// \brief Construct with a graph that will be shared.
  ///
  /// Allocate the values array for subsquent fill.
  CrsMatrix (const std::string& arg_label,
             const StaticCrsGraphType& arg_graph) :
    graph (arg_graph),
    values (arg_label, arg_graph.entries.dimension_0 ()),
    numCols_ (maximum_entry (arg_graph) + 1)
  {}

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
             size_type annz,
             ScalarType* val,
             OrdinalType* rows,
             OrdinalType* cols,
             bool pad = false)
  {
    (void) pad;
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

  /// \brief Constructor that accepts a row map, column indices, and
  ///   values.
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
  CrsMatrix (const std::string& label,
             const OrdinalType nrows,
             const OrdinalType ncols,
             const size_type annz,
             const values_type& vals,
             const row_map_type& rows,
             const index_type& cols) :
    graph (cols, rows),
    values (vals),
    numCols_ (ncols)
  {
    const ordinal_type actualNumRows = (rows.dimension_0 () != 0) ?
      static_cast<ordinal_type> (rows.dimension_0 () - static_cast<size_type> (1)) :
      static_cast<ordinal_type> (0);
    if (nrows != actualNumRows) {
      std::ostringstream os;
      os << "Input argument nrows = " << nrows << " != the actual number of "
        "rows " << actualNumRows << " according to the 'rows' input argument.";
      throw std::invalid_argument (os.str ());
    }
    if (annz != nnz ()) {
      std::ostringstream os;
      os << "Input argument annz = " << annz
         << " != this->nnz () = " << nnz () << ".";
      throw std::invalid_argument (os.str ());
    }

#ifdef KOKKOS_USE_CUSPARSE
    cusparseCreate (&cusparse_handle);
    cusparseCreateMatDescr (&cusparse_descr);
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
  CrsMatrix (const std::string& label,
             const OrdinalType& ncols,
             const values_type& vals,
             const StaticCrsGraphType& graph_) :
    graph (graph_),
    values (vals),
    numCols_ (ncols)
  {
#ifdef KOKKOS_USE_CUSPARSE
    cusparseCreate (&cusparse_handle);
    cusparseCreateMatDescr (&cusparse_descr);
#endif // KOKKOS_USE_CUSPARSE
  }

  void
  import (const std::string &label,
          const OrdinalType nrows,
          const OrdinalType ncols,
          const size_type annz,
          ScalarType* val,
          OrdinalType* rows,
          OrdinalType* cols);

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
  sumIntoValues (const OrdinalType rowi,
                 const OrdinalType cols[],
                 const OrdinalType ncol,
                 ScalarType vals[],
                 const bool force_atomic = false) const
  {
    SparseRowView<CrsMatrix> row_view = this->row (rowi);
    const size_type length = row_view.length;
    for (OrdinalType i = 0; i < ncol; ++i) {
      for (size_type j = 0; j < length; ++j) {
        if (row_view.colidx(j) == cols[i]) {
          if (force_atomic) {
            atomic_add(&row_view.value(j), vals[i]);
          } else {
            row_view.value(j) += vals[i];
          }
          break;
        }
      }
    }
  }

  // FIXME (mfh 29 Sep 2013) See above notes on sumIntoValues.
  KOKKOS_INLINE_FUNCTION
  void
  replaceValues (const OrdinalType rowi,
                 const OrdinalType cols[],
                 const OrdinalType ncol,
                 ScalarType vals[],
                 const bool force_atomic = false) const
  {
    SparseRowView<CrsMatrix> row_view = this->row (rowi);
    const int length = row_view.length;
    for (OrdinalType i = 0; i < ncol; ++i) {
      for (int j = 0; j < length; ++j) {
        if (row_view.colidx(j) == cols[i]) {
          if (force_atomic) {
            atomic_assign(&row_view.value(j), vals[i]);
          } else {
            row_view.value(j) = vals[i];
          }
        }
      }
    }
  }

  //! Attempt to assign the input matrix to \c *this.
  template<typename aScalarType, typename aOrdinalType, class aDevice, class aMemoryTraits,typename aSizeType>
  CrsMatrix&
  operator= (const CrsMatrix<aScalarType, aOrdinalType, aDevice, aMemoryTraits, aSizeType>& mtx)
  {
    numCols_ = mtx.numCols ();
    graph = mtx.graph;
    values = mtx.values;
    dev_config = mtx.dev_config;
    return *this;
  }

  //! The number of rows in the sparse matrix.
  KOKKOS_INLINE_FUNCTION ordinal_type numRows () const {
    return graph.numRows ();
  }

  //! The number of columns in the sparse matrix.
  KOKKOS_INLINE_FUNCTION ordinal_type numCols () const {
    return numCols_;
  }

  //! The number of stored entries in the sparse matrix.
  KOKKOS_INLINE_FUNCTION size_type nnz () const {
    return graph.entries.dimension_0 ();
  }

  friend struct SparseRowView<CrsMatrix>;

  /// \brief Return a view of row i of the matrix.
  ///
  /// If row i does not belong to the matrix, return an empty view.
  KOKKOS_INLINE_FUNCTION
  SparseRowView<CrsMatrix> row (const ordinal_type i) const {
    const size_type start = graph.row_map(i);
    const size_type count = graph.row_map(i+1) - start;

    if (count == 0) {
      return SparseRowView<CrsMatrix> (NULL, NULL, 1, 0);
    } else {
      return SparseRowView<CrsMatrix> (values, graph.entries, 1, count, start);
    }
  }

  /// \brief Return a const view of row i of the matrix.
  ///
  /// If row i does not belong to the matrix, return an empty view.
  KOKKOS_INLINE_FUNCTION
  SparseRowViewConst<CrsMatrix> rowConst (const ordinal_type i) const {
    const size_type start = graph.row_map(i);
    const size_type count = graph.row_map(i+1) - start;

    if (count == 0) {
      return SparseRowViewConst<CrsMatrix> (NULL, NULL, 1, 0);
    } else {
      return SparseRowViewConst<CrsMatrix> (values, graph.entries, 1, count, start);
    }
  }

private:
  ordinal_type numCols_;
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< typename ScalarType , typename OrdinalType, class Device, class MemoryTraits, typename SizeType >
void
CrsMatrix<ScalarType , OrdinalType, Device, MemoryTraits, SizeType >::
import (const std::string &label,
        const OrdinalType nrows,
        const OrdinalType ncols,
        const size_type annz,
        ScalarType* val,
        OrdinalType* rows,
        OrdinalType* cols)
{
  std::string str = label;
  values = values_type (str.append (".values"), annz);

  numCols_ = ncols;

  // FIXME (09 Aug 2013) CrsArray only takes std::vector for now.
  // We'll need to fix that.
  std::vector<int> row_lengths (nrows, 0);

  // FIXME (mfh 21 Jun 2013) This calls for a parallel_for kernel.
  for (OrdinalType i = 0; i < nrows; ++i) {
    row_lengths[i] = rows[i + 1] - rows[i];
  }

  str = label;
  graph = Kokkos::create_staticcrsgraph<StaticCrsGraphType> (str.append (".graph"), row_lengths);
  typename values_type::HostMirror h_values = Kokkos::create_mirror_view (values);
  typename index_type::HostMirror h_entries = Kokkos::create_mirror_view (graph.entries);

  // FIXME (mfh 21 Jun 2013) This needs to be a parallel copy.
  // Furthermore, why are the arrays copied twice? -- once here, to a
  // host view, and once below, in the deep copy?
  for (size_type i = 0; i < annz; ++i) {
    if (val) {
      h_values(i) = val[i];
    }
    h_entries(i) = cols[i];
  }

  Kokkos::deep_copy (values, h_values);
  Kokkos::deep_copy (graph.entries, h_entries);
}

// FIXME (mfh 20 Mar 2015) Shouldn't this use ordinal_type instead of int?
template<class DeviceType>
inline int RowsPerThread(const int NNZPerRow) {
  if(NNZPerRow == 0) return 1;
  int result = 2;
  while(result*NNZPerRow <= 2048) {
    result*=2;
  }
  return result/2;
}
#ifdef KOKKOS_HAVE_CUDA
template<>
inline int RowsPerThread<Kokkos::Cuda>(const int NNZPerRow) {
  return 1;
}
#endif

//----------------------------------------------------------------------------

template<class DeviceType, typename ScalarType, int NNZPerRow=27>
struct MV_MultiplyShflThreadsPerRow {
private:
  typedef typename Kokkos::Impl::remove_const< ScalarType >::type value_type;

  // The shuffle operation only works with CUDA, and only works for
  // certain ScalarType types.
#ifdef KOKKOS_HAVE_CUDA
  enum { shfl_possible =
    Kokkos::Impl::is_same< DeviceType , Kokkos::Cuda >::value &&
    (
      Kokkos::Impl::is_same< value_type , unsigned int >::value ||
      Kokkos::Impl::is_same< value_type , int >::value ||
      Kokkos::Impl::is_same< value_type , float >::value ||
      Kokkos::Impl::is_same< value_type , double >::value
    )};
#else // NOT KOKKOS_HAVE_CUDA
  enum { shfl_possible = 0 };
#endif // KOKKOS_HAVE_CUDA

public:

#if defined( __CUDA_ARCH__ )
  enum { device_value = shfl_possible && ( 300 <= __CUDA_ARCH__ ) ?
         (NNZPerRow<8?2:
         (NNZPerRow<16?4:
         (NNZPerRow<32?8:
         (NNZPerRow<64?16:
         32))))
         :1 };
#else
  enum { device_value = 1 };
#endif

#ifdef KOKKOS_HAVE_CUDA
  inline static int host_value()
    { return shfl_possible && ( 300 <= Kokkos::Cuda::device_arch() ) ?
         (NNZPerRow<8?2:
         (NNZPerRow<16?4:
         (NNZPerRow<32?8:
         (NNZPerRow<64?16:
         32))))
         :1; }
#else // NOT KOKKOS_HAVE_CUDA
  inline static int host_value() { return 1; }
#endif // KOKKOS_HAVE_CUDA
};

//----------------------------------------------------------------------------

template<class RangeVector,
         class CrsMatrix,
         class DomainVector,
         class CoeffVector1,
         class CoeffVector2,
         int doalpha,
         int dobeta>
struct MV_MultiplyFunctor {
  typedef typename CrsMatrix::execution_space execution_space;
  typedef typename CrsMatrix::size_type       size_type;
  typedef typename CrsMatrix::ordinal_type    ordinal_type;
  typedef typename CrsMatrix::non_const_value_type value_type;
  typedef typename Kokkos::View<value_type*, execution_space> range_values;
  typedef typename Kokkos::TeamPolicy<execution_space>        team_policy;
  typedef typename team_policy::member_type                   team_member;

  CoeffVector1 beta;
  CoeffVector2 alpha;
  CrsMatrix  m_A;
  DomainVector  m_x;
  RangeVector  m_y;
  /// \brief The number of columns in the input and output MultiVectors.
  ///
  /// Its approxpriate type is therefore size_type, but we don't
  /// expect the input and output MultiVectors to have more columns
  /// than the sparse matrix has rows or columns.  Thus, we prefer the
  /// (likely both smaller and signed, vs. the larger and likely
  /// unsigned) size_type.
  ordinal_type n;
  int rows_per_thread;

  MV_MultiplyFunctor (const CoeffVector1 beta_,
                      const CoeffVector2 alpha_,
                      const CrsMatrix m_A_,
                      const DomainVector m_x_,
                      const RangeVector m_y_,
                      const ordinal_type n_,
                      const int rows_per_thread_) :
      beta (beta_), alpha (alpha_),
      m_A (m_A_), m_x (m_x_), m_y (m_y_), n (n_),
      rows_per_thread (rows_per_thread_)
  {}

  template<int UNROLL>
  KOKKOS_INLINE_FUNCTION void
  strip_mine (const team_member& dev, const size_type& iRow, const size_type& kk) const
  {
    value_type sum[UNROLL];

#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
    for (int k = 0; k < UNROLL; ++k) {
      // NOTE (mfh 09 Aug 2013) This requires that assignment from int
      // (in this case, 0) to value_type be defined.  It's not for
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

    const SparseRowView<CrsMatrix> row = m_A.row(iRow);

    // NOTE (mfh 20 Mar 2015) Unfortunately, Kokkos::Vectorization
    // lacks a typedef for determining the type of the return value of
    // begin().  I know that it returns int now, but this may change
    // at some point.
    //
    // The correct type of iEntry is ordinal_type.  This is because we
    // assume that rows have no duplicate entries.  As a result, a row
    // cannot have more entries than the number of columns in the
    // matrix.

#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
#ifdef KOKKOS_HAVE_PRAGMA_LOOPCOUNT
#pragma loop count (15)
#endif
#ifdef __CUDA_ARCH__
        for (ordinal_type iEntry = static_cast<ordinal_type> (threadIdx.x);
             iEntry < static_cast<ordinal_type> (row.length);
             iEntry += static_cast<ordinal_type> (blockDim.x)) {
#else
        for (ordinal_type iEntry = 0;
             iEntry < static_cast<ordinal_type> (row.length);
             iEntry ++) {
#endif
      const value_type val = row.value(iEntry);
      const ordinal_type ind = row.colidx(iEntry);

#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        sum[k] +=  val * m_x(ind, kk + k);
      }
    }

    if (doalpha == -1) {
      for (int ii=0; ii < UNROLL; ++ii) {
        value_type sumt=sum[ii];
        #ifdef __CUDA_ARCH__
        if (blockDim.x > 1)
          sumt += shfl_down(sumt, 1,blockDim.x);
        if (blockDim.x > 2)
          sumt += shfl_down(sumt, 2,blockDim.x);
        if (blockDim.x > 4)
          sumt += shfl_down(sumt, 4,blockDim.x);
        if (blockDim.x > 8)
          sumt += shfl_down(sumt, 8,blockDim.x);
        if (blockDim.x > 16)
          sumt += shfl_down(sumt, 16,blockDim.x);
        #endif
        sum[ii] = - sumt;
      }
    }
    else {
      for (int ii=0; ii < UNROLL; ++ii) {
        value_type sumt = sum[ii];
        #ifdef __CUDA_ARCH__
        if (blockDim.x > 1)
          sumt += shfl_down(sumt, 1,blockDim.x);
        if (blockDim.x > 2)
          sumt += shfl_down(sumt, 2,blockDim.x);
        if (blockDim.x > 4)
          sumt += shfl_down(sumt, 4,blockDim.x);
        if (blockDim.x > 8)
          sumt += shfl_down(sumt, 8,blockDim.x);
        if (blockDim.x > 16)
          sumt += shfl_down(sumt, 16,blockDim.x);
        #endif
        sum[ii] = sumt;
      }
    }

#ifdef __CUDA_ARCH__
    if (threadIdx.x==0) {
#else
    if (true) {
#endif
      if (doalpha * doalpha != 1) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (int k = 0; k < UNROLL; ++k) {
          sum[k] *= alpha(kk + k);
        }
      }

      if (dobeta == 0) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (int k = 0; k < UNROLL; ++k) {
          m_y(iRow, kk + k) = sum[k];
        }
      } else if (dobeta == 1) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (int k = 0; k < UNROLL; ++k) {
          m_y(iRow, kk + k) += sum[k];
        }
      } else if (dobeta == -1) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (int k = 0; k < UNROLL; ++k) {
          m_y(iRow, kk + k) = -m_y(iRow, kk + k) +  sum[k];
        }
      } else {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (int k = 0; k < UNROLL; ++k) {
          m_y(iRow, kk + k) = beta(kk + k) * m_y(iRow, kk + k) + sum[k] ;
        }
      }
    }
  }

  KOKKOS_INLINE_FUNCTION void
  strip_mine_1 (const team_member& dev, const size_type& iRow) const
  {
    value_type sum = 0;

    const SparseRowViewConst<CrsMatrix> row = m_A.rowConst(iRow);

    // NOTE (mfh 20 Mar 2015) Unfortunately, Kokkos::Vectorization
    // lacks a typedef for determining the type of the return value of
    // begin().  I know that it returns int now, but this may change
    // at some point.
    //
    // The correct type of iEntry is ordinal_type.  This is because we
    // assume that rows have no duplicate entries.  As a result, a row
    // cannot have more entries than the number of columns in the
    // matrix.

#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
#ifdef KOKKOS_HAVE_PRAGMA_LOOPCOUNT
#pragma loop count (15)
#endif
#ifdef __CUDA_ARCH__
    for (ordinal_type iEntry = static_cast<ordinal_type> (threadIdx.x);
       iEntry < static_cast<ordinal_type> (row.length);
       iEntry += static_cast<ordinal_type> (blockDim.x)) {
#else
    for (ordinal_type iEntry = 0;
       iEntry < static_cast<ordinal_type> (row.length);
       iEntry ++) {
#endif
      sum += row.value(iEntry) * m_x(row.colidx(iEntry),0);
    }
    #ifdef __CUDA_ARCH__
    if (blockDim.x > 1)
      sum += shfl_down(sum, 1,blockDim.x);
    if (blockDim.x > 2)
      sum += shfl_down(sum, 2,blockDim.x);
    if (blockDim.x > 4)
      sum += shfl_down(sum, 4,blockDim.x);
    if (blockDim.x > 8)
      sum += shfl_down(sum, 8,blockDim.x);
    if (blockDim.x > 16)
      sum += shfl_down(sum, 16,blockDim.x);
    #endif

#ifdef __CUDA_ARCH__
    if (threadIdx.x==0) {
#else
    if (true) {
#endif
      if (doalpha == -1) {
        sum *= value_type(-1);
      } else if (doalpha * doalpha != 1) {
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


  KOKKOS_INLINE_FUNCTION void
  operator() (const team_member& dev) const
  {
    for (int loop = 0; loop < rows_per_thread; ++loop) {

      // NOTE (mfh 20 Mar 2015) Unfortunately, Kokkos::Vectorization
      // lacks a typedef for determining the type of the return value
      // of global_thread_rank().  I know that it returns int now, but
      // this may change at some point.
      //
      // iRow represents a row of the matrix, so its correct type is
      // ordinal_type.

      const ordinal_type iRow = (dev.league_rank() * dev.team_size() + dev.team_rank())
                                * rows_per_thread + loop;
      if (iRow >= m_A.numRows ()) {
        return;
      }

      // mfh 20 Mar 2015: This relates to n, so its correct type is
      // ordinal_type.  Once we can use C++11 without protection, the
      // right thing to do would be to use decltype to pick up n's
      // type here, rather than assuming that it's ordinal_type.
      ordinal_type kk = 0;

#ifdef KOKKOS_FAST_COMPILE
      for (; kk + 4 <= n; kk += 4) {
        strip_mine<4>(dev, iRow, kk);
      }
      for( ; kk < n; ++kk) {
        strip_mine<1>(dev, iRow, kk);
      }
#else
#  ifdef __CUDA_ARCH__
      if ((n > 8) && (n % 8 == 1)) {
        strip_mine<9>(dev, iRow, kk);
        kk += 9;
      }
      for(; kk + 8 <= n; kk += 8)
        strip_mine<8>(dev, iRow, kk);
      if(kk < n)
        switch(n - kk) {
#  else // NOT a CUDA device
          if ((n > 16) && (n % 16 == 1)) {
            strip_mine<17>(dev, iRow, kk);
            kk += 17;
          }

          for (; kk + 16 <= n; kk += 16) {
            strip_mine<16>(dev, iRow, kk);
          }

          if(kk < n)
            switch(n - kk) {
            case 15:
              strip_mine<15>(dev, iRow, kk);
              break;

            case 14:
              strip_mine<14>(dev, iRow, kk);
              break;

            case 13:
              strip_mine<13>(dev, iRow, kk);
              break;

            case 12:
              strip_mine<12>(dev, iRow, kk);
              break;

            case 11:
              strip_mine<11>(dev, iRow, kk);
              break;

            case 10:
              strip_mine<10>(dev, iRow, kk);
              break;

            case 9:
              strip_mine<9>(dev, iRow, kk);
              break;

            case 8:
              strip_mine<8>(dev, iRow, kk);
              break;
#  endif // __CUDA_ARCH__
            case 7:
              strip_mine<7>(dev, iRow, kk);
              break;

            case 6:
              strip_mine<6>(dev, iRow, kk);
              break;

            case 5:
              strip_mine<5>(dev, iRow, kk);
              break;

            case 4:
              strip_mine<4>(dev, iRow, kk);
              break;

            case 3:
              strip_mine<3>(dev, iRow, kk);
              break;

            case 2:
              strip_mine<2>(dev, iRow, kk);
              break;

            case 1:
              strip_mine_1(dev, iRow);
              break;
            }
#endif // KOKKOS_FAST_COMPILE
        }
    }
  };


  template<class RangeVector,
           class CrsMatrix,
           class DomainVector,
           class CoeffVector1,
           class CoeffVector2,
           int doalpha,
           int dobeta>
  struct MV_MultiplySingleFunctor {
    typedef typename CrsMatrix::execution_space      execution_space;
    typedef typename CrsMatrix::ordinal_type         ordinal_type;
    typedef typename CrsMatrix::non_const_value_type value_type;
    typedef typename Kokkos::View<value_type*, typename CrsMatrix::execution_space> range_values;
    typedef typename Kokkos::TeamPolicy<execution_space> team_policy;
    typedef typename team_policy::member_type            team_member;

    CoeffVector1 beta;
    CoeffVector2 alpha;
    CrsMatrix  m_A;
    DomainVector m_x;
    RangeVector m_y;

    ordinal_type rows_per_thread;

    MV_MultiplySingleFunctor (const CoeffVector1 beta_,
                              const CoeffVector2 alpha_,
                              const CrsMatrix m_A_,
                              const DomainVector m_x_,
                              const RangeVector m_y_,
                              const int rows_per_thread_) :
      beta (beta_), alpha (alpha_),
      m_A (m_A_), m_x (m_x_), m_y (m_y_),
      rows_per_thread (rows_per_thread_)
    {}

    KOKKOS_INLINE_FUNCTION void
    operator() (const team_member& dev) const
    {
      // This should be a thread loop as soon as we can use C++11
      // FIXME (mfh 20 Mar 2015) The correct type of 'loop' should be
      // ordinal_type, not int.  Ditto for rows_per_thread.
      for (int loop = 0; loop < rows_per_thread; ++loop) {
        // NOTE (mfh 20 Mar 2015) Unfortunately, Kokkos::Vectorization
        // lacks a typedef for determining the type of the return
        // value of global_thread_rank().  I know that it returns int
        // now, but this may change at some point.
        //
        // iRow represents a row of the matrix, so its correct type is
        // ordinal_type.
        const ordinal_type iRow = (dev.league_rank() * dev.team_size() + dev.team_rank())
                                  * rows_per_thread + loop;
        if (iRow >= m_A.numRows ()) {
          return;
        }
        const SparseRowViewConst<CrsMatrix> row = m_A.rowConst(iRow);
        const ordinal_type row_length = static_cast<ordinal_type> (row.length);
        value_type sum = 0;

        // Use explicit Cuda below to avoid C++11 for now. This should be a vector reduce loop !
        #ifdef KOKKOS_HAVE_PRAGMA_IVDEP
        #pragma ivdep
        #endif
        #ifdef KOKKOS_HAVE_PRAGMA_UNROLL
        #pragma unroll
        #endif
        #ifdef KOKKOS_HAVE_PRAGMA_LOOPCOUNT
        #pragma loop count (15)
        #endif
#ifdef __CUDA_ARCH__
        for (ordinal_type iEntry = static_cast<ordinal_type> (threadIdx.x);
             iEntry < static_cast<ordinal_type> (row_length);
             iEntry += static_cast<ordinal_type> (blockDim.x)) {
#else
        for (ordinal_type iEntry = 0;
             iEntry < static_cast<ordinal_type> (row_length);
             iEntry ++) {
#endif
          sum += row.value(iEntry) * m_x(row.colidx(iEntry));
        }

#ifdef __CUDA_ARCH__
        if (blockDim.x > 1)
          sum += shfl_down(sum, 1,blockDim.x);
        if (blockDim.x > 2)
          sum += shfl_down(sum, 2,blockDim.x);
        if (blockDim.x > 4)
          sum += shfl_down(sum, 4,blockDim.x);
        if (blockDim.x > 8)
          sum += shfl_down(sum, 8,blockDim.x);
        if (blockDim.x > 16)
          sum += shfl_down(sum, 16,blockDim.x);

        if (threadIdx.x==0) {
#else
        if (true) {
#endif
          if (doalpha == -1) {
            sum *= value_type(-1);
          } else if (doalpha * doalpha != 1) {
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
        msg << "\t Labels are: y(" << y.tracker().label() << ") b("
            << betav.tracker().label() << ") a("
            << alphav.tracker().label() << ") x("
            << A.values.tracker().label() << ") x("
            << x.tracker().label() << ")\n";
        msg << "\t Dimensions are: y(" << y.dimension_0() << "," << y.dimension_1() << ") x(" << x.dimension_0() << "," << x.dimension_1() << ")\n";
        Impl::throw_runtime_exception( msg.str() );
      }
      if (numRows > y.dimension_0()) {
        std::ostringstream msg;
        msg << "Error in CRSMatrix - Vector Multiply (y = by + aAx): dimensions of y and A do not match\n";
        msg << "\t Labels are: y(" << y.tracker().label() << ") b("
            << betav.tracker().label() << ") a("
            << alphav.tracker().label() << ") x("
            << A.values.tracker().label() << ") x("
            << x.tracker().label() << ")\n";
        msg << "\t Dimensions are: y(" << y.dimension_0() << "," << y.dimension_1() << ") A(" << A.numRows() << "," << A.numCols() << ")\n";
        Impl::throw_runtime_exception( msg.str() );
      }
      if (numCols > x.dimension_0()) {
        std::ostringstream msg;
        msg << "Error in CRSMatrix - Vector Multiply (y = by + aAx): dimensions of x and A do not match\n";
        msg << "\t Labels are: y(" << y.tracker().label() << ") b("
            << betav.tracker().label() << ") a("
            << alphav.tracker().label() << ") x("
            << A.values.tracker().label() << ") x("
            << x.tracker().label() << ")\n";
        msg << "\t Dimensions are: x(" << x.dimension_0() << "," << x.dimension_1() << ") A(" << A.numRows() << "," << A.numCols() << ")\n";
        Impl::throw_runtime_exception( msg.str() );
      }
      if (dobeta==2) {
        if (betav.dimension_0()!=numVecs) {
          std::ostringstream msg;
          msg << "Error in CRSMatrix - Vector Multiply (y = by + aAx): 2nd dimensions of y and b do not match\n";
          msg << "\t Labels are: y(" << y.tracker().label() << ") b("
              << betav.tracker().label() << ") a("
              << alphav.tracker().label() << ") x("
              << A.values.tracker().label() << ") x("
              << x.tracker().label() << ")\n";
          msg << "\t Dimensions are: y(" << y.dimension_0() << "," << y.dimension_1() << ") b(" << betav.dimension_0() << ")\n";
          Impl::throw_runtime_exception( msg.str() );
        }
      }
      if(doalpha==2) {
        if(alphav.dimension_0()!=numVecs) {
          std::ostringstream msg;
          msg << "Error in CRSMatrix - Vector Multiply (y = by + aAx): 2nd dimensions of x and b do not match\n";
          msg << "\t Labels are: y(" << y.tracker().label() << ") b("
              << betav.tracker().label() << ") a("
              << alphav.tracker().label() << ") x("
              << A.values.tracker().label() << ") x("
              << x.tracker().label() << ")\n";
          msg << "\t Dimensions are: x(" << x.dimension_0() << "," << x.dimension_1() << ") b(" << betav.dimension_0() << ")\n";
          Impl::throw_runtime_exception( msg.str() );
        }
      }
    }
  } // namespace Impl

  // This TransposeFunctor is functional, but not necessarily performant.
  template<class RangeVector,
           class CrsMatrix,
           class DomainVector,
           class CoeffVector1,
           class CoeffVector2,
           int doalpha,
           int dobeta,
           bool conjugate = false,
           int NNZPerRow = 27>
  struct MV_MultiplyTransposeFunctor {
    typedef typename CrsMatrix::execution_space                 execution_space;
    typedef typename CrsMatrix::ordinal_type                    ordinal_type;
    typedef typename CrsMatrix::size_type                       size_type;
    typedef typename CrsMatrix::non_const_value_type            value_type;
    typedef typename Kokkos::View<value_type*, execution_space> range_values;

    typedef MV_MultiplyShflThreadsPerRow<execution_space, value_type, NNZPerRow> ShflThreadsPerRow;

    CoeffVector1 beta;
    CoeffVector2 alpha;
    CrsMatrix  m_A ;
    DomainVector  m_x ;
    RangeVector  m_y ;
    ordinal_type n;

    // This is an iteration over rows of the matrix (modulo the
    // shuffle width), so the correct type of i is ordinal_type.
    KOKKOS_INLINE_FUNCTION
    void operator() (const ordinal_type i) const {
      typedef Kokkos::Details::ArithTraits<value_type> ATV;

      const ordinal_type iRow = i / ShflThreadsPerRow::device_value;
      const int lane = static_cast<int> (i) % ShflThreadsPerRow::device_value;
      const SparseRowViewConst<CrsMatrix> row = m_A.rowConst(iRow);

      for (ordinal_type iEntry = static_cast<ordinal_type> (lane);
           iEntry < static_cast<ordinal_type> (row.length);
           iEntry += static_cast<ordinal_type> (ShflThreadsPerRow::device_value)) {
        const value_type val = conjugate ?
          ATV::conj (row.value(iEntry)) :
          row.value(iEntry);
        const ordinal_type ind = row.colidx(iEntry);

        if (doalpha != 1) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
          for (ordinal_type k = 0; k < n; ++k) {
            atomic_add (&m_y(ind,k), value_type(alpha(k) * val * m_x(iRow, k)));
          }
        } else {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
          for (ordinal_type k = 0; k < n; ++k) {
            atomic_add (&m_y(ind,k), value_type(val * m_x(iRow, k)));
          }
        }
      }
    }
  };

  // This TansposeFunctor is functional, but not necessarily performant.
  template<class RangeVector,
           class CrsMatrix,
           class DomainVector,
           class CoeffVector1,
           class CoeffVector2,
           int doalpha,
           int dobeta,
           bool conjugate = false,
           int NNZPerRow = 27 >
  struct MV_MultiplyTransposeSingleFunctor {
    typedef typename CrsMatrix::execution_space      execution_space;
    typedef typename CrsMatrix::ordinal_type         ordinal_type;
    typedef typename CrsMatrix::non_const_value_type value_type;
    typedef typename Kokkos::View<value_type*, execution_space> range_values;

    typedef MV_MultiplyShflThreadsPerRow< execution_space , value_type , NNZPerRow > ShflThreadsPerRow ;

    CoeffVector1 beta;
    CoeffVector2 alpha;
    CrsMatrix  m_A ;
    DomainVector  m_x ;
    RangeVector  m_y ;
    ordinal_type n;

    KOKKOS_INLINE_FUNCTION
    void operator() (const ordinal_type i) const {
      typedef Kokkos::Details::ArithTraits<value_type> ATV;

      const ordinal_type iRow = i / ShflThreadsPerRow::device_value;
      const int lane = static_cast<int> (i) % ShflThreadsPerRow::device_value;
      const SparseRowViewConst<CrsMatrix> row = m_A.rowConst(iRow);

      for (ordinal_type iEntry = static_cast<ordinal_type> (lane);
           iEntry < static_cast<ordinal_type> (row.length);
           iEntry += static_cast<ordinal_type> (ShflThreadsPerRow::device_value)) {
        const value_type val = conjugate ?
          ATV::conj (row.value(iEntry)) :
          row.value(iEntry);
        const ordinal_type ind = row.colidx(iEntry);

        if (doalpha != 1) {
          atomic_add (&m_y(ind), value_type(alpha(0) * val * m_x(iRow)));
        } else {
          atomic_add (&m_y(ind), value_type(val * m_x(iRow)));
        }
      }
    }
  };

template <class RangeVector,
          class TCrsMatrix,
          class DomainVector,
          class CoeffVector1,
          class CoeffVector2,
          int doalpha,
          int dobeta>
void
MV_MultiplyTranspose (typename Kokkos::Impl::enable_if<DomainVector::Rank == 2, const CoeffVector1>::type& betav,
                      const RangeVector &y,
                      const CoeffVector2 &alphav,
                      const TCrsMatrix &A,
                      const DomainVector &x,
                      const bool conjugate = false)
{
  typedef typename TCrsMatrix::ordinal_type ordinal_type;

  // FIXME (mfh 02 Jan 2015) Is numRows() always signed?  More
  // importantly, if the calling process owns zero rows in the row
  // Map, numRows() should return 0, not -1.
  //
  //Special case for zero Rows RowMap
  if (A.numRows () == static_cast<ordinal_type> (-1)) {
    return;
  }

  if (doalpha == 0) {
    if (dobeta == 2) {
      MV_MulScalar (y, betav, y);
    } else {
      MV_MulScalar (y, static_cast<typename RangeVector::const_value_type> (dobeta), y);
    }
    return;
  } else {
    typedef View< typename RangeVector::non_const_data_type ,
                  typename RangeVector::array_layout ,
                  typename RangeVector::execution_space ,
                  typename RangeVector::memory_traits >
    RangeVectorType;

    typedef View< typename DomainVector::const_data_type ,
                  typename DomainVector::array_layout ,
                  typename DomainVector::execution_space ,
                  Kokkos::MemoryRandomAccess >
    DomainVectorType;

    typedef View< typename CoeffVector1::const_data_type ,
                  typename CoeffVector1::array_layout ,
                  typename CoeffVector1::execution_space ,
                  Kokkos::MemoryRandomAccess >
    CoeffVector1Type;

    typedef View< typename CoeffVector2::const_data_type ,
                  typename CoeffVector2::array_layout ,
                  typename CoeffVector2::execution_space ,
                  Kokkos::MemoryRandomAccess >
    CoeffVector2Type;

    typedef CrsMatrix<typename TCrsMatrix::const_value_type,
                      typename TCrsMatrix::ordinal_type,
                      typename TCrsMatrix::execution_space,
                      typename TCrsMatrix::memory_traits,
                      typename TCrsMatrix::size_type> CrsMatrixType;

    //Impl::MV_Multiply_Check_Compatibility(betav,y,alphav,A,x,doalpha,dobeta);
/*
#ifndef KOKKOS_FAST_COMPILE

    if(x.dimension_1()==1) {
      typedef View<typename DomainVectorType::const_value_type*,typename DomainVector::array_layout ,typename DomainVectorType::execution_space,Kokkos::MemoryRandomAccess> DomainVector1D;
      typedef View<typename DomainVectorType::const_value_type*,typename DomainVector::array_layout ,typename DomainVectorType::execution_space> DomainVector1DPlain;
      typedef View<typename RangeVectorType::value_type*,typename RangeVector::array_layout ,typename RangeVectorType::execution_space,typename RangeVector::memory_traits> RangeVector1D;

       Kokkos::subview( y , ALL(),0 );

       if (conjugate) {
         typedef MV_MultiplySingleFunctor<RangeVector1D, CrsMatrixType, DomainVector1D,
           CoeffVector1Type, CoeffVector2Type, doalpha, dobeta, true> OpType;
         OpType op;
         const typename CrsMatrixType::ordinal_type nrow = A.numRows ();
         op.m_A = A;
         op.m_x = Kokkos::subview (x, ALL (), 0);
         op.m_y = Kokkos::subview (y, ALL(), 0);
         op.beta = betav;
         op.alpha = alphav;
         op.n = x.dimension(1);
         Kokkos::parallel_for (nrow * OpType::ShflThreadsPerRow::host_value(), op);
      }
      else {
         typedef MV_MultiplySingleFunctor<RangeVector1D, CrsMatrixType, DomainVector1D,
           CoeffVector1Type, CoeffVector2Type, doalpha, dobeta, false> OpType;
         OpType op;
         const typename CrsMatrixType::ordinal_type nrow = A.numRows ();
         op.m_A = A;
         op.m_x = Kokkos::subview (x, ALL (), 0);
         op.m_y = Kokkos::subview (y, ALL(), 0);
         op.beta = betav;
         op.alpha = alphav;
         op.n = x.dimension(1);
         Kokkos::parallel_for (nrow * OpType::ShflThreadsPerRow::host_value(), op);
      }
    }
    else {
      if (conjugate) {
        typedef MV_MultiplyFunctor<RangeVectorType, CrsMatrixType, DomainVectorType,
          CoeffVector1Type, CoeffVector2Type, doalpha, dobeta, true> OpType ;
        OpType op ;
        const typename CrsMatrixType::ordinal_type nrow = A.numRows();
        op.m_A = A ;
        op.m_x = x ;
        op.m_y = y ;
        op.beta = betav;
        op.alpha = alphav;
        op.n = x.dimension(1);
        Kokkos::parallel_for(nrow*OpType::ShflThreadsPerRow::host_value() , op);
      }
      else {
        typedef MV_MultiplyFunctor<RangeVectorType, CrsMatrixType, DomainVectorType,
          CoeffVector1Type, CoeffVector2Type, doalpha, dobeta, false> OpType ;
        OpType op ;
        const typename CrsMatrixType::ordinal_type nrow = A.numRows();
        op.m_A = A ;
        op.m_x = x ;
        op.m_y = y ;
        op.beta = betav;
        op.alpha = alphav;
        op.n = x.dimension(1);
        Kokkos::parallel_for(nrow*OpType::ShflThreadsPerRow::host_value() , op);
      }
    }

#else // NOT KOKKOS_FAST_COMPILE
*/

    // FIXME (mfh 20 Mar 2015) The dimension doesn't have type int.
    // On the other hand, this is the number of columns in the input
    // (and output) MultiVectors, so it's not likely to be large (the
    // typical case is < 100).
    int numVecs = x.dimension_1();
    CoeffVector1 beta = betav;
    CoeffVector2 alpha = alphav;

    if (doalpha != 2) {
      alpha = CoeffVector2("CrsMatrix::auto_a", numVecs);
      typename CoeffVector2::HostMirror h_a = Kokkos::create_mirror_view(alpha);
      typename CoeffVector2::value_type s_a = (typename CoeffVector2::value_type) doalpha;

      for (int i = 0; i < numVecs; ++i) {
        h_a(i) = s_a;
      }

      Kokkos::deep_copy (alpha, h_a);
    }

    if (dobeta != 2) {
      beta = CoeffVector1("CrsMatrix::auto_b", numVecs);
      typename CoeffVector1::HostMirror h_b = Kokkos::create_mirror_view(beta);
      typename CoeffVector1::value_type s_b = (typename CoeffVector1::value_type) dobeta;

      for (int i = 0; i < numVecs; ++i) {
        h_b(i) = s_b;
      }

      Kokkos::deep_copy (beta, h_b);
    }

    if (dobeta == 2) {
      MV_MulScalar (y, betav, y);
    } else {
      if (dobeta != 1) {
        MV_MulScalar (y, static_cast<typename RangeVector::const_value_type> (dobeta), y);
      }
    }

    const ordinal_type nrow = A.numRows ();

    if (conjugate) {
      typedef MV_MultiplyTransposeFunctor<RangeVectorType, CrsMatrixType,
                                          DomainVectorType, CoeffVector1Type,
                                          CoeffVector2Type, 2, 2, true> OpType;
      OpType op ;
      op.m_A = A;
      op.m_x = x;
      op.m_y = y;
      op.beta = beta;
      op.alpha = alpha;
      op.n = x.dimension_1();
      Kokkos::parallel_for (nrow * OpType::ShflThreadsPerRow::host_value (), op);
    }
    else {
      typedef MV_MultiplyTransposeFunctor<RangeVectorType, CrsMatrixType,
                                          DomainVectorType, CoeffVector1Type,
                                          CoeffVector2Type, 2, 2, false> OpType;
      OpType op ;
      op.m_A = A;
      op.m_x = x;
      op.m_y = y;
      op.beta = beta;
      op.alpha = alpha;
      op.n = x.dimension_1();
      Kokkos::parallel_for (nrow * OpType::ShflThreadsPerRow::host_value (), op);
    }

//#endif // KOKKOS_FAST_COMPILE
  }
}

template<class RangeVector,
         class CrsMatrix,
         class DomainVector,
         class CoeffVector1,
         class CoeffVector2>
void
MV_MultiplyTranspose (const CoeffVector1& betav,
                      const RangeVector& y,
                      const CoeffVector2& alphav,
                      const CrsMatrix& A,
                      const DomainVector& x,
                      int beta,
                      int alpha,
                      const bool conjugate = false)
{
  if (beta == 0) {
    if (alpha == 0) {
      MV_MultiplyTranspose<RangeVector, CrsMatrix, DomainVector, CoeffVector1,
                           CoeffVector2, 0, 0 > (betav, y, alphav, A, x, conjugate);
    }
    else if (alpha == 1) {
      MV_MultiplyTranspose<RangeVector, CrsMatrix, DomainVector, CoeffVector1,
                           CoeffVector2, 1, 0 > (betav, y, alphav, A, x, conjugate);
    }
    else if (alpha == -1) {
      MV_MultiplyTranspose<RangeVector, CrsMatrix, DomainVector, CoeffVector1,
                           CoeffVector2, -1, 0 > (betav, y, alphav, A, x, conjugate);
    }
    else {
      MV_MultiplyTranspose<RangeVector, CrsMatrix, DomainVector, CoeffVector1,
                           CoeffVector2, 2, 0 > (betav, y, alphav, A, x, conjugate);
    }
  } else if (beta == 1) {
    if (alpha == 0) {
      return;
    }
    else if (alpha == 1) {
      MV_MultiplyTranspose<RangeVector, CrsMatrix, DomainVector, CoeffVector1,
                           CoeffVector2, 1, 1 > (betav, y, alphav, A, x, conjugate);
    }
    else if (alpha == -1) {
      MV_MultiplyTranspose<RangeVector, CrsMatrix, DomainVector, CoeffVector1,
                           CoeffVector2, -1, 1 > (betav, y, alphav, A, x, conjugate);
    }
    else {
      MV_MultiplyTranspose<RangeVector, CrsMatrix, DomainVector, CoeffVector1,
                           CoeffVector2, 2, 1 > (betav, y, alphav, A, x, conjugate);
    }
  } else if (beta == -1) {
    if (alpha == 0) {
      MV_MultiplyTranspose<RangeVector, CrsMatrix, DomainVector, CoeffVector1,
                           CoeffVector2, 0, -1 > (betav, y, alphav, A, x, conjugate);
    }
    else if (alpha == 1) {
      MV_MultiplyTranspose<RangeVector, CrsMatrix, DomainVector, CoeffVector1,
                           CoeffVector2, 1, -1 > (betav, y, alphav, A, x, conjugate);
    }
    else if (alpha == -1) {
      MV_MultiplyTranspose<RangeVector, CrsMatrix, DomainVector, CoeffVector1,
                           CoeffVector2, -1, -1 > (betav, y, alphav, A, x, conjugate);
    }
    else {
      MV_MultiplyTranspose<RangeVector, CrsMatrix, DomainVector, CoeffVector1,
                           CoeffVector2, 2, -1 > (betav, y, alphav, A, x, conjugate);
    }
  } else {
    if (alpha == 0) {
      MV_MultiplyTranspose<RangeVector, CrsMatrix, DomainVector, CoeffVector1,
                           CoeffVector2, 0, 2 > (betav, y, alphav, A, x, conjugate);
    }
    else if (alpha == 1) {
      MV_MultiplyTranspose<RangeVector, CrsMatrix, DomainVector, CoeffVector1,
                           CoeffVector2, 1, 2 > (betav, y, alphav, A, x, conjugate);
    }
    else if (alpha == -1) {
      MV_MultiplyTranspose<RangeVector, CrsMatrix, DomainVector, CoeffVector1,
                           CoeffVector2, -1, 2 > (betav, y, alphav, A, x, conjugate);
    }
    else {
      MV_MultiplyTranspose<RangeVector, CrsMatrix, DomainVector, CoeffVector1,
                           CoeffVector2, 2, 2> (betav, y, alphav, A, x, conjugate);
    }
  }
}

template<class RangeVector, class CrsMatrix, class DomainVector>
void
MV_MultiplyTranspose (typename RangeVector::const_value_type s_b,
                      const RangeVector& y,
                      typename DomainVector::const_value_type s_a,
                      const CrsMatrix& A,
                      const DomainVector& x,
                      const bool conjugate = false)
{
/*#ifdef KOKKOS_USE_CUSPARSE
  if (MV_Multiply_Try_CuSparse (s_b, y, s_a, A, x, conjugate)) {
    return;
  }
#endif // KOKKOSE_USE_CUSPARSE
#ifdef KOKKOS_USE_MKL
  if (MV_Multiply_Try_MKL (s_b, y, s_a, A, x, conjugate)) {
    return;
  }
#endif // KOKKOS_USE_MKL*/
  typedef Kokkos::View<typename RangeVector::value_type*,
                       typename RangeVector::execution_space> aVector;
  aVector a;
  aVector b;
  int numVecs = x.dimension_1();

  if (s_b == 0) {
    if (s_a == 0)
      return MV_MultiplyTranspose (a, y, a, A, x, 0, 0, conjugate);
    else if (s_a == 1)
      return MV_MultiplyTranspose (a, y, a, A, x, 0, 1, conjugate);
    else if (s_a == static_cast<typename DomainVector::const_value_type> (-1))
      return MV_MultiplyTranspose (a, y, a, A, x, 0, -1, conjugate);
    else {
      a = aVector("a", numVecs);
      typename aVector::HostMirror h_a = Kokkos::create_mirror_view (a);
      for (int i = 0; i < numVecs; ++i) {
        h_a(i) = s_a;
      }
      Kokkos::deep_copy (a, h_a);
      return MV_MultiplyTranspose (a, y, a, A, x, 0, 2, conjugate);
    }
  } else if (s_b == 1) {
    if (s_a == 0)
      return MV_MultiplyTranspose (a, y, a, A, x, 1, 0, conjugate);
    else if (s_a == 1)
      return MV_MultiplyTranspose (a, y, a, A, x, 1, 1, conjugate);
    else if (s_a == static_cast<typename DomainVector::const_value_type> (-1))
      return MV_MultiplyTranspose (a, y, a, A, x, 1, -1, conjugate);
    else {
      a = aVector("a", numVecs);
      typename aVector::HostMirror h_a = Kokkos::create_mirror_view (a);
      for (int i = 0; i < numVecs; ++i) {
        h_a(i) = s_a;
      }
      Kokkos::deep_copy (a, h_a);
      return MV_MultiplyTranspose (a, y, a, A, x, 1, 2, conjugate);
    }
  } else if (s_b == static_cast<typename RangeVector::const_value_type> (-1)) {
    if (s_a == 0)
      return MV_MultiplyTranspose (a, y, a, A, x, -1, 0, conjugate);
    else if (s_a == 1)
      return MV_MultiplyTranspose (a, y, a, A, x, -1, 1, conjugate);
    else if (s_a == static_cast<typename DomainVector::const_value_type> (-1))
      return MV_MultiplyTranspose (a, y, a, A, x, -1, -1, conjugate);
    else {
      a = aVector("a", numVecs);
      typename aVector::HostMirror h_a = Kokkos::create_mirror_view (a);
      for (int i = 0; i < numVecs; ++i) {
        h_a(i) = s_a;
      }
      Kokkos::deep_copy (a, h_a);
      return MV_MultiplyTranspose (a, y, a, A, x, -1, 2, conjugate);
    }
  } else {
    b = aVector("b", numVecs);
    typename aVector::HostMirror h_b = Kokkos::create_mirror_view (b);
    for (int i = 0; i < numVecs; ++i) {
      h_b(i) = s_b;
    }
    Kokkos::deep_copy (b, h_b);

    if (s_a == 0)
      return MV_MultiplyTranspose (b, y, a, A, x, 2, 0, conjugate);
    else if (s_a == 1)
      return MV_MultiplyTranspose (b, y, a, A, x, 2, 1, conjugate);
    else if (s_a == static_cast<typename DomainVector::const_value_type> (-1))
      return MV_MultiplyTranspose (b, y, a, A, x, 2, -1, conjugate);
    else {
      a = aVector("a", numVecs);
      typename aVector::HostMirror h_a = Kokkos::create_mirror_view (a);
      for (int i = 0; i < numVecs; ++i) {
        h_a(i) = s_a;
      }
      Kokkos::deep_copy (a, h_a);
      return MV_MultiplyTranspose (b, y, a, A, x, 2, 2, conjugate);
    }
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
  MV_MultiplySingle (typename Kokkos::Impl::enable_if<DomainVector::Rank == 1, const CoeffVector1>::type& betav,
               const RangeVector &y,
               const CoeffVector2 &alphav,
               const TCrsMatrix& A,
               const DomainVector& x)
  {
    typedef typename TCrsMatrix::ordinal_type ordinal_type;

    if (A.numRows () <= static_cast<ordinal_type> (0)) {
      return;
    }
    if (doalpha == 0) {
      if (dobeta==2) {
        V_MulScalar (y, betav, y);
      }
      else {
        V_MulScalar (y, typename RangeVector::value_type (dobeta), y);
      }
      return;
    } else {
      typedef View< typename RangeVector::non_const_data_type ,
                    typename RangeVector::array_layout ,
                    typename RangeVector::execution_space ,
                    typename RangeVector::memory_traits >
      RangeVectorType;

      typedef View< typename DomainVector::const_data_type ,
                    typename DomainVector::array_layout ,
                    typename DomainVector::execution_space ,
                    //typename DomainVector::memory_traits >
                    Kokkos::MemoryRandomAccess >
      DomainVectorType;

      typedef View< typename CoeffVector1::const_data_type ,
                    typename CoeffVector1::array_layout ,
                    typename CoeffVector1::execution_space ,
                    Kokkos::MemoryRandomAccess >
      CoeffVector1Type;

      typedef View< typename CoeffVector2::const_data_type ,
                    typename CoeffVector2::array_layout ,
                    typename CoeffVector2::execution_space ,
                    Kokkos::MemoryRandomAccess >
      CoeffVector2Type;

      typedef CrsMatrix<typename TCrsMatrix::const_value_type,
                        typename TCrsMatrix::ordinal_type,
                        typename TCrsMatrix::execution_space,
                        typename TCrsMatrix::memory_traits,
                        typename TCrsMatrix::size_type>
      CrsMatrixType;
      typedef typename CrsMatrixType::size_type size_type;

      Impl::MV_Multiply_Check_Compatibility(betav,y,alphav,A,x,doalpha,dobeta);

      // NNZPerRow could be anywhere from 0, to A.numRows()*A.numCols().
      // Thus, the appropriate type is size_type.
      const size_type NNZPerRow = A.nnz () / A.numRows ();

      int vector_length = 1;
      while( (static_cast<size_type> (vector_length*2*3) <= NNZPerRow) && (vector_length<32) ) vector_length*=2;

#ifndef KOKKOS_FAST_COMPILE // This uses templated fucntions on doalpha and dobeta and will produce 16 kernels

      typedef MV_MultiplySingleFunctor<RangeVectorType, CrsMatrixType, DomainVectorType,
                                       CoeffVector1Type, CoeffVector2Type, doalpha, dobeta > OpType ;

      const typename CrsMatrixType::ordinal_type nrow = A.numRows();

      OpType op(betav,alphav,A,x,y,RowsPerThread<typename RangeVector::execution_space >(NNZPerRow)) ;

      Kokkos::parallel_for( Kokkos::TeamPolicy< typename RangeVector::execution_space >
         ( nrow ,RowsPerThread<typename RangeVector::execution_space >(NNZPerRow), vector_length ) , op );

#else // NOT KOKKOS_FAST_COMPILE

      typedef MV_MultiplySingleFunctor<RangeVectorType, CrsMatrixType, DomainVectorType,
                               CoeffVector1Type, CoeffVector2Type, 2, 2> OpType ;

      int numVecs = x.dimension_1(); // == 1
      CoeffVector1 beta = betav;
      CoeffVector2 alpha = alphav;

      if(doalpha!=2) {
              alpha = CoeffVector2("CrsMatrix::auto_a", numVecs);
              typename CoeffVector2::HostMirror h_a = Kokkos::create_mirror_view(alpha);
              typename CoeffVector2::value_type s_a = (typename CoeffVector2::value_type) doalpha;

              for(int i = 0; i < numVecs; i++)
                h_a(i) = s_a;

              Kokkos::deep_copy(alpha, h_a);
      }
      if(dobeta!=2) {
              beta = CoeffVector1("CrsMatrix::auto_b", numVecs);
              typename CoeffVector1::HostMirror h_b = Kokkos::create_mirror_view(beta);
              typename CoeffVector1::value_type s_b = (typename CoeffVector1::value_type) dobeta;

              for(int i = 0; i < numVecs; i++)
                h_b(i) = s_b;

              Kokkos::deep_copy(beta, h_b);
      }

      const typename CrsMatrixType::ordinal_type nrow = A.numRows();

      OpType op(beta,alpha,A,x,y,RowsPerThread<typename RangeVector::execution_space >(NNZPerRow)) ;
      Kokkos::parallel_for( Kokkos::TeamPolicy< typename RangeVector::execution_space >
           ( nrow ,RowsPerThread<typename RangeVector::execution_space >(NNZPerRow), vector_length ) , op );

#endif // KOKKOS_FAST_COMPILE
    }
  }

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
    //Special case for zero Rows RowMap
    if(A.numRows() <= 0) return;

    if (doalpha == 0) {
      if (dobeta==2) {
              MV_MulScalar(y,betav,y);
      } else {
              MV_MulScalar(y,static_cast<typename RangeVector::const_value_type> (dobeta),y);
      }
      return;
    } else {
      typedef View< typename RangeVector::non_const_data_type ,
                    typename RangeVector::array_layout ,
                    typename RangeVector::execution_space ,
                    typename RangeVector::memory_traits >
      RangeVectorType;

      typedef View< typename DomainVector::const_data_type ,
                    typename DomainVector::array_layout ,
                    typename DomainVector::execution_space ,
                    Kokkos::MemoryRandomAccess >
      DomainVectorType;

      typedef View< typename CoeffVector1::const_data_type ,
                    typename CoeffVector1::array_layout ,
                    typename CoeffVector1::execution_space ,
                    Kokkos::MemoryRandomAccess >
      CoeffVector1Type;

      typedef View< typename CoeffVector2::const_data_type ,
                    typename CoeffVector2::array_layout ,
                    typename CoeffVector2::execution_space ,
                    Kokkos::MemoryRandomAccess >
      CoeffVector2Type;

      typedef CrsMatrix<typename TCrsMatrix::const_value_type,
                        typename TCrsMatrix::ordinal_type,
                        typename TCrsMatrix::execution_space,
                        typename TCrsMatrix::memory_traits,
                        typename TCrsMatrix::size_type> CrsMatrixType;

      Impl::MV_Multiply_Check_Compatibility(betav,y,alphav,A,x,doalpha,dobeta);

      const int NNZPerRow = A.nnz()/A.numRows();

#ifndef KOKKOS_FAST_COMPILE

      if(x.dimension_1()==1) {
        typedef View<typename DomainVectorType::const_value_type*,typename DomainVector::array_layout ,typename DomainVectorType::execution_space> DomainVector1D;
        typedef View<typename RangeVectorType::value_type*,typename RangeVector::array_layout ,typename RangeVectorType::execution_space,typename RangeVector::memory_traits> RangeVector1D;
        RangeVector1D y_sub = RangeVector1D(y.ptr_on_device(),y.dimension_0());
        DomainVector1D x_sub = DomainVector1D(x.ptr_on_device(),x.dimension_0());

        return MV_MultiplySingle<RangeVector1D,TCrsMatrix,DomainVector1D,CoeffVector1,CoeffVector2,doalpha,dobeta>
          (betav,y_sub,alphav,A,x_sub);

      } else {

        //Currently for multiple right hand sides its not worth it to use more than 8 threads per row on GPUs
        int vector_length = 1;
        while( (vector_length*2*3 <= NNZPerRow) && (vector_length<8) ) vector_length*=2;
        typedef MV_MultiplyFunctor<RangeVectorType, CrsMatrixType, DomainVectorType,
          CoeffVector1Type, CoeffVector2Type, doalpha, dobeta> OpType ;

        OpType op(betav,alphav,A,x,y,x.dimension_1(),RowsPerThread<typename RangeVector::execution_space >(NNZPerRow)) ;
        Kokkos::parallel_for( Kokkos::TeamPolicy< typename RangeVector::execution_space >
             ( A.numRows() ,RowsPerThread<typename RangeVector::execution_space >(NNZPerRow), vector_length ) , op );

      }

#else // NOT KOKKOS_FAST_COMPILE

      //Currently for multiple right hand sides its not worth it to use more than 8 threads per row on GPUs
      int vector_length = 1;
      while( (vector_length*2*3 <= NNZPerRow) && (vector_length<8) ) vector_length*=2;
      typedef MV_MultiplyFunctor<RangeVectorType, CrsMatrixType, DomainVectorType,
        CoeffVector1Type, CoeffVector2Type, 2, 2> OpType ;

      int numVecs = x.dimension_1();
      CoeffVector1 beta = betav;
      CoeffVector2 alpha = alphav;

      if (doalpha != 2) {
              alpha = CoeffVector2("CrsMatrix::auto_a", numVecs);
              typename CoeffVector2::HostMirror h_a = Kokkos::create_mirror_view(alpha);
              typename CoeffVector2::value_type s_a = (typename CoeffVector2::value_type) doalpha;

              for (int i = 0; i < numVecs; ++i)
                h_a(i) = s_a;

              Kokkos::deep_copy(alpha, h_a);
      }

      if (dobeta != 2) {
              beta = CoeffVector1("CrsMatrix::auto_b", numVecs);
              typename CoeffVector1::HostMirror h_b = Kokkos::create_mirror_view(beta);
              typename CoeffVector1::value_type s_b = (typename CoeffVector1::value_type) dobeta;

              for(int i = 0; i < numVecs; i++)
                h_b(i) = s_b;

              Kokkos::deep_copy(beta, h_b);
      }

      const typename CrsMatrixType::ordinal_type nrow = A.numRows();

      OpType op(betav,alphav,A,x,y,x.dimension_1(),RowsPerThread<typename RangeVector::execution_space >(NNZPerRow)) ;
      Kokkos::parallel_for( Kokkos::TeamPolicy< typename RangeVector::execution_space >
           ( A.numRows() ,RowsPerThread<typename RangeVector::execution_space >(NNZPerRow), vector_length ) , op );

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
               const DomainVector& x) {
    return MV_MultiplySingle<RangeVector,TCrsMatrix,DomainVector,CoeffVector1,CoeffVector2,doalpha,dobeta>
           (betav,y,alphav,A,x);
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
    if (CuSparse::MV_Multiply_Try_CuSparse (0.0, y, 1.0, A, x)) {
      return;
    }
#endif // KOKKOS_USE_CUSPARSE
#ifdef KOKKOS_USE_MKL
    if (MV_Multiply_Try_MKL (0.0, y, 1.0, A, x)) {
      return;
    }
#endif // KOKKOS_USE_MKL
    typedef Kokkos::View<typename DomainVector::value_type*, typename DomainVector::execution_space> aVector;
    aVector a;

    return MV_Multiply (a, y, a, A, x, 0, 1);
  }

  template<class RangeVector, class CrsMatrix, class DomainVector>
  void
  MV_Multiply (const RangeVector& y,
               typename DomainVector::const_value_type s_a,
               const CrsMatrix& A,
               const DomainVector& x)
  {
#ifdef KOKKOS_USE_CUSPARSE
    if (CuSparse::MV_Multiply_Try_CuSparse (0.0, y, s_a, A, x)) {
      return;
    }
#endif // KOKKOS_USE_CUSPARSE
#ifdef KOKKOS_USE_MKL
    if (MV_Multiply_Try_MKL (0.0, y, s_a, A, x)) {
      return;
    }
#endif // KOKKOS_USE_MKL
    typedef Kokkos::View<typename RangeVector::value_type*, typename RangeVector::execution_space> aVector;
    aVector a;
    const int numVecs = x.dimension_1();

    //if ((s_a < 1) && (s_a != 0)) {
    if (s_a == -1.0) {
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
  MV_Multiply (typename RangeVector::const_value_type s_b,
               const RangeVector& y,
               typename DomainVector::const_value_type s_a,
               const CrsMatrix& A,
               const DomainVector& x)
  {
#ifdef KOKKOS_USE_CUSPARSE
    if (CuSparse::MV_Multiply_Try_CuSparse (s_b, y, s_a, A, x)) {
      return;
    }
#endif // KOKKOSE_USE_CUSPARSE
#ifdef KOKKOS_USE_MKL
    if (MV_Multiply_Try_MKL (s_b, y, s_a, A, x)) {
      return;
    }
#endif // KOKKOS_USE_MKL
    typedef Kokkos::View<typename RangeVector::value_type*, typename RangeVector::execution_space> aVector;
    aVector a;
    aVector b;
    int numVecs = x.dimension_1();

    // [HCE 2013/12/09] Following 'if' appears to be a mistake and has been commented out
    // if(numVecs==1)
    if (s_b == 0) {
      if (s_a == 0)
        return MV_Multiply (a, y, a, A, x, 0, 0);
      else if (s_a == 1)
        return MV_Multiply (a, y, a, A, x, 0, 1);
      else if (s_a == static_cast<typename DomainVector::const_value_type> (-1))
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
      else if (s_a == static_cast<typename DomainVector::const_value_type> (-1))
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
    } else if (s_b == static_cast<typename RangeVector::const_value_type> (-1)) {
      if (s_a == 0)
        return MV_Multiply (a, y, a, A, x, -1, 0);
      else if (s_a == 1)
        return MV_Multiply (a, y, a, A, x, -1, 1);
      else if (s_a == static_cast<typename DomainVector::const_value_type> (-1))
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
      else if (s_a == static_cast<typename DomainVector::const_value_type> (-1))
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
