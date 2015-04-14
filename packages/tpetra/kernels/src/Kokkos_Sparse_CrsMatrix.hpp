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

#ifndef KOKKOS_SPARSE_CRSMATRIX_HPP_
#define KOKKOS_SPARSE_CRSMATRIX_HPP_

#include <Kokkos_Core.hpp>
#include <Kokkos_StaticCrsGraph.hpp>

#ifdef KOKKOS_HAVE_CXX11
#include <type_traits>
#endif // KOKKOS_HAVE_CXX11


namespace KokkosSparse {

// Blas like Transpose / Conjugate attributes for sparse kernels
static char Transpose[] = "T";
static char Conjugate[] = "C";
static char ConjugateTranspose[] = "H";
static char NoTranspose[] = "N";

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
  //! Type of the matrix's memory space.
  typedef typename Device::memory_space memory_space;
  //! Type of the matrix's device type.
  typedef Kokkos::Device<execution_space, memory_space> device_type;

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
  typedef Kokkos::StaticCrsGraph<OrdinalType, Kokkos::LayoutLeft, execution_space, SizeType> StaticCrsGraphType;
  //! Type of column indices in the sparse matrix.
  typedef typename StaticCrsGraphType::entries_type index_type;
  //! Const version of the type of column indices in the sparse matrix.
  typedef typename index_type::non_const_value_type const_ordinal_type;
  //! Nonconst version of the type of column indices in the sparse matrix.
  typedef typename index_type::non_const_value_type non_const_ordinal_type;
  //! Type of the "row map" (which contains the offset for each row's data).
  typedef typename StaticCrsGraphType::row_map_type row_map_type;
  //! Const version of the type of row offsets in the sparse matrix.
  typedef typename row_map_type::non_const_value_type const_size_type;
  //! Nonconst version of the type of row offsets in the sparse matrix.
  typedef typename row_map_type::non_const_value_type non_const_size_type;
  //! Kokkos Array type of the entries (values) in the sparse matrix.
  typedef Kokkos::View<value_type*, Kokkos::LayoutRight, device_type, MemoryTraits> values_type;
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
    dev_config (B.dev_config),
#ifdef KOKKOS_USE_CUSPARSE
    cusparse_handle (B.cusparse_handle),
    cusparse_descr (B.cusparse_descr),
#endif // KOKKOS_USE_CUSPARSE
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
  template<typename SType>
  KOKKOS_INLINE_FUNCTION
  SparseRowView<CrsMatrix,SType> row (const ordinal_type i) const {
    const size_type start = graph.row_map(i);
    const size_type count = graph.row_map(i+1) - start;

    if (count == 0) {
      return SparseRowView<CrsMatrix,SType> (NULL, NULL, 1, 0);
    } else {
      return SparseRowView<CrsMatrix,SType> (values, graph.entries, 1, count, start);
    }
  }

  /// \brief Return a const view of row i of the matrix.
  ///
  /// If row i does not belong to the matrix, return an empty view.
  template<typename SType>
  KOKKOS_INLINE_FUNCTION
  SparseRowViewConst<CrsMatrix,SType> rowConst (const ordinal_type i) const {
    const size_type start = graph.row_map(i);
    const size_type count = graph.row_map(i+1) - start;

    if (count == 0) {
      return SparseRowViewConst<CrsMatrix,SType> (NULL, NULL, 1, 0);
    } else {
      return SparseRowViewConst<CrsMatrix,SType> (values, graph.entries, 1, count, start);
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
}
#endif
