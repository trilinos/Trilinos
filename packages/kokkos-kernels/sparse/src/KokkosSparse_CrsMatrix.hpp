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

/// \file KokkosSparse_CrsMatrix.hpp
/// \brief Local sparse matrix interface
///
/// This file provides KokkosSparse::CrsMatrix.  This implements a
/// local (no MPI) sparse matrix stored in compressed row sparse
/// ("Crs") format.

#ifndef KOKKOS_SPARSE_CRSMATRIX_HPP_
#define KOKKOS_SPARSE_CRSMATRIX_HPP_

#include "Kokkos_Core.hpp"
#include "Kokkos_StaticCrsGraph.hpp"
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include "KokkosSparse_findRelOffset.hpp"
#include "KokkosKernels_default_types.hpp"
#include "KokkosKernels_Macros.hpp"

namespace KokkosSparse {
//! String that tells sparse kernels to use the transpose of the matrix.
static char KOKKOSKERNELS_UNUSED_ATTRIBUTE Transpose[] = "T";
/// \brief String that tells sparse kernels to use the conjugate (NOT
///   transpose) of the matrix.
static char KOKKOSKERNELS_UNUSED_ATTRIBUTE Conjugate[] = "C";
/// \brief String that tells sparse kernels to use the conjugate
///   transpose of the matrix.
static char KOKKOSKERNELS_UNUSED_ATTRIBUTE ConjugateTranspose[] = "H";
/// \brief String that tells sparse kernels not to use the transpose
///   or conjugate of the matrix.
static char KOKKOSKERNELS_UNUSED_ATTRIBUTE NoTranspose[] = "N";

template <class DeviceType>
inline int RowsPerThread(const int NNZPerRow) {
  if (NNZPerRow == 0) return 1;
  int result = 2;
  while (result * NNZPerRow <= 2048) {
    result *= 2;
  }
  return result / 2;
}
#ifdef KOKKOS_ENABLE_CUDA
template <>
inline int RowsPerThread<Kokkos::Cuda>(const int /*NNZPerRow*/) {
  return 1;
}
#endif
#ifdef KOKKOS_ENABLE_HIP
template <>
inline int RowsPerThread<Kokkos::HIP>(const int /*NNZPerRow*/) {
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
    KOKKOS_INLINE_FUNCTION
    Dim3(const size_t x_, const size_t y_ = 1, const size_t z_ = 1) : x(x_), y(y_), z(z_) {}
  };

  Dim3 block_dim;
  size_t num_blocks;
  size_t num_threads_per_block;

  KOKKOS_INLINE_FUNCTION
  DeviceConfig(const size_t num_blocks_ = 0, const size_t threads_per_block_x_ = 0,
               const size_t threads_per_block_y_ = 0, const size_t threads_per_block_z_ = 1)
      : block_dim(threads_per_block_x_, threads_per_block_y_, threads_per_block_z_),
        num_blocks(num_blocks_),
        num_threads_per_block(block_dim.x * block_dim.y * block_dim.z) {}
};

/// \class SparseRowView
/// \brief View of a row of a sparse matrix.
/// \tparam MatrixType Sparse matrix type, such as (but not limited to)
/// CrsMatrix.
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
/// const ordinal_type numEntries = A_i.length;
/// for (ordinal_type k = 0; k < numEntries; ++k) {
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
template <class MatrixType>
struct SparseRowView {
  //! The type of the values in the row.
  typedef typename MatrixType::value_type value_type;
  //! The type of the column indices in the row.
  typedef typename MatrixType::ordinal_type ordinal_type;

 private:
  //! Array of values in the row.
  value_type* values_;
  //! Array of (local) column indices in the row.
  ordinal_type* colidx_;
  /// \brief Stride between successive entries in the row.
  ///
  /// For compressed sparse row (CSR) storage, this is always one.
  /// This might be greater than one for storage formats like ELLPACK
  /// or Jagged Diagonal.  Nevertheless, the stride can never be
  /// greater than the number of rows or columns in the matrix.  Thus,
  /// \c ordinal_type is the correct type.
  const ordinal_type stride_;

 public:
  /// \brief Constructor
  ///
  /// \param values   [in] Array of the row's values.
  /// \param colidx__ [in] Array of the row's column indices.
  /// \param stride   [in] (Constant) stride between matrix entries in
  ///   each of the above arrays.
  /// \param count    [in] Number of entries in the row.
  KOKKOS_INLINE_FUNCTION
  SparseRowView(value_type* const values, ordinal_type* const colidx__, const ordinal_type& stride,
                const ordinal_type& count)
      : values_(values), colidx_(colidx__), stride_(stride), length(count) {}

  /// \brief Constructor with offset into \c colidx array
  /// \tparam OffsetType The type of \c idx (see above).  Must be a
  ///   built-in integer type.  This may differ from ordinal_type.
  ///   For example, the matrix may have dimensions that fit in int,
  ///   but a number of entries that does not fit in int.
  ///
  /// \param values   [in] Array of the row's values.
  /// \param colidx__ [in] Array of the row's column indices.
  /// \param stride   [in] (Constant) stride between matrix entries in
  ///                 each of the above arrays.
  /// \param count    [in] Number of entries in the row.
  /// \param idx      [in] Start offset into \c colidx array
  template <class OffsetType>
  KOKKOS_INLINE_FUNCTION SparseRowView(
      const typename MatrixType::values_type& values, const typename MatrixType::index_type& colidx__,
      const ordinal_type& stride, const ordinal_type& count, const OffsetType& idx,
      const typename std::enable_if<std::is_integral<OffsetType>::value, int>::type& = 0)
      : values_(&values(idx)), colidx_(&colidx__(idx)), stride_(stride), length(count) {}

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
  value_type& value(const ordinal_type& i) const { return values_[i * stride_]; }

  /// \brief Reference to the column index of entry i in this row of the sparse
  /// matrix.
  ///
  /// "Entry i" is not necessarily the entry with column index i, nor
  /// does i necessarily correspond to the (local) row index.
  KOKKOS_INLINE_FUNCTION
  ordinal_type& colidx(const ordinal_type& i) const { return colidx_[i * stride_]; }
};

/// \class SparseRowViewConst
/// \brief Const view of a row of a sparse matrix.
/// \tparam MatrixType Sparse matrix type, such as (but not limited to)
/// CrsMatrix.
///
/// This class is like SparseRowView, except that it provides a const
/// view.  This class exists in order to let users get a const view of
/// a row of a nonconst matrix.
template <class MatrixType>
struct SparseRowViewConst {
  //! The type of the values in the row.
  typedef const typename MatrixType::non_const_value_type value_type;
  //! The type of the column indices in the row.
  typedef const typename MatrixType::non_const_ordinal_type ordinal_type;

 private:
  //! Array of values in the row.
  value_type* values_;
  //! Array of (local) column indices in the row.
  ordinal_type* colidx_;
  /// \brief Stride between successive entries in the row.
  ///
  /// For compressed sparse row (CSR) storage, this is always one.
  /// This might be greater than one for storage formats like ELLPACK
  /// or Jagged Diagonal.  Nevertheless, the stride can never be
  /// greater than the number of rows or columns in the matrix.  Thus,
  /// \c ordinal_type is the correct type.
  const ordinal_type stride_;

 public:
  /// \brief Constructor
  ///
  /// \param values   [in] Array of the row's values.
  /// \param colidx__ [in] Array of the row's column indices.
  /// \param stride   [in] (Constant) stride between matrix entries in
  ///                 each of the above arrays.
  /// \param count    [in] Number of entries in the row.
  KOKKOS_INLINE_FUNCTION
  SparseRowViewConst(value_type* const values, ordinal_type* const colidx__, const ordinal_type& stride,
                     const ordinal_type& count)
      : values_(values), colidx_(colidx__), stride_(stride), length(count) {}

  /// \brief Constructor with offset into \c colidx array
  ///
  /// \tparam OffsetType The type of \c idx (see above).  Must be a
  ///   built-in integer type.  This may differ from ordinal_type.
  ///   For example, the matrix may have dimensions that fit in int,
  ///   but a number of entries that does not fit in int.
  /// \param values   [in] Array of the row's values.
  /// \param colidx__ [in] Array of the row's column indices.
  /// \param stride   [in] (Constant) stride between matrix entries in
  ///                 each of the above arrays.
  /// \param count    [in] Number of entries in the row.
  /// \param idx      [in] Start offset into \c colidx array
  template <class OffsetType>
  KOKKOS_INLINE_FUNCTION SparseRowViewConst(
      const typename MatrixType::values_type& values, const typename MatrixType::index_type& colidx__,
      const ordinal_type& stride, const ordinal_type& count, const OffsetType& idx,
      const typename std::enable_if<std::is_integral<OffsetType>::value, int>::type& = 0)
      : values_(&values(idx)), colidx_(&colidx__(idx)), stride_(stride), length(count) {}

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
  value_type& value(const ordinal_type& i) const { return values_[i * stride_]; }

  /// \brief (Const) reference to the column index of entry i in this
  ///   row of the sparse matrix.
  ///
  /// "Entry i" is not necessarily the entry with column index i, nor
  /// does i necessarily correspond to the (local) row index.
  KOKKOS_INLINE_FUNCTION
  ordinal_type& colidx(const ordinal_type& i) const { return colidx_[i * stride_]; }
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
template <class ScalarType, class OrdinalType, class Device, class MemoryTraits = void,
          class SizeType = default_size_type>
class CrsMatrix {
  static_assert(std::is_signed<OrdinalType>::value, "CrsMatrix requires that OrdinalType is a signed integer type.");

 private:
  typedef typename Kokkos::ViewTraits<ScalarType*, Device, void, MemoryTraits>::host_mirror_space host_mirror_space;

 public:
  //! Type of the matrix's execution space.
  typedef typename Device::execution_space execution_space;
  //! Type of the matrix's memory space.
  typedef typename Device::memory_space memory_space;
  //! Canonical device type
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
  typedef CrsMatrix<ScalarType, OrdinalType, host_mirror_space, MemoryTraits, SizeType> HostMirror;
  //! Type of the graph structure of the sparse matrix.
  typedef Kokkos::StaticCrsGraph<ordinal_type, default_layout, device_type, memory_traits, size_type>
      StaticCrsGraphType;
  //! Type of the graph structure of the sparse matrix - consistent with Kokkos.
  typedef Kokkos::StaticCrsGraph<ordinal_type, default_layout, device_type, memory_traits, size_type>
      staticcrsgraph_type;
  //! Type of column indices in the sparse matrix.
  typedef typename staticcrsgraph_type::entries_type index_type;
  //! Const version of the type of column indices in the sparse matrix.
  typedef typename index_type::const_value_type const_ordinal_type;
  //! Nonconst version of the type of column indices in the sparse matrix.
  typedef typename index_type::non_const_value_type non_const_ordinal_type;
  //! Type of the "row map" (which contains the offset for each row's data).
  typedef typename staticcrsgraph_type::row_map_type row_map_type;
  //! Const version of the type of row offsets in the sparse matrix.
  typedef typename row_map_type::const_value_type const_size_type;
  //! Nonconst version of the type of row offsets in the sparse matrix.
  typedef typename row_map_type::non_const_value_type non_const_size_type;
  //! Kokkos Array type of the entries (values) in the sparse matrix.
  typedef Kokkos::View<value_type*, Kokkos::LayoutRight, device_type, MemoryTraits> values_type;
  //! Const version of the type of the entries in the sparse matrix.
  typedef typename values_type::const_value_type const_value_type;
  //! Nonconst version of the type of the entries in the sparse matrix.
  typedef typename values_type::non_const_value_type non_const_value_type;

  typedef CrsMatrix<const_value_type, ordinal_type, device_type, memory_traits, size_type> const_type;

  /// \name Storage of the actual sparsity structure and values.
  ///
  /// CrsMatrix uses the compressed sparse row (CSR) storage format to
  /// store the sparse matrix.  CSR is also called "compressed row
  /// storage"; hence the name, which it inherits from Tpetra and from
  /// Epetra before it.
  //@{
  //! The graph (sparsity structure) of the sparse matrix.
  staticcrsgraph_type graph;
  //! The 1-D array of values of the sparse matrix.
  values_type values;
  //@}

 private:
  /// \brief The number of distinct column indices used by the matrix
  ///
  /// This value might not be exact but rather an upper bound of the
  /// number of distinct column indices used by the matrix.
  /// It provides multiple sparse algorithms to allocate appropriate
  /// amount of temporary work space or to allocate memory for the
  /// output of the kernel.
  ordinal_type numCols_;

 public:
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
  KOKKOS_INLINE_FUNCTION
  CrsMatrix() : numCols_(0) {}

  //! Copy constructor (shallow copy).
  template <typename InScalar, typename InOrdinal, class InDevice, class InMemTraits, typename InSizeType>
  KOKKOS_INLINE_FUNCTION CrsMatrix(const CrsMatrix<InScalar, InOrdinal, InDevice, InMemTraits, InSizeType>& B)
      : graph(B.graph.entries, B.graph.row_map), values(B.values), numCols_(B.numCols()), dev_config(B.dev_config) {
    graph.row_block_offsets = B.graph.row_block_offsets;
    // TODO: MD 07/2017: Changed the copy constructor of graph
    // as the constructor of StaticCrsGraph does not allow copy from non const
    // version.
  }

  //! Deep copy constructor (can cross spaces)
  template <typename InScalar, typename InOrdinal, typename InDevice, typename InMemTraits, typename InSizeType>
  CrsMatrix(const std::string&, const CrsMatrix<InScalar, InOrdinal, InDevice, InMemTraits, InSizeType>& mat_) {
    typename row_map_type::non_const_type rowmap(Kokkos::view_alloc(Kokkos::WithoutInitializing, "rowmap"),
                                                 mat_.graph.row_map.extent(0));
    index_type cols(Kokkos::view_alloc(Kokkos::WithoutInitializing, "cols"), mat_.nnz());
    values = values_type(Kokkos::view_alloc(Kokkos::WithoutInitializing, "values"), mat_.nnz());
    Kokkos::deep_copy(rowmap, mat_.graph.row_map);
    Kokkos::deep_copy(cols, mat_.graph.entries);
    Kokkos::deep_copy(values, mat_.values);

    numCols_ = mat_.numCols();
    graph    = StaticCrsGraphType(cols, rowmap);
  }

  /// \brief Construct with a graph that will be shared.
  ///
  /// Allocate the values array for subsquent fill.
  template <typename InOrdinal, typename InLayout, typename InDevice, typename InMemTraits, typename InSizeType>
  [[deprecated(
      "Use the constructor that accepts ncols as input "
      "instead.")]] CrsMatrix(const std::string& label,
                              const Kokkos::StaticCrsGraph<InOrdinal, InLayout, InDevice, InMemTraits, InSizeType>&
                                  graph_)
      : graph(graph_.entries, graph_.row_map),
        values(label, graph_.entries.extent(0)),
        numCols_(maximum_entry(graph_) + 1) {}

  /// \brief Constructor that accepts a a static graph, and numCols.
  ///
  /// The matrix will store and use the row map, indices
  /// (by view, not by deep copy) and allocate the values view.
  ///
  /// \param label  [in] The sparse matrix's label.
  /// \param graph_ [in] The graph for storing the rowmap and col ids.
  /// \param ncols  [in] The number of columns.
  template <typename InOrdinal, typename InLayout, typename InDevice, typename InMemTraits, typename InSizeType>
  CrsMatrix(const std::string& label,
            const Kokkos::StaticCrsGraph<InOrdinal, InLayout, InDevice, InMemTraits, InSizeType>& graph_,
            const OrdinalType& ncols)
      : graph(graph_.entries, graph_.row_map), values(label, graph_.entries.extent(0)), numCols_(ncols) {}

  /// \brief Constructor that accepts a a static graph, and values.
  ///
  /// The matrix will store and use the row map, indices, and values
  /// directly (by view, not by deep copy).
  ///
  /// \param ncols [in] The number of columns.
  /// \param vals [in/out] The entries.
  /// \param graph_ The graph for storing the rowmap and col ids.
  template <typename InOrdinal, typename InLayout, typename InDevice, typename InMemTraits, typename InSizeType>
  CrsMatrix(const std::string&, const OrdinalType& ncols, const values_type& vals,
            const Kokkos::StaticCrsGraph<InOrdinal, InLayout, InDevice, InMemTraits, InSizeType>& graph_)
      : graph(graph_.entries, graph_.row_map), values(vals), numCols_(ncols) {}

  /// \brief Constructor that copies raw arrays of host data in
  ///   3-array CRS (compresed row storage) format.
  ///
  /// On input, the entries must be sorted by row. \c rowmap determines where
  /// each row begins and ends. For each entry k (0 <= k < annz), \c cols[k]
  /// gives the adjacent column, and \c val[k] gives the corresponding matrix
  /// value.
  ///
  /// This constructor is mainly useful for benchmarking or for
  /// reading the sparse matrix's data from a file.
  ///
  /// \param nrows [in] The number of rows.
  /// \param ncols [in] The number of columns.
  /// \param annz [in] The number of entries.
  /// \param val [in] The values.
  /// \param rowmap [in] The row offsets. The values/columns in row k begin at
  /// index
  ///   \c rowmap[k] and end at \c rowmap[k+1]-1 (inclusive). This means the
  ///   array must have length \c nrows+1.
  /// \param cols [in] The column indices. \c cols[k] is the column
  ///   index of entry k, with a corresponding value of \c val[k] .
  CrsMatrix(const std::string& /*label*/, OrdinalType nrows, OrdinalType ncols, size_type annz, ScalarType* val,
            OrdinalType* rowmap, OrdinalType* cols) {
    using Kokkos::Unmanaged;
    using HostRowmap       = Kokkos::View<SizeType*, Kokkos::HostSpace>;
    using UnmanagedRowmap  = Kokkos::View<const SizeType*, Kokkos::HostSpace, Kokkos::MemoryTraits<Unmanaged>>;
    using UnmanagedEntries = Kokkos::View<const OrdinalType*, Kokkos::HostSpace, Kokkos::MemoryTraits<Unmanaged>>;
    using UnmanagedValues  = Kokkos::View<const ScalarType*, Kokkos::HostSpace, Kokkos::MemoryTraits<Unmanaged>>;
    // Allocate device rowmap, entries, values views
    typename row_map_type::non_const_type rowmapDevice(Kokkos::view_alloc(Kokkos::WithoutInitializing, "rowmap"),
                                                       nrows + 1);
    index_type entriesDevice(Kokkos::view_alloc(Kokkos::WithoutInitializing, "entries"), annz);
    // given rowmap in ordinal_type, so may need to convert to size_type
    // explicitly
    HostRowmap rowmapConverted;
    UnmanagedRowmap rowmapRaw;
    if (!std::is_same<OrdinalType, SizeType>::value) {
      rowmapConverted = HostRowmap(Kokkos::view_alloc(Kokkos::WithoutInitializing, "rowmap raw"), nrows + 1);
      for (OrdinalType i = 0; i <= nrows; i++) rowmapConverted(i) = rowmap[i];
      rowmapRaw = rowmapConverted;
    } else {
      rowmapRaw = UnmanagedRowmap((const SizeType*)rowmap, nrows + 1);
    }
    Kokkos::deep_copy(rowmapDevice, rowmapRaw);
    UnmanagedEntries entriesRaw(cols, annz);
    Kokkos::deep_copy(entriesDevice, entriesRaw);
    // Construct graph and populate all members
    this->numCols_ = ncols;
    this->graph    = StaticCrsGraphType(entriesDevice, rowmapDevice);
    this->values   = values_type(Kokkos::view_alloc(Kokkos::WithoutInitializing, "values"), annz);
    UnmanagedValues valuesRaw(val, annz);
    Kokkos::deep_copy(this->values, valuesRaw);
  }

  /// \brief Constructor that accepts a row map, column indices, and
  ///   values.
  ///
  /// The matrix will store and use the row map, indices, and values
  /// directly (by view, not by deep copy).
  ///
  /// \param nrows [in] The number of rows.
  /// \param ncols [in] The number of columns.
  /// \param annz [in] The number of entries.
  /// \param vals [in] The entries.
  /// \param rowmap [in] The row map (containing the offsets to the
  ///   data in each row).
  /// \param cols [in] The column indices.
  CrsMatrix(const std::string& /* label */, const OrdinalType nrows, const OrdinalType ncols, const size_type annz,
            const values_type& vals, const row_map_type& rowmap, const index_type& cols)
      : graph(cols, rowmap), values(vals), numCols_(ncols) {
    const ordinal_type actualNumRows = (rowmap.extent(0) != 0)
                                           ? static_cast<ordinal_type>(rowmap.extent(0) - static_cast<size_type>(1))
                                           : static_cast<ordinal_type>(0);
    if (nrows != actualNumRows) {
      std::ostringstream os;
      os << "Input argument nrows = " << nrows
         << " != the actual number of "
            "rows "
         << actualNumRows << " according to the 'rows' input argument.";
      throw std::invalid_argument(os.str());
    }
    if (annz != nnz()) {
      std::ostringstream os;
      os << "Input argument annz = " << annz << " != this->nnz () = " << nnz() << ".";
      throw std::invalid_argument(os.str());
    }
  }

  KOKKOS_INLINE_FUNCTION
  OrdinalType sumIntoValues(const OrdinalType rowi, const OrdinalType cols[], const OrdinalType ncol,
                            const ScalarType vals[], const bool is_sorted = false,
                            const bool force_atomic = false) const {
    SparseRowView<CrsMatrix> row_view = this->row(rowi);
    const ordinal_type length         = row_view.length;

    ordinal_type hint     = 0;  // Guess for offset of current column index in row
    ordinal_type numValid = 0;  // number of valid local column indices

    for (ordinal_type i = 0; i < ncol; ++i) {
      // NOTE (mfh 19 Sep 2017) This assumes that row_view stores
      // column indices contiguously.  It does, but one could imagine
      // changing that at some point.
      const ordinal_type offset = findRelOffset(&(row_view.colidx(0)), length, cols[i], hint, is_sorted);
      if (offset != length) {
        if (force_atomic) {
          Kokkos::atomic_add(&(row_view.value(offset)), vals[i]);
        } else {
          row_view.value(offset) += vals[i];
        }
        ++numValid;
        // If the hint is out of range, findRelOffset will ignore it.
        // Thus, while it's harmless to have a hint out of range, it
        // may slow down searches for subsequent valid input column
        // indices.
        hint = offset + 1;
      }
    }
    return numValid;
  }

  KOKKOS_INLINE_FUNCTION
  OrdinalType replaceValues(const OrdinalType rowi, const OrdinalType cols[], const OrdinalType ncol,
                            const ScalarType vals[], const bool is_sorted = false,
                            const bool force_atomic = false) const {
    SparseRowView<CrsMatrix> row_view = this->row(rowi);
    const ordinal_type length         = row_view.length;

    ordinal_type hint     = 0;  // Guess for offset of current column index in row
    ordinal_type numValid = 0;  // number of valid local column indices

    for (ordinal_type i = 0; i < ncol; ++i) {
      // NOTE (mfh 19 Sep 2017) This assumes that row_view stores
      // column indices contiguously.  It does, but one could imagine
      // changing that at some point.
      const ordinal_type offset = findRelOffset(&(row_view.colidx(0)), length, cols[i], hint, is_sorted);
      if (offset != length) {
        if (force_atomic) {
          Kokkos::atomic_assign(&(row_view.value(offset)), vals[i]);
        } else {
          row_view.value(offset) = vals[i];
        }
        ++numValid;
        // If the hint is out of range, findRelOffset will ignore it.
        // Thus, while it's harmless to have a hint out of range, it
        // may slow down searches for subsequent valid input column
        // indices.
        hint = offset + 1;
      }
    }
    return numValid;
  }

  //! Attempt to assign the input matrix to \c *this.
  template <typename aScalarType, typename aOrdinalType, class aDevice, class aMemoryTraits, typename aSizeType>
  CrsMatrix& operator=(const CrsMatrix<aScalarType, aOrdinalType, aDevice, aMemoryTraits, aSizeType>& mtx) {
    numCols_   = mtx.numCols();
    graph      = mtx.graph;
    values     = mtx.values;
    dev_config = mtx.dev_config;
    return *this;
  }

  //! The number of rows in the sparse matrix.
  KOKKOS_INLINE_FUNCTION ordinal_type numRows() const { return graph.numRows(); }

  //! The number of columns in the sparse matrix.
  KOKKOS_INLINE_FUNCTION ordinal_type numCols() const { return numCols_; }

  //! The number of "point" (non-block) rows in the matrix. Since Crs is not
  //! blocked, this is just the number of regular rows.
  KOKKOS_INLINE_FUNCTION ordinal_type numPointRows() const { return numRows(); }

  //! The number of "point" (non-block) columns in the matrix. Since Crs is not
  //! blocked, this is just the number of regular columns.
  KOKKOS_INLINE_FUNCTION ordinal_type numPointCols() const { return numCols(); }

  //! The number of stored entries in the sparse matrix.
  KOKKOS_INLINE_FUNCTION size_type nnz() const { return graph.entries.extent(0); }

  friend struct SparseRowView<CrsMatrix>;

  /// \brief Return a view of row i of the matrix.
  ///
  /// If row i does not belong to the matrix, return an empty view.
  ///
  /// The returned object \c view implements the following interface:
  /// <ul>
  /// <li> \c view.length is the number of entries in the row </li>
  /// <li> \c view.value(k) returns a nonconst reference
  ///      to the value of the k-th entry in the row </li>
  /// <li> \c view.colidx(k) returns a nonconst reference to
  ///      the column index of the k-th entry in the row </li>
  /// </ul>
  /// k is not a column index; it just counts from 0 to
  /// <tt>view.length - 1</tt>.
  ///
  /// Users should not rely on the return type of this method.  They
  /// should instead assign to 'auto'.  That allows compile-time
  /// polymorphism for different kinds of sparse matrix formats (e.g.,
  /// ELLPACK or Jagged Diagonal) that we may wish to support in the
  /// future.
  ///
  /// Both row() and rowConst() used to take a "SizeType" template
  /// parameter, which was the type to use for row offsets.  This is
  /// unnecessary, because the CrsMatrix specialization already has
  /// the row offset type available, via the <tt>size_type</tt>
  /// typedef.  Our sparse matrix-vector multiply implementation for
  /// CrsMatrix safely uses <tt>ordinal_type</tt> rather than
  /// <tt>size_type</tt> to iterate over all the entries in a row of
  /// the sparse matrix.  Since <tt>ordinal_type</tt> may be smaller
  /// than <tt>size_type</tt>, compilers may generate more efficient
  /// code.  The row() and rowConst() methods first compute the
  /// difference of consecutive row offsets as <tt>size_type</tt>, and
  /// then cast to <tt>ordinal_type</tt>.  If you want to do this
  /// yourself, here is an example:
  ///
  /// \code
  /// for (ordinal_type lclRow = 0; lclRow < A.numRows (); ++lclRow) {
  ///   const ordinal_type numEnt =
  ///     static_cast<ordinal_type> (A.graph.row_map(i+1) - A.graph.row_map(i));
  ///   for (ordinal_type k = 0; k < numEnt; ++k) {
  ///     // etc.
  ///   }
  /// }
  /// \endcode
  KOKKOS_INLINE_FUNCTION
  SparseRowView<CrsMatrix> row(const ordinal_type i) const {
    const size_type start = graph.row_map(i);
    // count is guaranteed to fit in ordinal_type, as long as no row
    // has duplicate entries.
    const ordinal_type count = static_cast<ordinal_type>(graph.row_map(i + 1) - start);

    if (count == 0) {
      return SparseRowView<CrsMatrix>(NULL, NULL, 1, 0);
    } else {
      return SparseRowView<CrsMatrix>(values, graph.entries, 1, count, start);
    }
  }

  /// \brief Return a const view of row i of the matrix.
  ///
  /// If row i does not belong to the matrix, return an empty view.
  ///
  /// The returned object \c view implements the following interface:
  /// <ul>
  /// <li> \c view.length is the number of entries in the row </li>
  /// <li> \c view.value(k) returns a const reference to
  ///      the value of the k-th entry in the row </li>
  /// <li> \c view.colidx(k) returns a const reference to the
  ///      column index of the k-th entry in the row </li>
  /// </ul>
  /// k is not a column index; it just counts from 0 to
  /// <tt>view.length - 1</tt>.
  ///
  /// Users should not rely on the return type of this method.  They
  /// should instead assign to 'auto'.  That allows compile-time
  /// polymorphism for different kinds of sparse matrix formats (e.g.,
  /// ELLPACK or Jagged Diagonal) that we may wish to support in the
  /// future.
  ///
  /// Both row() and rowConst() used to take a "SizeType" template
  /// parameter, which was the type to use for row offsets.  This is
  /// unnecessary, because the CrsMatrix specialization already has
  /// the row offset type available, via the <tt>size_type</tt>
  /// typedef.  Our sparse matrix-vector multiply implementation for
  /// CrsMatrix safely uses <tt>ordinal_type</tt> rather than
  /// <tt>size_type</tt> to iterate over all the entries in a row of
  /// the sparse matrix.  Since <tt>ordinal_type</tt> may be smaller
  /// than <tt>size_type</tt>, compilers may generate more efficient
  /// code.  The row() and rowConst() methods first compute the
  /// difference of consecutive row offsets as <tt>size_type</tt>, and
  /// then cast to <tt>ordinal_type</tt>.  If you want to do this
  /// yourself, here is an example:
  ///
  /// \code
  /// for (ordinal_type lclRow = 0; lclRow < A.numRows (); ++lclRow) {
  ///   const ordinal_type numEnt =
  ///     static_cast<ordinal_type> (A.graph.row_map(i+1) - A.graph.row_map(i));
  ///   for (ordinal_type k = 0; k < numEnt; ++k) {
  ///     // etc.
  ///   }
  /// }
  /// \endcode
  KOKKOS_INLINE_FUNCTION
  SparseRowViewConst<CrsMatrix> rowConst(const ordinal_type i) const {
    const size_type start = graph.row_map(i);
    // count is guaranteed to fit in ordinal_type, as long as no row
    // has duplicate entries.
    const ordinal_type count = static_cast<ordinal_type>(graph.row_map(i + 1) - start);

    if (count == 0) {
      return SparseRowViewConst<CrsMatrix>(NULL, NULL, 1, 0);
    } else {
      return SparseRowViewConst<CrsMatrix>(values, graph.entries, 1, count, start);
    }
  }
};

/// \class is_crs_matrix
/// \brief is_crs_matrix<T>::value is true if T is a CrsMatrix<...>, false
/// otherwise
template <typename>
struct is_crs_matrix : public std::false_type {};
template <typename... P>
struct is_crs_matrix<CrsMatrix<P...>> : public std::true_type {};
template <typename... P>
struct is_crs_matrix<const CrsMatrix<P...>> : public std::true_type {};

/// \brief Equivalent to is_crs_matrix<T>::value.
template <typename T>
inline constexpr bool is_crs_matrix_v = is_crs_matrix<T>::value;

}  // namespace KokkosSparse
#endif
