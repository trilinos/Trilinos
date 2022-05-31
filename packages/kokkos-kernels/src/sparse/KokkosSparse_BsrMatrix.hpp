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
// Questions? Contact Luc Berger-Vergiat (lberge@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

/// \file Kokkos_Sparse_BsrMatrix.hpp
/// \brief Local sparse matrix interface
///
/// This file provides KokkosSparse::Experimental::BsrMatrix.
/// This implements a local (no MPI) sparse matrix stored in block-by-block
/// compressed row sparse format.

#ifndef KOKKOS_SPARSE_BSRMATRIX_HPP_
#define KOKKOS_SPARSE_BSRMATRIX_HPP_

#include <set>
#include <sstream>
#include <stdexcept>
#include <type_traits>

#include "Kokkos_Core.hpp"
#include "Kokkos_StaticCrsGraph.hpp"
#include "Kokkos_ArithTraits.hpp"
#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosKernels_Error.hpp"

namespace KokkosSparse {

namespace Experimental {

template <class MatrixType>
struct BsrRowView {
  //! The type of the values in the row.
  typedef typename MatrixType::value_type value_type;
  //! The type of the column indices in the row.
  typedef typename MatrixType::ordinal_type ordinal_type;
  //! The type for returned block of values.
  typedef Kokkos::View<value_type**, Kokkos::LayoutRight,
                       typename MatrixType::device_type,
                       Kokkos::MemoryUnmanaged>
      block_values_type;

 private:
  //! Array of values in the row.
  value_type* values_;
  //! Array of (local) column indices in the row.
  ordinal_type* colidx_;
  /// \brief Stride between successive rows in a block.
  ///
  const ordinal_type blockDim_;

 public:
  /// \brief Constructor
  ///
  /// \param values [in] Array of the row's values.
  /// \param colidx [in] Array of the row's column indices.
  /// \param blockDim [in] (Constant) stride between block rows
  ///   within a block-row in the above arrays.
  /// \param count [in] Number of blocks in the desired block-row.
  //
  // Assumes values and colidx already offset to the correct location
  KOKKOS_INLINE_FUNCTION
  BsrRowView(value_type* const values, ordinal_type* const colidx,
             const ordinal_type& blockDim, const ordinal_type& count)
      : values_(values), colidx_(colidx), blockDim_(blockDim), length(count) {}

  /// \brief Constructor with offset into \c colidx array
  ///
  /// \param values [in] Array of the row's values.
  /// \param colidx [in] Array of the row's column indices.
  /// \param blockDim [in] (Constant) stride between rows in
  ///   within a block in the above arrays.
  /// \param count [in] Number of blocks in the desired block-row
  /// \param start [in] Offset into values and colidx of the desired block-row
  /// start.
  ///   Note: The offset into the values array for a block-row equals
  ///           num_blocks_prior_to_block-row*blockDim*blockDim
  ///
  /// \tparam OffsetType The type of \c start (see above).  Must be a
  ///   built-in integer type.  This may differ from ordinal_type.
  ///   For example, the matrix may have dimensions that fit in int,
  ///   but a number of entries that does not fit in int.
  template <class OffsetType>
  KOKKOS_INLINE_FUNCTION BsrRowView(
      const typename MatrixType::values_type& values,
      const typename MatrixType::index_type& colidx,
      const ordinal_type& blockDim, const ordinal_type& count,
      const OffsetType& start,
      const typename std::enable_if<std::is_integral<OffsetType>::value,
                                    int>::type& = 0)
      : values_(&values(start * blockDim * blockDim)),
        colidx_(&colidx(start)),
        blockDim_(blockDim),
        length(count) {}

  /// \brief Number of entries (i.e. blocks) in the row.
  ///
  /// This is a public const field rather than a public const method,
  /// in order to avoid possible overhead of a method call if the
  /// compiler is unable to inline that method call.
  ///
  /// We assume that rows contain no duplicate entries (i.e., entries
  /// with the same column index).  Thus, a row may have up to
  /// A.numCols() entries.  This means that the correct type of
  /// 'length' is ordinal_type.
  /// Here, length refers to the number of blocks in a block-row
  const ordinal_type length;

  /// /brief Return a pointer offset to local row i of block K of values_ array;
  ///        user responsible for indexing into this pointer correctly
  /// \param K [in] must be the LOCAL block index within this block-row
  /// \param i [in] must be the LOCAL row index offset within this block-row
  ///
  /// Output: pointer to values_ array at start of local row within block K
  ///
  /// Pointer interfaces are NOT guaranteed for backward compatibility
  /// This interface is intended for performant kernels, not common usage
  KOKKOS_INLINE_FUNCTION
  value_type* local_row_in_block(const ordinal_type& K,
                                 const ordinal_type& i) const {
    return (values_ + (K * blockDim_ * blockDim_ + i * blockDim_));
  }

  /// \brief Return the value at a specified block K of block-row
  ///        with local row and col offset (i,j)
  /// \param K [in] must be the LOCAL block index within this block-row
  /// \param i [in] must be the LOCAL row index offset within this block-row
  /// \param j [in] must be the LOCAL col index offset within this block-row
  ///
  /// Output: reference to value_type at the given (K, i, j) offset into values_
  KOKKOS_INLINE_FUNCTION
  value_type& local_block_value(const ordinal_type& K, const ordinal_type& i,
                                const ordinal_type& j) const {
    return values_[K * blockDim_ * blockDim_ + i * blockDim_ + j];
  }

  /// \brief Return unmanaged 2D strided View wrapping local block K from this
  /// block-row
  /// \param K [in] must be the LOCAL block index within this
  /// block-row
  KOKKOS_INLINE_FUNCTION
  block_values_type block(const ordinal_type& K) const {
    return block_values_type(&(values_[K * blockDim_ * blockDim_]),
                             Kokkos::LayoutRight(blockDim_, blockDim_));
  }

  /// \brief Return offset into colidx_ for the requested block idx
  ///        If none found, return Kokkos::Details::ArithTraits::max
  /// \param idx_to_match [in] local block idx within block-row
  /// \param is_sorted [in] defaulted to false; no usage at this time
  KOKKOS_INLINE_FUNCTION
  ordinal_type findRelBlockOffset(const ordinal_type idx_to_match,
                                  bool /*is_sorted*/ = false) const {
    ordinal_type offset = Kokkos::Details::ArithTraits<ordinal_type>::max();
    for (ordinal_type blk_offset = 0; blk_offset < length; ++blk_offset) {
      ordinal_type idx = colidx_[blk_offset];
      if (idx == idx_to_match) {
        offset = blk_offset;
        break;
      }  // return relative offset
    }
    return offset;
  }
};

template <class MatrixType>
struct BsrRowViewConst {
  //! The type of the values in the row.
  typedef const typename MatrixType::non_const_value_type value_type;
  //! The type of the column indices in the row.
  typedef const typename MatrixType::non_const_ordinal_type ordinal_type;
  //! The type for returned block of values.
  typedef Kokkos::View<value_type**, Kokkos::LayoutRight,
                       typename MatrixType::device_type,
                       Kokkos::MemoryUnmanaged>
      block_values_type;

 private:
  //! Array of values in the row.
  value_type* values_;
  //! Array of (local) column indices in the row.
  ordinal_type* colidx_;
  /// \brief Stride between successive rows in a block
  const ordinal_type blockDim_;

 public:
  /// \brief Constructor
  ///
  /// \param values [in] Array of the row's values.
  /// \param colidx [in] Array of the row's column indices.
  /// \param blockDim [in] (Constant) stride between block rows
  ///   within a block-row in the above arrays.
  /// \param count [in] Number of entries in the row.
  //
  // Assumes values and colidx already offset to the correct location
  KOKKOS_INLINE_FUNCTION
  BsrRowViewConst(value_type* const values, ordinal_type* const colidx,
                  const ordinal_type& blockDim, const ordinal_type& count)
      : values_(values), colidx_(colidx), blockDim_(blockDim), length(count) {}

  /// \brief Constructor with offset into \c colidx array
  ///
  /// \param values [in] Array of the row's values.
  /// \param colidx [in] Array of the row's column indices.
  /// \param count [in] Number of entries in the row.
  /// \param start [in] Offset into values and colidx of the desired block-row
  /// start.
  ///   Note: The offset into the values array for a block-row equals
  ///           num_blocks_prior_to_block-row*blockDim*blockDim
  ///
  /// \tparam OffsetType The type of \c start (see above).  Must be a
  ///   built-in integer type.  This may differ from ordinal_type.
  ///   For example, the matrix may have dimensions that fit in int,
  ///   but a number of entries that does not fit in int.
  template <class OffsetType>
  KOKKOS_INLINE_FUNCTION BsrRowViewConst(
      const typename MatrixType::values_type& values,
      const typename MatrixType::index_type& colidx,
      const ordinal_type& blockDim, const ordinal_type& count,
      const OffsetType& start,
      const typename std::enable_if<std::is_integral<OffsetType>::value,
                                    int>::type& = 0)
      : values_(&values(start * blockDim * blockDim)),
        colidx_(&colidx(start)),
        blockDim_(blockDim),
        length(count) {}

  /// \brief Number of entries (i.e. blocks) in the row.
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

  /// /brief Return a pointer offset to local row i of block K of values_ array;
  ///        user responsible for indexing into this pointer correctly
  /// \param K [in] must be the LOCAL block index within this block-row
  /// \param i [in] must be the LOCAL row index offset within this block-row
  ///
  /// Output: pointer to values_ array at start of local row within block K
  ///
  /// Pointer interfaces are NOT guaranteed for backward compatibility
  /// This interface is intended for performant kernels, not common usage
  KOKKOS_INLINE_FUNCTION
  value_type* local_row_in_block(const ordinal_type& K,
                                 const ordinal_type& i) const {
    return (values_ + (K * blockDim_ * blockDim_ + i * blockDim_));
  }

  /// \brief Return the value at a specified block K with local row and col ids
  /// (i,j) \param K [in] must be the LOCAL block index within this block-row
  /// \param i [in] must be the LOCAL row index offset within this block-row
  /// \param j [in] must be the LOCAL col index offset within this block-row
  ///
  /// Output: reference to value_type at the given (K, i, j) offset into values_
  KOKKOS_INLINE_FUNCTION
  value_type& local_block_value(const ordinal_type& K, const ordinal_type& i,
                                const ordinal_type& j) const {
    return values_[K * blockDim_ * blockDim_ + i * blockDim_ + j];
  }

  /// \brief Return the block column index for a specified block K
  ///
  /// \param K [in] must be the LOCAL block index within this block-row
  /// \return Block column index for "uncompressed" block row
  KOKKOS_INLINE_FUNCTION
  ordinal_type block_colidx(const ordinal_type K) const { return colidx_[K]; }

  /// \brief Return unmanaged 2D strided View wrapping local block K from this
  /// block-row \param K [in] must be the LOCAL block index within this
  /// block-row
  KOKKOS_INLINE_FUNCTION
  block_values_type block(const ordinal_type& K) const {
    return block_values_type(&(values_[K * blockDim_ * blockDim_]),
                             Kokkos::LayoutRight(blockDim_, blockDim_));
  }

  /// \brief Return offset into colidx_ for the requested block idx
  ///        If none found, return Kokkos::Details::ArithTraits::max
  /// \param idx_to_match [in] local block idx within block-row
  /// \param is_sorted [in] defaulted to false; no usage at this time
  KOKKOS_INLINE_FUNCTION
  ordinal_type findRelBlockOffset(const ordinal_type& idx_to_match,
                                  bool /*is_sorted*/ = false) const {
    typedef typename std::remove_cv<ordinal_type>::type non_const_ordinal_type;
    non_const_ordinal_type offset =
        Kokkos::Details::ArithTraits<non_const_ordinal_type>::max();
    for (non_const_ordinal_type blk_offset = 0; blk_offset < length;
         ++blk_offset) {
      ordinal_type idx = colidx_[blk_offset];
      if (idx == idx_to_match) {
        offset = blk_offset;
        break;
      }  // return relative offset
    }
    return offset;
  }
};

/// \class BsrMatrix
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
template <class ScalarType, class OrdinalType, class Device,
          class MemoryTraits = void,
          class SizeType     = typename Kokkos::ViewTraits<OrdinalType*, Device,
                                                       void, void>::size_type>
class BsrMatrix {
  static_assert(
      std::is_signed<OrdinalType>::value,
      "BsrMatrix requires that OrdinalType is a signed integer type.");

 private:
  typedef
      typename Kokkos::ViewTraits<ScalarType*, Device, void,
                                  void>::host_mirror_space host_mirror_space;

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
  typedef BsrMatrix<ScalarType, OrdinalType, host_mirror_space, MemoryTraits>
      HostMirror;
  //! Type of the graph structure of the sparse matrix.
  typedef Kokkos::StaticCrsGraph<ordinal_type, Kokkos::LayoutLeft,
                                 execution_space, memory_traits, size_type>
      StaticCrsGraphType;
  //! Type of the graph structure of the sparse matrix - consistent with Kokkos.
  typedef Kokkos::StaticCrsGraph<ordinal_type, Kokkos::LayoutLeft,
                                 execution_space, memory_traits, size_type>
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
  typedef Kokkos::View<value_type*, Kokkos::LayoutRight, device_type,
                       MemoryTraits>
      values_type;
  //! Const version of the type of the entries in the sparse matrix.
  typedef typename values_type::const_value_type const_value_type;
  //! Nonconst version of the type of the entries in the sparse matrix.
  typedef typename values_type::non_const_value_type non_const_value_type;

  /// \name Storage of the actual sparsity structure and values.
  ///
  /// BsrMatrix uses the compressed sparse row (CSR) storage format to
  /// store the sparse matrix.  CSR is also called "compressed row
  /// storage"; hence the name, which it inherits from Tpetra and from
  /// Epetra before it.
  //@{
  //! The graph (sparsity structure) of the sparse matrix.
  staticcrsgraph_type graph;
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
  /// mfh: numCols and nnz should be properties of the graph, not the matrix.
  /// Then BsrMatrix needs methods to get these from the graph.
  BsrMatrix() : graph(), values(), dev_config(), numCols_(0), blockDim_(1) {}

  //! Copy constructor (shallow copy).
  template <typename SType, typename OType, class DType, class MTType,
            typename IType>
  explicit BsrMatrix(const BsrMatrix<SType, OType, DType, MTType, IType>& B)
      : graph(B.graph.entries, B.graph.row_map),
        values(B.values),
        dev_config(B.dev_config),
        numCols_(B.numCols()),
        blockDim_(B.blockDim()) {
    graph.row_block_offsets = B.graph.row_block_offsets;
    // MD: Changed the copy constructor of graph
    // as the constructor of StaticCrsGraph does not allow copy from non const
    // version.
  }

  /// \brief Construct with a graph that will be shared.
  ///
  /// \param[in] arg_label   The sparse matrix's label.
  /// \param[in] arg_graph   The graph between the blocks.
  /// \param[in] blockDimIn  The block size.
  ///
  /// Allocate the values array for subsequent fill.
  BsrMatrix(const std::string& arg_label, const staticcrsgraph_type& arg_graph,
            const OrdinalType& blockDimIn)
      : graph(arg_graph),
        values(arg_label,
               arg_graph.entries.extent(0) * blockDimIn * blockDimIn),
        numCols_(maximum_entry(arg_graph) + 1),
        blockDim_(blockDimIn) {
    if (blockDim_ < 1) {
      std::ostringstream os;
      os << "KokkosSparse::BsrMatrix: Inappropriate block size: " << blockDim_;
      KokkosKernels::Impl::throw_runtime_exception(os.str());
    }
  }

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
  /// The \c pad argument is currently not used.
  BsrMatrix(const std::string& label, OrdinalType nrows, OrdinalType ncols,
            size_type annz, ScalarType* val, OrdinalType* rows,
            OrdinalType* cols, OrdinalType blockdim, bool pad = false) {
    (void)label;
    (void)pad;
    blockDim_ = blockdim;

    if (blockDim_ < 1) {
      std::ostringstream os;
      os << "KokkosSparse::BsrMatrix: Inappropriate block size: " << blockDim_;
      KokkosKernels::Impl::throw_runtime_exception(os.str());
    }

    if ((ncols % blockDim_ != 0) || (nrows % blockDim_ != 0)) {
      assert(
          (ncols % blockDim_ == 0) &&
          "BsrMatrix: input CrsMatrix columns is not a multiple of block size");
      assert((nrows % blockDim_ == 0) &&
             "BsrMatrix: input CrsMatrix rows is not a multiple of block size");
    }

    numCols_                  = ncols / blockDim_;
    ordinal_type tmp_num_rows = nrows / blockDim_;

    //
    // Wrap the raw pointers in unmanaged host Views
    // Note that the inputs are in coordinate format.
    // So unman_rows and unman_cols have the same type.
    //
    typename values_type::HostMirror unman_val(val, annz);
    typename index_type::HostMirror unman_rows(rows, annz);
    typename index_type::HostMirror unman_cols(cols, annz);

    typename row_map_type::non_const_type tmp_row_map(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "rowmap"),
        tmp_num_rows + 1);
    auto row_map_host = Kokkos::create_mirror_view(tmp_row_map);
    Kokkos::deep_copy(row_map_host, 0);

    if (annz > 0) {
      ordinal_type iblock = 0;
      std::set<ordinal_type> set_blocks;
      for (size_type ii = 0; ii <= annz; ++ii) {
        if ((ii == annz) || ((unman_rows(ii) / blockDim_) > iblock)) {
          // Flush the stored entries
          row_map_host(iblock + 1) = set_blocks.size();
          if (ii == annz) break;
          set_blocks.clear();
          iblock = unman_rows(ii) / blockDim_;
        }
        ordinal_type tmp_jblock = unman_cols(ii) / blockDim_;
        set_blocks.insert(tmp_jblock);
      }
    }

    for (size_type ii = 0; ii < annz; ++ii)
      row_map_host(ii + 1) += row_map_host(ii);

    Kokkos::deep_copy(tmp_row_map, row_map_host);

    // Create temporary Views for row_map and entries
    // because the StaticCrsGraph ctor requires View inputs
    index_type tmp_entries("tmp_entries", row_map_host(tmp_num_rows));
    auto tmp_entries_host = Kokkos::create_mirror_view(tmp_entries);

    Kokkos::resize(values, row_map_host(tmp_num_rows) * blockDim_ * blockDim_);
    auto values_host = Kokkos::create_mirror_view(values);
    Kokkos::deep_copy(values_host, 0);

    if (annz > 0) {
      //--- Fill tmp_entries
      ordinal_type cur_block = 0;
      std::set<ordinal_type> set_blocks;
      for (size_type ii = 0; ii <= annz; ++ii) {
        if ((ii == annz) || ((unman_rows(ii) / blockDim_) > cur_block)) {
          // Flush the stored entries
          ordinal_type ipos = row_map_host(cur_block);
          for (auto jblock : set_blocks) tmp_entries_host(ipos++) = jblock;
          if (ii == annz) break;
          set_blocks.clear();
          cur_block = unman_rows(ii) / blockDim_;
        }
        ordinal_type tmp_jblock = unman_cols(ii) / blockDim_;
        set_blocks.insert(tmp_jblock);
      }
      //--- Fill numerical values
      for (size_type ii = 0; ii < annz; ++ii) {
        const auto ilocal = unman_rows(ii) % blockDim_;
        const auto jblock = unman_cols(ii) / blockDim_;
        const auto jlocal = unman_cols(ii) % blockDim_;
        for (auto jj = row_map_host(jblock); jj < row_map_host(jblock + 1);
             ++jj) {
          if (tmp_entries_host(jj) == jblock) {
            const auto shift =
                jj * blockDim_ * blockDim_ + ilocal * blockDim_ + jlocal;
            values_host(shift) = unman_val(ii);
            break;
          }
        }
      }
    }

    Kokkos::deep_copy(tmp_entries, tmp_entries_host);
    Kokkos::deep_copy(values, values_host);

    // Initialize graph using the temp entries and row_map Views
    graph = staticcrsgraph_type(tmp_entries, tmp_row_map);
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
  BsrMatrix(const std::string& /*label*/, const OrdinalType nrows,
            const OrdinalType ncols, const size_type /*annz*/,
            const values_type& vals, const row_map_type& rows,
            const index_type& cols, const OrdinalType blockDimIn)
      : graph(cols, rows),
        values(vals),
        numCols_(ncols),
        blockDim_(blockDimIn) {
    if (blockDim_ < 1) {
      std::ostringstream os;
      os << "KokkosSparse::BsrMatrix: Inappropriate block size: " << blockDim_;
      KokkosKernels::Impl::throw_runtime_exception(os.str());
    }

    const ordinal_type actualNumRows =
        (rows.extent(0) != 0) ? static_cast<ordinal_type>(
                                    rows.extent(0) - static_cast<size_type>(1))
                              : static_cast<ordinal_type>(0);
    if (nrows != actualNumRows) {
      std::ostringstream os;
      os << "Input argument nrows = " << nrows
         << " != the actual number of "
            "rows "
         << actualNumRows << " according to the 'rows' input argument.";
      throw std::invalid_argument(os.str());
    }
    // nnz returns graph.entries.extent(0) i.e. ptr[ nrows + 1 ] nnz entry
    // input annz is nnz of values, not comparable with block ptr 'nnz' i.e.
    // numBlocks
    if (blockDim_ <= 0) {
      std::ostringstream os;
      os << "Input argument blockDim = " << blockDim_
         << " is not larger than 0.";
      throw std::invalid_argument(os.str());
    }
  }

  /// \brief Constructor that accepts a a static graph, and values.
  ///
  /// The matrix will store and use the row map, indices, and values
  /// directly (by view, not by deep copy).
  ///
  /// \param[in] label  The sparse matrix's label.
  /// \param[in] ncols  The number of columns.
  /// \param[in] vals   The entries.
  /// \param[in] graph_ The graph between the blocks.
  /// \param[in] blockDimIn  The block size.
  BsrMatrix(const std::string& /*label*/, const OrdinalType& ncols,
            const values_type& vals, const staticcrsgraph_type& graph_,
            const OrdinalType& blockDimIn)
      : graph(graph_), values(vals), numCols_(ncols), blockDim_(blockDimIn) {
    if (blockDim_ < 1) {
      std::ostringstream os;
      os << "KokkosSparse::BsrMatrix: Inappropriate block size: " << blockDim_;
      KokkosKernels::Impl::throw_runtime_exception(os.str());
    }
  }

  /// \brief Constructor that accepts a CrsMatrix and block dimension,
  ///        assuming the provided CrsMatrix has appropriate block structure.
  template <typename SType, typename OType, class DType, class MTType,
            typename IType>
  BsrMatrix(const KokkosSparse::CrsMatrix<SType, OType, DType, MTType, IType>&
                crs_mtx,
            const OrdinalType blockDimIn) {
    typedef typename KokkosSparse::CrsMatrix<SType, OType, DType, MTType, IType>
        crs_matrix_type;
    typedef typename crs_matrix_type::staticcrsgraph_type crs_graph_type;
    typedef typename crs_graph_type::entries_type crs_graph_entries_type;
    typedef typename crs_graph_type::row_map_type crs_graph_row_map_type;

    blockDim_ = blockDimIn;
    if (blockDim_ < 1) {
      std::ostringstream os;
      os << "KokkosSparse::BsrMatrix: Inappropriate block size: " << blockDim_;
      KokkosKernels::Impl::throw_runtime_exception(os.str());
    }

    assert(
        (crs_mtx.numCols() % blockDim_ == 0) &&
        "BsrMatrix: input CrsMatrix columns is not a multiple of block size");
    assert((crs_mtx.numRows() % blockDim_ == 0) &&
           "BsrMatrix: input CrsMatrix rows is not a multiple of block size");

    numCols_ = crs_mtx.numCols() / blockDim_;

    OrdinalType nbrows =
        crs_mtx.numRows() /
        blockDim_;  // actual number of block rows; add 1 for ptr length

    // block_rows will accumulate the number of blocks per row - this is NOT the
    // row_map with cum sum!!
    std::vector<OrdinalType> block_rows(nbrows, 0);

    typename crs_graph_row_map_type::HostMirror h_crs_row_map =
        Kokkos::create_mirror_view(crs_mtx.graph.row_map);
    Kokkos::deep_copy(h_crs_row_map, crs_mtx.graph.row_map);
    typename crs_graph_entries_type::HostMirror h_crs_entries =
        Kokkos::create_mirror_view(crs_mtx.graph.entries);
    Kokkos::deep_copy(h_crs_entries, crs_mtx.graph.entries);

    // determine size of block cols indices == number of blocks,
    // i.e. nnz for the block CRS graph
    OrdinalType numBlocks = 0;
    for (OrdinalType i = 0; i < crs_mtx.numRows(); i += blockDim_) {
      std::set<OrdinalType> col_set;
      for (auto ie = h_crs_row_map(i); ie < h_crs_row_map(i + blockDim_);
           ++ie) {
        col_set.insert(h_crs_entries(ie) / blockDim_);
      }
      numBlocks += col_set.size();                 // cum sum
      block_rows[i / blockDim_] = col_set.size();  // frequency counts
    }

    // create_staticcrsgraph takes the frequency of blocks per row
    // and returns the cum sum pointer row_map with nbrows+1 size, and total
    // numBlocks in the final entry
    graph = Kokkos::create_staticcrsgraph<staticcrsgraph_type>("blockgraph",
                                                               block_rows);
    typename row_map_type::HostMirror h_row_map =
        Kokkos::create_mirror_view(graph.row_map);
    Kokkos::deep_copy(h_row_map, graph.row_map);

    typename index_type::HostMirror h_entries =
        Kokkos::create_mirror_view(graph.entries);

    OrdinalType ientry = 0;
    for (OrdinalType ib = 0; ib < nbrows; ++ib) {
      auto ir_start = ib * blockDim_;
      auto ir_stop  = (ib + 1) * blockDim_;
      std::set<OrdinalType> col_set;
      for (auto jk = h_crs_row_map(ir_start); jk < h_crs_row_map(ir_stop);
           ++jk) {
        col_set.insert(h_crs_entries(jk) / blockDim_);
      }
      for (auto col_block : col_set) {
        h_entries(ientry++) = col_block;
      }
    }
    Kokkos::deep_copy(graph.entries, h_entries);

    // Copy the numerical values

    typename values_type::HostMirror h_crs_values =
        Kokkos::create_mirror_view(crs_mtx.values);
    Kokkos::deep_copy(h_crs_values, crs_mtx.values);

    typename values_type::HostMirror h_values =
        Kokkos::create_mirror_view(values);
    if (h_values.extent(0) < size_t(numBlocks * blockDim_ * blockDim_)) {
      Kokkos::resize(h_values, numBlocks * blockDim_ * blockDim_);
      Kokkos::resize(values, numBlocks * blockDim_ * blockDim_);
    }
    Kokkos::deep_copy(h_values, 0);

    for (OrdinalType ir = 0; ir < crs_mtx.numRows(); ++ir) {
      const auto iblock = ir / blockDim_;
      const auto ilocal = ir % blockDim_;
      for (auto jk = h_crs_row_map(ir); jk < h_crs_row_map(ir + 1); ++jk) {
        const auto jc     = h_crs_entries(jk);
        const auto jblock = jc / blockDim_;
        const auto jlocal = jc % blockDim_;
        for (auto jkb = h_row_map(iblock); jkb < h_row_map(iblock + 1); ++jkb) {
          if (h_entries(jkb) == jblock) {
            OrdinalType shift = jkb * blockDim_ * blockDim_;
            h_values(shift + ilocal * blockDim_ + jlocal) = h_crs_values(jk);
            break;
          }
        }
      }
    }
    Kokkos::deep_copy(values, h_values);
  }

  /// \brief Given an array of blocks, sum the values into corresponding
  ///        block in BsrMatrix
  /// \param[in] rowi    is a block-row index
  /// \param[in] ncol  is number of blocks referenced in cols[] array
  /// \param[in] cols[] are block colidxs within the block-row to be summed
  /// into ncol entries
  /// \param[in] vals[] array containing 'block' of values
  ///        ncol*block_size*block_size entries
  ///        assume vals block is provided in 'LayoutRight' or 'Row Major'
  ///        format, that is e.g. 2x2 block [ a b ; c d ] provided as flattened
  ///        1d array as [a b c d] Assume that each block is stored contiguously
  ///        in vals: [a b; c d] [e f; g h] -> [a b c d e f g h] If so, then i
  ///        in [0, ncols) for cols[] maps to i*block_size*block_size in vals[]
  KOKKOS_INLINE_FUNCTION
  OrdinalType sumIntoValues(const OrdinalType rowi, const OrdinalType cols[],
                            const OrdinalType ncol, const ScalarType vals[],
                            const bool is_sorted    = false,
                            const bool force_atomic = false) const {
    return operateValues(BsrMatrix::valueOperation::ADD, rowi, cols, ncol, vals,
                         is_sorted, force_atomic);
  }

  /// \brief Given an array of blocks, replace the values of corresponding
  ///        blocks in BsrMatrix
  /// \param[in] rowi    is a block-row index
  /// \param[in] ncol is number of blocks referenced in cols[] array
  /// \param[in] cols[] are block colidxs within the block-row to be summed
  /// into ncol entries
  /// \param vals[] [in] array containing 'block' of values
  //        ncol*block_size*block_size entries
  //        assume vals block is provided in 'LayoutRight' or 'Row Major'
  //        format, that is e.g. 2x2 block [ a b ; c d ] provided as flattened
  //        1d array as [a b c d] Assume that each block is stored contiguously
  //        in vals: [a b; c d] [e f; g h] -> [a b c d e f g h] If so, then i in
  //        [0, ncols) for cols[] maps to i*block_size*block_size in vals[]
  KOKKOS_INLINE_FUNCTION
  OrdinalType replaceValues(const OrdinalType rowi, const OrdinalType cols[],
                            const OrdinalType ncol, const ScalarType vals[],
                            const bool is_sorted    = false,
                            const bool force_atomic = false) const {
    return operateValues(BsrMatrix::valueOperation::ASSIGN, rowi, cols, ncol,
                         vals, is_sorted, force_atomic);
  }

  //! Attempt to assign the input matrix to \c *this.
  // Are the CUDA sparse handles needed to be copied here??
  template <typename aScalarType, typename aOrdinalType, class aDevice,
            class aMemoryTraits, typename aSizeType>
  BsrMatrix& operator=(const BsrMatrix<aScalarType, aOrdinalType, aDevice,
                                       aMemoryTraits, aSizeType>& mtx) {
    numCols_   = mtx.numCols();
    blockDim_  = mtx.blockDim();
    graph      = mtx.graph;
    values     = mtx.values;
    dev_config = mtx.dev_config;
    return *this;
  }

  //! The number of rows in the sparse matrix.
  KOKKOS_INLINE_FUNCTION ordinal_type numRows() const {
    return graph.numRows();
  }

  //! The number of columns in the sparse matrix.
  KOKKOS_INLINE_FUNCTION ordinal_type numCols() const { return numCols_; }

  //! The block dimension in the sparse block matrix.
  KOKKOS_INLINE_FUNCTION ordinal_type blockDim() const { return blockDim_; }

  //! The number of "point" (non-block) rows in the matrix.
  //  This is the dimension of the range of this matrix as a linear operator.
  KOKKOS_INLINE_FUNCTION ordinal_type numPointRows() const {
    return numRows() * blockDim();
  }

  //! The number of "point" (non-block) columns in the matrix.
  //  This is the dimension of the domain of this matrix as a linear operator.
  KOKKOS_INLINE_FUNCTION ordinal_type numPointCols() const {
    return numCols() * blockDim();
  }

  //! The number of stored entries in the sparse matrix.
  KOKKOS_INLINE_FUNCTION size_type nnz() const {
    return graph.entries.extent(0);
  }

  friend struct BsrRowView<BsrMatrix>;

  /// \brief Return a BsrRowView of block-row i of the matrix.
  ///
  /// If row i does not belong to the matrix, return an empty view.
  ///
  /// The returned object \c view implements the following interface:
  /// <ul>
  /// <li> \c view.length is the number of entries (i.e. blocks)
  ///      in the block row </li>
  /// <li> \c view.local_row_in_block_row(K, i) returns a nonconst pointer
  ///      to the values of the ith local row in the k-th block of the block-row
  ///      </li>
  /// <li> \c view.full_row_in_block_row(i) returns a nonconst pointer
  ///      to the values of the ith local row of the block-row </li>
  /// <li> \c view.local_block_value(K, i, j) returns a nonconst reference
  ///      to the value in the ith local row and jth local col
  ///      of the k-th block of the block-row </li>
  /// <li> \c view.block(K) returns an unmanaged 2D strided Kokkos::View
  ///      of the values of the k-th block of the block-row </li>
  /// </ul>
  ///
  /// Users should not rely on the return type of this method.  They
  /// should instead assign to 'auto'.
  ///
  KOKKOS_INLINE_FUNCTION
  BsrRowView<BsrMatrix> block_row(const ordinal_type i) const {
    const size_type start =
        graph.row_map(i);  // total num blocks prior to this block-row
    const auto count = static_cast<ordinal_type>(
        graph.row_map(i + 1) - start);  // num blocks in this row

    if (count == 0) {
      return BsrRowView<BsrMatrix>(nullptr, nullptr, 1, 0);
    } else {
      return BsrRowView<BsrMatrix>(values, graph.entries, blockDim(), count,
                                   start);
    }
  }

  /// \brief Return a BsrRowViewConst of block-row i of the matrix.
  ///
  /// If row i does not belong to the matrix, return an empty view.
  ///
  /// The returned object \c view implements the following interface:
  /// <ul>
  /// <li> \c view.length is the number of entries (i.e. blocks)
  ///      in the block row </li>
  /// <li> \c view.local_row_in_block_row(K, i) returns a nonconst pointer
  ///      to the values of the ith local row in the k-th block of the block-row
  ///      </li>
  /// <li> \c view.full_row_in_block_row(i) returns a nonconst pointer
  ///      to the values of the ith local row of the block-row </li>
  /// <li> \c view.local_block_value(K, i, j) returns a nonconst reference
  ///      to the value in the ith local row and jth local col
  ///      of the k-th block of the block-row </li>
  /// <li> \c view.block(K) returns an unmanaged 2D strided Kokkos::View
  ///      of the values of the k-th block of the block-row </li>
  /// </ul>
  ///
  /// Users should not rely on the return type of this method.  They
  /// should instead assign to 'auto'.
  ///
  KOKKOS_INLINE_FUNCTION
  BsrRowViewConst<BsrMatrix> block_row_Const(const ordinal_type i) const {
    const size_type start =
        graph.row_map(i);  // total num blocks prior to this block-row
    const auto count = static_cast<ordinal_type>(
        graph.row_map(i + 1) - start);  // num blocks in this row

    if (count == 0) {
      return BsrRowViewConst<BsrMatrix>(nullptr, nullptr, 1, 0);
    } else {
      return BsrRowViewConst<BsrMatrix>(values, graph.entries, blockDim(),
                                        count, start);
    }
  }

 protected:
  enum class valueOperation { ADD, ASSIGN };

  /// \brief Given an array of blocks, operate on the values of corresponding
  ///        blocks in BsrMatrix
  /// \param[in] rowi    is a block-row index
  /// \param[in] ncol is number of blocks referenced in cols[] array
  /// \param[in] cols[] are block colidxs within the block-row to be op-ed
  /// into ncol entries
  /// \param vals[] [in] array containing 'block' of values
  //        ncol*block_size*block_size entries
  //        assume vals block is provided in 'LayoutRight' or 'Row Major'
  //        format, that is e.g. 2x2 block [ a b ; c d ] provided as flattened
  //        1d array as [a b c d] Assume that each block is stored contiguously
  //        in vals: [a b; c d] [e f; g h] -> [a b c d e f g h] If so, then i in
  //        [0, ncols) for cols[] maps to i*block_size*block_size in vals[]
  KOKKOS_INLINE_FUNCTION
  OrdinalType operateValues(const BsrMatrix::valueOperation op,
                            const OrdinalType rowi, const OrdinalType cols[],
                            const OrdinalType ncol, const ScalarType vals[],
                            const bool is_sorted    = false,
                            const bool force_atomic = false) const {
    BsrRowView<BsrMatrix> row_view = this->block_row(rowi);
    const ordinal_type block_size  = this->blockDim();

    ordinal_type numValid = 0;  // number of valid local column indices

    for (ordinal_type i = 0; i < ncol; ++i) {
      // Find offset into values for block-row rowi and colidx cols[i]
      // cols[i] is the index to match
      // blk_offset is the offset for block colidx from bptr[rowi] to bptr[rowi
      // + 1] (not global offset) colidx_ and values_ are already offset to the
      // beginning of blockrow rowi
      auto blk_offset = row_view.findRelBlockOffset(cols[i], is_sorted);
      if (blk_offset != Kokkos::Details::ArithTraits<ordinal_type>::max()) {
        ordinal_type offset_into_vals =
            i * block_size *
            block_size;  // stride == 1 assumed between elements
        for (ordinal_type lrow = 0; lrow < block_size; ++lrow) {
          auto local_row_values = row_view.local_row_in_block(
              blk_offset, lrow);  // pointer to start of specified local row
          // within this block
          switch (op) {
            case BsrMatrix::valueOperation::ADD: {
              for (ordinal_type lcol = 0; lcol < block_size; ++lcol) {
                if (force_atomic) {
                  Kokkos::atomic_add(
                      &(local_row_values[lcol]),
                      vals[offset_into_vals + lrow * block_size + lcol]);
                } else {
                  local_row_values[lcol] +=
                      vals[offset_into_vals + lrow * block_size + lcol];
                }
              }
              break;
            }
            case BsrMatrix::valueOperation::ASSIGN: {
              for (ordinal_type lcol = 0; lcol < block_size; ++lcol) {
                if (force_atomic) {
                  Kokkos::atomic_assign(
                      &(local_row_values[lcol]),
                      vals[offset_into_vals + lrow * block_size + lcol]);
                } else {
                  local_row_values[lcol] =
                      vals[offset_into_vals + lrow * block_size + lcol];
                }
              }
              break;
            }
          }
        }
        ++numValid;
      }
    }  // end for ncol
    return numValid;
  }

 private:
  ordinal_type numCols_  = 0;
  ordinal_type blockDim_ = 1;  // TODO Assuming square blocks for now
};

//----------------------------------------------------------------------------
/// \class is_bsr_matrix
/// \brief is_bsr_matrix<T>::value is true if T is a BsrMatrix<...>, false
/// otherwise
template <typename>
struct is_bsr_matrix : public std::false_type {};
template <typename... P>
struct is_bsr_matrix<BsrMatrix<P...>> : public std::true_type {};
template <typename... P>
struct is_bsr_matrix<const BsrMatrix<P...>> : public std::true_type {};
//----------------------------------------------------------------------------

}  // namespace Experimental
}  // namespace KokkosSparse
#endif
