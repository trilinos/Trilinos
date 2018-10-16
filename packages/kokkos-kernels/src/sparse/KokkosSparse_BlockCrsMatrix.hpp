/*
//@HEADER
// ************************************************************************
//
//               KokkosKernels 0.9: Linear Algebra and Graph Kernels
//                 Copyright 2017 Sandia Corporation
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
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

/// \file Kokkos_Sparse_BlockCrsMatrix.hpp
/// \brief Local sparse matrix interface
///
/// This file provides KokkosSparse::BlockCrsMatrix.  This implements a
/// local (no MPI) sparse matrix stored in block compressed row sparse
/// ("BlockCrs") format.

#ifndef KOKKOS_SPARSE_BLOCKCRSMATRIX_HPP_
#define KOKKOS_SPARSE_BLOCKCRSMATRIX_HPP_

#include "Kokkos_Core.hpp"
#include "Kokkos_StaticCrsGraph.hpp"
#include "Kokkos_ArithTraits.hpp"
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include "KokkosSparse_CrsMatrix.hpp"

namespace KokkosSparse {

namespace Experimental {

/// \class SparseBlockRowView
/// \brief View of a block-row of a sparse matrix.
/// \tparam MatrixType BlockCrsMatrix Sparse matrix type
///
/// This class provides a generic view of a block-row of a sparse matrix.
///
/// Whether the view is const or not, depends on whether
/// MatrixType is a const or nonconst view of the matrix.  If
/// you always want a const view, use SparseBlockRowViewConst (see below).
///
/// MatrixType must provide the \c value_type and \c ordinal_type
/// typedefs.  In addition, it must make sense to use SparseBlockRowView to
/// view a block-row of MatrixType.  
template<class MatrixType>
struct SparseBlockRowView {
  //! The type of the values in the row.
  typedef typename MatrixType::value_type value_type;
  //! The type of the column indices in the row.
  typedef typename MatrixType::ordinal_type ordinal_type;
  //! The type for returned block of values. 
  typedef Kokkos::View< value_type**, Kokkos::LayoutStride, typename MatrixType::device_type, Kokkos::MemoryUnmanaged > block_values_type;

private:
  //! Array of values in the row.
  value_type* values_;
  //! Array of (local) column indices in the row.
  ordinal_type* colidx_;
  /// \brief Stride between successive rows in a block.
  ///
  /// For block compressed sparse row (BlockCSR) storage with row-major layout by full row,
  /// (i.e. consecutive rows within a block are NOT contiguous), this will be the stride 
  /// between rows within a block-row
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
  // Assumes values and colidx__ already offset to the correct location
  KOKKOS_INLINE_FUNCTION
  SparseBlockRowView (value_type* const values,
                      ordinal_type* const colidx__,
                      const ordinal_type& blockDim,
                      const ordinal_type& count) :
    values_ (values), colidx_ (colidx__), blockDim_(blockDim), length (count)
  {}

  /// \brief Constructor with offset into \c colidx array
  ///
  /// \param values [in] Array of the row's values.
  /// \param colidx [in] Array of the row's column indices.
  /// \param blockDim [in] (Constant) stride between rows in
  ///   within a block in the above arrays.
  /// \param count [in] Number of blocks in the desired block-row
  /// \param start [in] Offset into values and colidx of the desired block-row start.
  ///   Note: The offset into the values array for a block-row equals
  ///           num_blocks_prior_to_block-row*blockDim*blockDim
  ///
  /// \tparam OffsetType The type of \c start (see above).  Must be a
  ///   built-in integer type.  This may differ from ordinal_type.
  ///   For example, the matrix may have dimensions that fit in int,
  ///   but a number of entries that does not fit in int.
  template<class OffsetType>
  KOKKOS_INLINE_FUNCTION
  SparseBlockRowView (const typename MatrixType::values_type& values,
                      const typename MatrixType::index_type& colidx__,
                      const ordinal_type& blockDim,
                      const ordinal_type& count,
                      const OffsetType& start,
                      const typename std::enable_if<std::is_integral<OffsetType>::value, int>::type& = 0) :
    values_ (&values(start*blockDim*blockDim)), colidx_ (&colidx__(start)), blockDim_(blockDim), length (count)
  {}

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


  /// \brief Return a pointer offset to full-row i of values_ array; 
  ///        user responsible for indexing into this pointer correctly
  /// \param i [in] must be the LOCAL row index offset within this block-row
  ///
  /// Output: pointer to values_ array at start of full row with local index i
  ///
  /// Pointer interfaces are NOT guaranteed for backward compatibility
  /// This interface is intended for performant kernels, not common usage
  KOKKOS_INLINE_FUNCTION
  value_type* full_row_in_block_row (const ordinal_type& i) const {
    return values_+(i*length*blockDim_) ;
  }

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
  value_type* local_row_in_block (const ordinal_type& K, const ordinal_type& i) const {
    return (values_+(K*blockDim_ + i*length*blockDim_)) ;
  }

  /// \brief Return the value at a specified block K of block-row 
  ///        with local row and col offset (i,j)
  /// \param K [in] must be the LOCAL block index within this block-row
  /// \param i [in] must be the LOCAL row index offset within this block-row
  /// \param j [in] must be the LOCAL col index offset within this block-row
  ///
  /// Output: reference to value_type at the given (K, i, j) offset into values_
  KOKKOS_INLINE_FUNCTION
  value_type& local_block_value (const ordinal_type& K, const ordinal_type& i, const ordinal_type& j) const {
    return values_[K*blockDim_ + i*length*blockDim_ + j];
  }

  /// \brief Return unmanaged 2D strided View wrapping local block K from this block-row
  /// \param K [in] must be the LOCAL block index within this block-row
  KOKKOS_INLINE_FUNCTION
  block_values_type block(const ordinal_type& K) const {
    return block_values_type( &(values_[K*blockDim_]), Kokkos::LayoutStride(blockDim_,length*blockDim_,blockDim_,1) );
  }


  /// \brief Return offset into colidx_ for the requested block idx
  ///        If none found, return Kokkos::Details::ArithTraits::max
  /// \param idx_to_match [in] local block idx within block-row
  /// \param is_sorted [in] defaulted to false; no usage at this time
  KOKKOS_INLINE_FUNCTION
  ordinal_type findRelBlockOffset ( const ordinal_type idx_to_match, bool is_sorted = false ) const {
    ordinal_type offset = Kokkos::Details::ArithTraits< ordinal_type >::max();
    for ( ordinal_type blk_offset = 0; blk_offset < length; ++blk_offset ) {
      ordinal_type idx = colidx_[blk_offset];
      if ( idx == idx_to_match ) 
      { 
        offset = blk_offset; 
        break;
      } // return relative offset
    }
    return offset;
  }

};


/// \class SparseBlockRowViewConst
/// \brief Const view of a row of a sparse matrix.
/// \tparam MatrixType Sparse matrix type, such as BlockCrsMatrix.
///
/// This class is like SparseBlockRowView, except that it provides a const
/// view.  This class exists in order to let users get a const view of
/// a row of a nonconst matrix.
template<class MatrixType>
struct SparseBlockRowViewConst {
  //! The type of the values in the row.
  typedef const typename MatrixType::non_const_value_type value_type;
  //! The type of the column indices in the row.
  typedef const typename MatrixType::non_const_ordinal_type ordinal_type;
  //! The type for returned block of values. 
  typedef Kokkos::View< value_type**, Kokkos::LayoutStride, typename MatrixType::device_type, Kokkos::MemoryUnmanaged > block_values_type;

private:
  //! Array of values in the row.
  value_type* values_;
  //! Array of (local) column indices in the row.
  ordinal_type* colidx_;
  /// \brief Stride between successive rows in a block-row
  ///
  /// For block compressed sparse row (BlockCSR) storage with row-major layout,
  /// (i.e. consecutive rows within a block are NOT contiguous), this will be the stride 
  /// between rows within a block-row
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
  // Assumes values and colidx__ already offset to the correct location
  KOKKOS_INLINE_FUNCTION
  SparseBlockRowViewConst (value_type* const values,
                           ordinal_type* const colidx__,
                           const ordinal_type& blockDim,
                           const ordinal_type& count) :
    values_ (values), colidx_ (colidx__), blockDim_(blockDim), length (count)
  {}

  /// \brief Constructor with offset into \c colidx array
  ///
  /// \param values [in] Array of the row's values.
  /// \param colidx [in] Array of the row's column indices.
  /// \param count [in] Number of entries in the row.
  /// \param start [in] Offset into values and colidx of the desired block-row start.
  ///   Note: The offset into the values array for a block-row equals
  ///           num_blocks_prior_to_block-row*blockDim*blockDim
  ///
  /// \tparam OffsetType The type of \c start (see above).  Must be a
  ///   built-in integer type.  This may differ from ordinal_type.
  ///   For example, the matrix may have dimensions that fit in int,
  ///   but a number of entries that does not fit in int.
  template<class OffsetType>
  KOKKOS_INLINE_FUNCTION
  SparseBlockRowViewConst (const typename MatrixType::values_type& values,
                           const typename MatrixType::index_type& colidx__,
                           const ordinal_type& blockDim,
                           const ordinal_type& count,
                           const OffsetType& start,
                           const typename std::enable_if<std::is_integral<OffsetType>::value, int>::type& = 0) :
    values_ (&values(start*blockDim*blockDim)), colidx_ (&colidx__(start)), blockDim_(blockDim), length (count)
  {}

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


  /// \brief Return a pointer offset to full-row i of values_ array; 
  ///        user responsible for indexing into this pointer correctly
  /// \param i [in] must be the LOCAL row index offset within this block-row
  ///
  /// Output: pointer to values_ array at start of full row with local index i
  ///
  /// Pointer interfaces are NOT guaranteed for backward compatibility
  /// This interface is intended for performant kernels, not common usage
  KOKKOS_INLINE_FUNCTION
  value_type* full_row_in_block_row (const ordinal_type& i) const {
    return values_+(i*length*blockDim_) ;
  }

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
  value_type* local_row_in_block (const ordinal_type& K, const ordinal_type& i) const {
    return (values_+(K*blockDim_ + i*length*blockDim_)) ;
  }

  /// \brief Return the value at a specified block K with local row and col ids (i,j)
  /// \param K [in] must be the LOCAL block index within this block-row
  /// \param i [in] must be the LOCAL row index offset within this block-row
  /// \param j [in] must be the LOCAL col index offset within this block-row
  ///
  /// Output: reference to value_type at the given (K, i, j) offset into values_
  KOKKOS_INLINE_FUNCTION
  value_type& local_block_value (const ordinal_type& K, const ordinal_type& i, const ordinal_type& j) const {
    return values_[K*blockDim_ + i*length*blockDim_ + j];
  }

  /// \brief Return unmanaged 2D strided View wrapping local block K from this block-row
  /// \param K [in] must be the LOCAL block index within this block-row
  KOKKOS_INLINE_FUNCTION
  block_values_type block(const ordinal_type& K) const {
    return block_values_type( &(values_[K*blockDim_]), Kokkos::LayoutStride(blockDim_,length*blockDim_,blockDim_,1) );
  }


  /// \brief Return offset into colidx_ for the requested block idx
  ///        If none found, return Kokkos::Details::ArithTraits::max
  /// \param idx_to_match [in] local block idx within block-row
  /// \param is_sorted [in] defaulted to false; no usage at this time
  KOKKOS_INLINE_FUNCTION
  ordinal_type findRelBlockOffset ( const ordinal_type &idx_to_match, bool is_sorted = false ) const {
    typedef typename std::remove_cv<ordinal_type>::type non_const_ordinal_type;
    non_const_ordinal_type offset = Kokkos::Details::ArithTraits< non_const_ordinal_type >::max();
    for ( non_const_ordinal_type blk_offset = 0; blk_offset < length; ++blk_offset ) {
      ordinal_type idx = colidx_[blk_offset];
      if ( idx == idx_to_match ) 
      { 
        offset = blk_offset; 
        break;
      } // return relative offset
    }
    return offset;
  }

};

/// \class BlockCrsMatrix
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
class BlockCrsMatrix {
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
  typedef BlockCrsMatrix<ScalarType, OrdinalType, host_mirror_space, MemoryTraits> HostMirror;
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
  //! Type of the graph structure of the sparse matrix.
  typedef Kokkos::StaticCrsGraph<ordinal_type, Kokkos::LayoutLeft, execution_space, size_type, memory_traits> StaticCrsGraphType;
  //! Type of the graph structure of the sparse matrix - consistent with Kokkos.
  typedef Kokkos::StaticCrsGraph<ordinal_type, Kokkos::LayoutLeft, execution_space, size_type, memory_traits> staticcrsgraph_type;
#else
  //! Type of the graph structure of the sparse matrix.
  typedef Kokkos::StaticCrsGraph<ordinal_type, Kokkos::LayoutLeft, execution_space, memory_traits, size_type> StaticCrsGraphType;
  //! Type of the graph structure of the sparse matrix - consistent with Kokkos.
  typedef Kokkos::StaticCrsGraph<ordinal_type, Kokkos::LayoutLeft, execution_space, memory_traits, size_type> staticcrsgraph_type;
#endif
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

  /// \name Storage of the actual sparsity structure and values.
  ///
  /// BlockCrsMatrix uses the compressed sparse row (CSR) storage format to
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
  /// Then BlockCrsMatrix needs methods to get these from the graph.
  BlockCrsMatrix () :
    numCols_ (0),
    blockDim_ (0)
  {}

  //! Copy constructor (shallow copy).
  template<typename SType,
           typename OType,
           class DType,
           class MTType,
           typename IType>
  BlockCrsMatrix (const BlockCrsMatrix<SType,OType,DType,MTType,IType> & B) :
    graph (B.graph.entries, B.graph.row_map),
    values (B.values),
    dev_config (B.dev_config),
    numCols_ (B.numCols ()),
    blockDim_ (B.blockDim ())
  {
    graph.row_block_offsets = B.graph.row_block_offsets;
    //MD: Changed the copy constructor of graph
    //as the constructor of StaticCrsGraph does not allow copy from non const version.
  }

  /// \brief Construct with a graph that will be shared.
  ///
  /// Allocate the values array for subsequent fill.
  BlockCrsMatrix (const std::string& arg_label,
                  const staticcrsgraph_type& arg_graph, 
                  const OrdinalType& blockDimIn) :
    graph (arg_graph),
    values (arg_label, arg_graph.entries.extent(0)),
    numCols_ (maximum_entry (arg_graph) + 1),
    blockDim_ (blockDimIn)
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
  /// The \c pad argument is currently not used.
  BlockCrsMatrix (const std::string &label,
                  OrdinalType nrows,
                  OrdinalType ncols,
                  size_type annz,
                  ScalarType* val,
                  OrdinalType* rows,
                  OrdinalType* cols,
                  OrdinalType blockdim,
                  bool pad = false)
  {
    (void) pad;
    ctor_impl (label, nrows, ncols, annz, val, rows, cols, blockdim);
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
  BlockCrsMatrix (const std::string& label,
                  const OrdinalType nrows,
                  const OrdinalType ncols,
                  const size_type annz,
                  const values_type& vals,
                  const row_map_type& rows,
                  const index_type& cols,
                  const OrdinalType blockDimIn) :
    graph (cols, rows),
    values (vals),
    numCols_ (ncols),
    blockDim_ (blockDimIn)
  {

    const ordinal_type actualNumRows = (rows.extent (0) != 0) ?
      static_cast<ordinal_type> (rows.extent (0) - static_cast<size_type> (1)) :
      static_cast<ordinal_type> (0);
    if (nrows != actualNumRows) {
      std::ostringstream os;
      os << "Input argument nrows = " << nrows << " != the actual number of "
        "rows " << actualNumRows << " according to the 'rows' input argument.";
      throw std::invalid_argument (os.str ());
    }
    // nnz returns graph.entries.extent(0) i.e. ptr[ nrows + 1 ] nnz entry
    // input annz is nnz of values, not comparable with block ptr 'nnz' i.e. numBlocks
    if (blockDim_ <= 0) {
      std::ostringstream os;
      os << "Input argument blockDim = " << blockDim_
         << " is not larger than 0.";
      throw std::invalid_argument (os.str ());
    }
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
  BlockCrsMatrix (const std::string& label,
                  const OrdinalType& ncols,
                  const values_type& vals,
                  const staticcrsgraph_type& graph_,
                  const OrdinalType& blockDimIn) :
    graph (graph_),
    values (vals),
    numCols_ (ncols),
    blockDim_ (blockDimIn)
  {}

  /// \brief Constructor that accepts a CrsMatrix and block dimension,
  ///        assuming the provided CrsMatrix has appropriate block structure.
  template<typename SType,
           typename OType,
           class DType,
           class MTType,
           typename IType>
  BlockCrsMatrix (const KokkosSparse::CrsMatrix<SType, OType, DType, MTType, IType> &crs_mtx,
                  const OrdinalType blockDimIn)
  {
    typedef typename KokkosSparse::CrsMatrix<SType, OType, DType, MTType, IType> crs_matrix_type;
    typedef typename crs_matrix_type::staticcrsgraph_type crs_graph_type;
    typedef typename crs_graph_type::entries_type crs_graph_entries_type;
    typedef typename crs_graph_type::row_map_type crs_graph_row_map_type;

    blockDim_ = blockDimIn;
    numCols_ = crs_mtx.numCols() / blockDim_;
    values = crs_mtx.values;

    OrdinalType nbrows = crs_mtx.numRows()/blockDim_; // actual number of block rows; add 1 for ptr length

    // block_rows will accumulate the number of blocks per row - this is NOT the row_map with cum sum!!
    std::vector<OrdinalType> block_rows( nbrows, 0 );

    typename crs_graph_row_map_type::HostMirror h_crs_row_map = Kokkos::create_mirror_view(crs_mtx.graph.row_map);
    Kokkos::deep_copy(h_crs_row_map, crs_mtx.graph.row_map);
    typename crs_graph_entries_type::HostMirror h_crs_entries = Kokkos::create_mirror_view(crs_mtx.graph.entries);
    Kokkos::deep_copy(h_crs_entries, crs_mtx.graph.entries);

    // determine size of block cols indices == number of blocks, i.e. nnz for the block CRS graph
    OrdinalType numBlocks = 0;
    for ( OrdinalType i = 0; i < crs_mtx.numRows(); i+=blockDim_ ) {
      numBlocks += ( h_crs_row_map(i+1) - h_crs_row_map(i) ) / blockDim_; // cum sum
      block_rows[ i/blockDim_ ] = ( h_crs_row_map(i+1) - h_crs_row_map(i) ) / blockDim_; // frequency counts
    }

    // create_staticcrsgraph takes the frequency of blocks per row
    // and returns the cum sum pointer row_map with nbrows+1 size, and total numBlocks in the final entry
    graph = Kokkos::create_staticcrsgraph<staticcrsgraph_type> ("blockgraph", block_rows);
    typename values_type::HostMirror h_values = Kokkos::create_mirror_view (values);
    typename index_type::HostMirror h_entries = Kokkos::create_mirror_view (graph.entries);

    for (OrdinalType i = 0; i < nbrows; ++i) {
      OrdinalType blks_in_row = block_rows[i];
      
      OrdinalType offset_into_blkcolidx_start = graph.row_map(i);
      OrdinalType offset_into_colidx_start = offset_into_blkcolidx_start*blockDim_*blockDim_;

      for ( OrdinalType lidx = 0; lidx < blks_in_row; ++lidx ) {
        h_entries( offset_into_blkcolidx_start+lidx ) = h_crs_entries( offset_into_colidx_start + blockDim_*lidx ) / blockDim_;
      }
    }

    Kokkos::deep_copy (graph.entries, h_entries);
  }


  /// Declaration for ctor_impl - this member function is not inlined
  void
  ctor_impl (const std::string &label,
          const OrdinalType nrows,
          const OrdinalType ncols,
          const size_type annz,
          ScalarType* val,
          OrdinalType* rows,
          OrdinalType* cols,
          const OrdinalType blockDimIn);


  /// \brief Given an array of blocks, sum the values into corresponding 
  ///        block in BlockCrsMatrix
  /// \param rowi [in]   is a block-row index
  /// \param ncol [in]   is number of blocks referenced in cols[] array
  /// \param cols[] [in] are block colidxs within the block-row to be summed into
  ///                    ncol entries
  /// \param vals[] [in] array containing 'block' of values
  ///        ncol*block_size*block_size entries
  ///        assume vals block is provided in 'LayoutRight' or 'Row Major' format, that is 
  ///        e.g. 2x2 block [ a b ; c d ] provided as flattened 1d array as [a b c d]
  ///        Assume that each block is stored contiguously in vals:
  ///        [a b; c d] [e f; g h] -> [a b c d e f g h]
  ///        If so, then i in [0, ncols) for cols[] 
  ///        maps to i*block_size*block_size in vals[]
  KOKKOS_INLINE_FUNCTION
  OrdinalType
  sumIntoValues (const OrdinalType rowi,
                 const OrdinalType cols[],
                 const OrdinalType ncol,
                 const ScalarType vals[],
                 const bool is_sorted = false,
                 const bool force_atomic = false) const
  {
    SparseBlockRowView<BlockCrsMatrix> row_view = this->block_row (rowi);
    const ordinal_type block_size = this->blockDim();

    ordinal_type numValid = 0; // number of valid local column indices

    for (ordinal_type i = 0; i < ncol; ++i) {

      // Find offset into values for block-row rowi and colidx cols[i]
      // cols[i] is the index to match
      // blk_offset is the offset for block colidx from bptr[rowi] to bptr[rowi + 1] (not global offset)
      // colidx_ and values_ are already offset to the beginning of blockrow rowi
      auto blk_offset = row_view.findRelBlockOffset(cols[i], is_sorted);
      if ( blk_offset != Kokkos::Details::ArithTraits<ordinal_type>::max() ) {
        ordinal_type offset_into_vals = i*block_size*block_size; //stride == 1 assumed between elements
        for ( ordinal_type lrow = 0; lrow < block_size; ++lrow ) {
          auto local_row_values = row_view.local_row_in_block(blk_offset, lrow); // pointer to start of specified local row within this block
          for ( ordinal_type lcol = 0; lcol < block_size; ++lcol ) {
            if (force_atomic) {
              Kokkos::atomic_add (&(local_row_values[lcol]), vals[ offset_into_vals + lrow*block_size + lcol ]);
            }
            else {
              local_row_values[lcol] += vals[ offset_into_vals + lrow*block_size + lcol];
            }
          }
        }
        ++numValid;
      }
    } // end for ncol
    return numValid;
  }


  /// \brief Given an array of blocks, replace the values of corresponding 
  ///        blocks in BlockCrsMatrix
  /// \param rowi [in]   is a block-row index
  /// \param ncol [in]   is number of blocks referenced in cols[] array
  /// \param cols[] [in] are block colidxs within the block-row to be summed into
  ///                    ncol entries
  /// \param vals[] [in] array containing 'block' of values
  //        ncol*block_size*block_size entries
  //        assume vals block is provided in 'LayoutRight' or 'Row Major' format, that is 
  //        e.g. 2x2 block [ a b ; c d ] provided as flattened 1d array as [a b c d]
  //        Assume that each block is stored contiguously in vals:
  //        [a b; c d] [e f; g h] -> [a b c d e f g h]
  //        If so, then i in [0, ncols) for cols[] 
  //        maps to i*block_size*block_size in vals[]
  KOKKOS_INLINE_FUNCTION
  OrdinalType
  replaceValues (const OrdinalType rowi,
                 const OrdinalType cols[],
                 const OrdinalType ncol,
                 const ScalarType vals[],
                 const bool is_sorted = false,
                 const bool force_atomic = false) const
  {
    SparseBlockRowView<BlockCrsMatrix> row_view = this->block_row (rowi);
    const ordinal_type block_size = this->blockDim();

    ordinal_type numValid = 0; // number of valid local column indices

    for (ordinal_type i = 0; i < ncol; ++i) {

      // Find offset into values for block-row rowi and colidx cols[i]
      // cols[i] is the index to match
      // blk_offset is the offset for block colidx from bptr[rowi] to bptr[rowi + 1] (not global offset)
      // colidx_ and values_ are already offset to the beginning of blockrow rowi
      auto blk_offset = row_view.findRelBlockOffset(cols[i], is_sorted);
      if ( blk_offset != Kokkos::Details::ArithTraits<ordinal_type>::max() ) {
        ordinal_type offset_into_vals = i*block_size*block_size; //stride == 1 assumed between elements
        for ( ordinal_type lrow = 0; lrow < block_size; ++lrow ) {
          auto local_row_values = row_view.local_row_in_block(blk_offset, lrow); // pointer to start of specified local row within this block
          for ( ordinal_type lcol = 0; lcol < block_size; ++lcol ) {
            if (force_atomic) {
              Kokkos::atomic_assign(&(local_row_values[lcol]), vals[ offset_into_vals + lrow*block_size + lcol ]);
            }
            else {
              local_row_values[lcol] = vals[ offset_into_vals + lrow*block_size + lcol];
            }
          }
        }
        ++numValid;
      }
    } // end for ncol
    return numValid;
  }

  //! Attempt to assign the input matrix to \c *this.
  // Are the CUDA sparse handles needed to be copied here??
  template<typename aScalarType, typename aOrdinalType, class aDevice, class aMemoryTraits,typename aSizeType>
  BlockCrsMatrix&
  operator= (const BlockCrsMatrix<aScalarType, aOrdinalType, aDevice, aMemoryTraits, aSizeType>& mtx)
  {
    numCols_ = mtx.numCols ();
    blockDim_ = mtx.blockDim ();
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

  //! The block dimension in the sparse block matrix.
  KOKKOS_INLINE_FUNCTION ordinal_type blockDim () const {
    return blockDim_ ;
  }

  //! The number of stored entries in the sparse matrix.
  KOKKOS_INLINE_FUNCTION size_type nnz () const {
    return graph.entries.extent (0);
  }

  friend struct SparseBlockRowView<BlockCrsMatrix>;

  /// \brief Return a SparseBlockRowView of block-row i of the matrix.
  ///
  /// If row i does not belong to the matrix, return an empty view.
  ///
  /// The returned object \c view implements the following interface:
  /// <ul>
  /// <li> \c view.length is the number of entries (i.e. blocks) 
  ///      in the block row </li>
  /// <li> \c view.local_row_in_block_row(K, i) returns a nonconst pointer
  ///      to the values of the ith local row in the k-th block of the block-row </li>
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
  SparseBlockRowView<BlockCrsMatrix> block_row (const ordinal_type i) const {

    const size_type start = graph.row_map(i); // total num blocks prior to this block-row
    const ordinal_type count = static_cast<ordinal_type> (graph.row_map(i+1) - start); // num blocks in this row

    if (count == 0) {
      return SparseBlockRowView<BlockCrsMatrix> (nullptr, nullptr, 1, 0);
    } else {
      return SparseBlockRowView<BlockCrsMatrix> (values, graph.entries, blockDim(), count, start);
    }
  }

  /// \brief Return a SparseBlockRowViewConst of block-row i of the matrix.
  ///
  /// If row i does not belong to the matrix, return an empty view.
  ///
  /// The returned object \c view implements the following interface:
  /// <ul>
  /// <li> \c view.length is the number of entries (i.e. blocks) 
  ///      in the block row </li>
  /// <li> \c view.local_row_in_block_row(K, i) returns a nonconst pointer
  ///      to the values of the ith local row in the k-th block of the block-row </li>
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
  SparseBlockRowViewConst<BlockCrsMatrix> block_row_Const (const ordinal_type i) const {

    const size_type start = graph.row_map(i); // total num blocks prior to this block-row
    const ordinal_type count = static_cast<ordinal_type> (graph.row_map(i+1) - start); // num blocks in this row

    if (count == 0) {
      return SparseBlockRowViewConst<BlockCrsMatrix> (nullptr, nullptr, 1, 0);
    } else {
      return SparseBlockRowViewConst<BlockCrsMatrix> (values, graph.entries, blockDim(), count, start);
    }
  }

private:
  ordinal_type numCols_;
  ordinal_type blockDim_; // TODO Assuming square blocks for now - add blockRowDim, blockColDim
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

// Input assumptions:
//   rows is pointer rep for the row_map member View of the BlockCrsMatrix graph (i.e. cum sum of number of blocks per block-row)
//   cols is pointer rep for the entries member View of the BlockCrsMatrix graph (colidx for block-row blocks)
//   annz is the total number of non-zeros in the CrsMatrix (equal to blockDim*blockDim*numBlocks)
template< typename ScalarType , typename OrdinalType, class Device, class MemoryTraits, typename SizeType >
void
BlockCrsMatrix<ScalarType , OrdinalType, Device, MemoryTraits, SizeType >::
ctor_impl (const std::string &label,
           const OrdinalType nrows,
           const OrdinalType ncols,
           const size_type annz,
           ScalarType* val,
           OrdinalType* rows,
           OrdinalType* cols,
           const OrdinalType blockDimIn)
{
  numCols_ = ncols;
  blockDim_ = blockDimIn;

  // Wrap the raw pointers in unmanaged host Views
  typename values_type::HostMirror unman_val( val, annz );
  typename row_map_type::HostMirror unman_rows( rows, nrows+1);
  typename index_type::HostMirror unman_cols( cols, ncols );

  // Create temporary Views for row_map and entries because the StaticCrsGraph ctor requires View inputs
  values_type tmp_row_map("tmp_row_map", nrows+1);
  values_type tmp_entries("tmp_entries", ncols);

  Kokkos::deep_copy( val, unman_val );
  Kokkos::deep_copy( tmp_row_map, unman_rows );
  Kokkos::deep_copy( tmp_entries, unman_cols );

  // Initialize graph using the temp entries and row_map Views
  graph = staticcrsgraph_type( tmp_entries, tmp_row_map );
}

}} // namespace KokkosSparse::Experimental
#endif
