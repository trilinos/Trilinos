// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_CRSMATRIX_HPP
#define STOKHOS_CRSMATRIX_HPP

#include <fstream>
#include <iomanip>

#include "Kokkos_Core.hpp"
#include "Kokkos_StaticCrsGraph.hpp"

#include "Stokhos_Multiply.hpp"
#include "Stokhos_MatrixMarket.hpp"

namespace Stokhos {

struct DeviceConfig {
  struct Dim3 {
    size_t x, y, z;
    Dim3(const size_t x_, const size_t y_ = 1, const size_t z_ = 1) :
      x(x_), y(y_), z(z_) {}
  };

  Dim3 block_dim;
  size_t num_blocks;
  size_t num_threads_per_block;

  DeviceConfig(const size_t num_blocks_,
               const size_t threads_per_block_x_,
               const size_t threads_per_block_y_ = 1,
               const size_t threads_per_block_z_ = 1) :
    block_dim(threads_per_block_x_,threads_per_block_y_,threads_per_block_z_),
    num_blocks(num_blocks_),
    num_threads_per_block(block_dim.x * block_dim.y * block_dim.z)
    {}
};

/** \brief  CRS matrix.  */
template <typename ValueType, typename Device,
          typename Layout = Kokkos::LayoutRight>
class CrsMatrix {
public:
  typedef Device execution_space;
  typedef ValueType value_type;
  typedef Kokkos::View< value_type[], Layout, execution_space > values_type;
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE // Don't remove this until Kokkos has removed the deprecated code path probably around September 2018
  typedef Kokkos::StaticCrsGraph< int , Layout, execution_space , int > graph_type;
#else
  typedef Kokkos::StaticCrsGraph< int , Layout, execution_space , void, int > graph_type;
#endif

  typedef CrsMatrix< ValueType, typename values_type::host_mirror_space, Layout> HostMirror;

  values_type values;
  graph_type graph;
  Stokhos::DeviceConfig dev_config;

  CrsMatrix() : dev_config(0, 0) {}
  CrsMatrix(Stokhos::DeviceConfig dev_config_) : dev_config(dev_config_) {}
};

// Generic matrix vector multiply kernel for CrsMatrix
template <typename MatrixValue,
          typename Layout,
          typename Device,
          typename InputVectorType,
          typename OutputVectorType>
class Multiply< CrsMatrix<MatrixValue,Device,Layout>,
                InputVectorType,
                OutputVectorType,
                void,
                IntegralRank<1> >
{
public:
  typedef CrsMatrix<MatrixValue,Device,Layout> matrix_type;
  typedef InputVectorType input_vector_type;
  typedef OutputVectorType output_vector_type;

  typedef Device execution_space;
  typedef typename execution_space::size_type size_type;
  typedef typename output_vector_type::value_type scalar_type;

  const matrix_type m_A;
  const input_vector_type m_x;
  output_vector_type m_y;

  Multiply( const matrix_type& A,
            const input_vector_type& x,
            output_vector_type& y )
  : m_A( A )
  , m_x( x )
  , m_y( y )
  {}

  //--------------------------------------------------------------------------

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type iRow ) const
  {
    const size_type iEntryBegin = m_A.graph.row_map[iRow];
    const size_type iEntryEnd   = m_A.graph.row_map[iRow+1];

    scalar_type sum = 0;

    for ( size_type iEntry = iEntryBegin; iEntry < iEntryEnd; ++iEntry ) {
      sum += m_A.values(iEntry) * m_x( m_A.graph.entries(iEntry) );
    }

    m_y(iRow) = sum;
  }

  static void apply( const matrix_type & A,
                     const input_vector_type & x,
                     output_vector_type & y )
  {
    const size_t row_count = A.graph.row_map.extent(0) - 1;
    Kokkos::parallel_for( row_count, Multiply(A,x,y) );
  }
};

// Generic matrix multi-vector multiply kernel for CrsMatrix
template <typename MatrixValue,
          typename Layout,
          typename Device,
          typename InputMultiVectorType,
          typename OutputMultiVectorType,
          typename OrdinalType >
class Multiply< CrsMatrix<MatrixValue,Device,Layout>,
                InputMultiVectorType,
                OutputMultiVectorType,
                std::vector<OrdinalType>,
                IntegralRank<2> >
{
public:
  typedef CrsMatrix<MatrixValue,Device,Layout> matrix_type;
  typedef InputMultiVectorType input_multi_vector_type;
  typedef OutputMultiVectorType output_multi_vector_type;
  typedef std::vector<OrdinalType> column_indices_type;

  typedef Device execution_space;
  typedef typename execution_space::size_type size_type;
  typedef typename output_multi_vector_type::value_type scalar_type;

  const matrix_type m_A;
  const input_multi_vector_type m_x;
  output_multi_vector_type m_y;
  const column_indices_type m_col_indices;
  const size_type m_num_vecs;

  Multiply( const matrix_type& A,
            const input_multi_vector_type& x,
            output_multi_vector_type& y,
            const column_indices_type& col_indices )
  : m_A( A )
  , m_x( x )
  , m_y( y )
  , m_col_indices( col_indices )
  , m_num_vecs( col_indices.size() )
  {}

  //--------------------------------------------------------------------------

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type iRow ) const
  {
    const size_type iEntryBegin = m_A.graph.row_map[iRow];
    const size_type iEntryEnd   = m_A.graph.row_map[iRow+1];

    for (size_type j=0; j<m_num_vecs; j++) {
      size_type iCol = m_col_indices[j];

      scalar_type sum = 0.0;

      for ( size_type iEntry = iEntryBegin ; iEntry < iEntryEnd ; ++iEntry ) {
        sum += m_A.values(iEntry) * m_x(  m_A.graph.entries(iEntry), iCol );
      }

      m_y( iRow, iCol ) = sum;

    }

  }

  static void apply( const matrix_type& A,
                     const input_multi_vector_type& x,
                     output_multi_vector_type& y,
                     const column_indices_type& col )
  {
    const size_t n = A.graph.row_map.extent(0) - 1 ;
    //Kokkos::parallel_for( n , Multiply(A,x,y,col) );

    const size_t block_size = 20;
    const size_t num_vecs = col.size();
    std::vector<OrdinalType> block_col;
    block_col.reserve(block_size);
    for (size_t block=0; block<num_vecs; block+=block_size) {
      const size_t bs =
        block+block_size <= num_vecs ? block_size : num_vecs-block;
      block_col.resize(bs);
      for (size_t i=0; i<bs; ++i)
        block_col[i] = col[block+i];
      Kokkos::parallel_for( n , Multiply(A,x,y,block_col) );
    }
  }
};

#define USE_NEW 1
#if USE_NEW
// Generic matrix multi-vector multiply kernel for CrsMatrix
// Experimenting with blocking of column and row loops to improve cache
// performance.  Seems to help signficantly on SandyBridge, little difference
// on MIC (although not extensive investigation of block sizes).
template <typename MatrixValue,
          typename Layout,
          typename Device,
          typename InputMultiVectorType,
          typename OutputMultiVectorType >
class Multiply< CrsMatrix<MatrixValue,Device,Layout>,
                InputMultiVectorType,
                OutputMultiVectorType,
                void,
                IntegralRank<2> >
{
public:
  typedef CrsMatrix<MatrixValue,Device,Layout> matrix_type;
  typedef InputMultiVectorType input_multi_vector_type;
  typedef OutputMultiVectorType output_multi_vector_type;

  typedef Device execution_space;
  typedef typename execution_space::size_type size_type;
  typedef typename output_multi_vector_type::value_type scalar_type;

  const matrix_type m_A;
  const input_multi_vector_type m_x;
  output_multi_vector_type m_y;
  const size_type m_num_row;
  const size_type m_num_col;

  static const size_type m_block_row_size = 32;
  static const size_type m_block_col_size = 20;

  Multiply( const matrix_type& A,
            const input_multi_vector_type& x,
            output_multi_vector_type& y )
  : m_A( A )
  , m_x( x )
  , m_y( y )
  , m_num_row( A.graph.row_map.extent(0)-1 )
  , m_num_col( m_y.extent(1) )
  {
  }

  //--------------------------------------------------------------------------

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type iBlockRow ) const
  {
    // Number of rows in this block
    const size_type num_row =
      iBlockRow+m_block_row_size <= m_num_row ?
      m_block_row_size : m_num_row-iBlockRow;

    // Loop over block columns of x
    for (size_type iBlockCol=0; iBlockCol<m_num_col; iBlockCol+=m_block_col_size) {
      // Number of columns in this block
      const size_type num_col =
        iBlockCol+m_block_col_size <= m_num_col ?
        m_block_col_size : m_num_col-iBlockCol;

      // Loop over rows in this block of A
      const size_type iRowEnd = iBlockRow + num_row;
      for (size_type iRow=iBlockRow; iRow<iRowEnd; ++iRow) {

        // Range of column entries for this row
        const size_type iEntryBegin = m_A.graph.row_map[iRow];
        const size_type iEntryEnd   = m_A.graph.row_map[iRow+1];

        // Loop over columns in this block of x
        const size_type iColEnd = iBlockCol + num_col;
        for (size_type iCol=iBlockCol; iCol<iColEnd; iCol++) {

          // Loop columns of A for this row
          scalar_type sum = 0.0;
          for (size_type iEntry = iEntryBegin; iEntry<iEntryEnd; ++iEntry) {
            sum += m_A.values(iEntry) * m_x(  m_A.graph.entries(iEntry), iCol );
          }
          m_y( iRow, iCol ) = sum;

        }

      }

    }

  }

  static void apply( const matrix_type & A,
                     const input_multi_vector_type& x,
                     output_multi_vector_type& y )
  {
    // Parallelize over row blocks of size m_block_row_size
    const size_type num_row = A.graph.row_map.extent(0) - 1;
    const size_type n = (num_row+m_block_row_size-1) / m_block_row_size;
    Kokkos::parallel_for( n , Multiply(A,x,y) );
  }
};
#else
// Generic matrix multi-vector multiply kernel for CrsMatrix
template <typename MatrixValue,
          typename Layout,
          typename Device,
          typename InputMultiVectorType,
          typename OutputMultiVectorType >
class Multiply< CrsMatrix<MatrixValue,Device,Layout>,
                InputMultiVectorType,
                OutputMultiVectorType,
                void,
                IntegralRank<2> >
{
public:
  typedef CrsMatrix<MatrixValue,Device,Layout> matrix_type;
  typedef InputMultiVectorType input_multi_vector_type;
  typedef OutputMultiVectorType output_multi_vector_type;

  typedef Device execution_space;
  typedef typename execution_space::size_type size_type;
  typedef typename output_multi_vector_type::value_type scalar_type;

  const matrix_type m_A;
  const input_multi_vector_type m_x;
  output_multi_vector_type m_y;
  const size_type m_num_vecs;

  Multiply( const matrix_type& A,
            const input_multi_vector_type& x,
            output_multi_vector_type& y)
  : m_A( A )
  , m_x( x )
  , m_y( y )
  , m_num_vecs( m_y.extent(1) )
  {}

  //--------------------------------------------------------------------------

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type iRow ) const
  {
    const size_type iEntryBegin = m_A.graph.row_map[iRow];
    const size_type iEntryEnd   = m_A.graph.row_map[iRow+1];

    for (size_type iCol=0; iCol<m_num_vecs; iCol++) {

      scalar_type sum = 0.0;

      for ( size_type iEntry = iEntryBegin ; iEntry < iEntryEnd ; ++iEntry ) {
        sum += m_A.values(iEntry) * m_x(  m_A.graph.entries(iEntry), iCol );
      }

      m_y( iRow, iCol ) = sum;

    }

  }

  static void apply( const matrix_type& A,
                     const input_multi_vector_type& x,
                     output_multi_vector_type& y )
  {
    const size_t n = A.graph.row_map.extent(0) - 1 ;
    Kokkos::parallel_for( n , Multiply(A,x,y) );

    // const size_t block_size = 20;
    // const size_t num_vecs = col.size();
    // std::vector<OrdinalType> block_col;
    // block_col.reserve(block_size);
    // for (size_t block=0; block<num_vecs; block+=block_size) {
    //   const size_t bs =
    //     block+block_size <= num_vecs ? block_size : num_vecs-block;
    //   block_col.resize(bs);
    //   for (size_t i=0; i<bs; ++i)
    //     block_col[i] = col[block+i];
    //   Kokkos::parallel_for( n , Multiply(A,x,y,block_col) );
    // }
  }
};
#endif

#if USE_NEW
// Generic matrix multi-vector multiply kernel for CrsMatrix
// Experimenting with blocking of column and row loops to improve cache
// performance.  Seems to help signficantly on SandyBridge, little difference
// on MIC (although not extensive investigation of block sizes).
template <typename MatrixValue,
          typename Layout,
          typename Device,
          typename InputViewType,
          typename OutputViewType>
class Multiply< CrsMatrix<MatrixValue,Device,Layout>,
                std::vector<InputViewType>,
                std::vector<OutputViewType>,
                void,
                IntegralRank<1> >
{
public:
  typedef CrsMatrix<MatrixValue,Device,Layout> matrix_type;
  typedef std::vector<InputViewType> input_multi_vector_type;
  typedef std::vector<OutputViewType> output_multi_vector_type;

  typedef Device execution_space;
  typedef typename execution_space::size_type size_type;
  typedef typename OutputViewType::value_type scalar_type;

  const matrix_type m_A;
  const input_multi_vector_type m_x;
  output_multi_vector_type m_y;
  const size_type m_num_row;
  const size_type m_num_col;

  static const size_type m_block_row_size = 32;
  static const size_type m_block_col_size = 20;

  Multiply( const matrix_type& A,
            const input_multi_vector_type& x,
            output_multi_vector_type& y )
  : m_A( A )
  , m_x( x )
  , m_y( y )
  , m_num_row( A.graph.row_map.extent(0)-1 )
  , m_num_col( x.size() )
  {
  }

  //--------------------------------------------------------------------------

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type iBlockRow ) const
  {
    // Number of rows in this block
    const size_type num_row =
      iBlockRow+m_block_row_size <= m_num_row ?
      m_block_row_size : m_num_row-iBlockRow;

    // Loop over block columns of x
    for (size_type iBlockCol=0; iBlockCol<m_num_col; iBlockCol+=m_block_col_size) {
      // Number of columns in this block
      const size_type num_col =
        iBlockCol+m_block_col_size <= m_num_col ?
        m_block_col_size : m_num_col-iBlockCol;

      // Loop over rows in this block of A
      const size_type iRowEnd = iBlockRow + num_row;
      for (size_type iRow=iBlockRow; iRow<iRowEnd; ++iRow) {

        // Range of column entries for this row
        const size_type iEntryBegin = m_A.graph.row_map[iRow];
        const size_type iEntryEnd   = m_A.graph.row_map[iRow+1];

        // Loop over columns in this block of x
        const size_type iColEnd = iBlockCol + num_col;
        for (size_type iCol=iBlockCol; iCol<iColEnd; iCol++) {

          // Loop columns of A for this row
          scalar_type sum = 0.0;
          for (size_type iEntry = iEntryBegin; iEntry<iEntryEnd; ++iEntry) {
            sum += m_A.values(iEntry) * m_x[iCol](m_A.graph.entries(iEntry));
          }
          m_y[iCol](iRow) = sum;

        }

      }

    }

  }

  static void apply( const matrix_type & A,
                     const input_multi_vector_type& x,
                     output_multi_vector_type& y )
  {
    // Parallelize over row blocks of size m_block_row_size
    const size_type num_row = A.graph.row_map.extent(0) - 1;
    const size_type n = (num_row+m_block_row_size-1) / m_block_row_size;
    Kokkos::parallel_for( n , Multiply(A,x,y) );
  }
};
#else
// Generic matrix multi-vector multiply kernel for CrsMatrix
template <typename MatrixValue,
          typename Layout,
          typename Device,
          typename InputViewType,
          typename OutputViewType>
class Multiply< CrsMatrix<MatrixValue,Device,Layout>,
                std::vector<InputViewType>,
                std::vector<OutputViewType>,
                void,
                IntegralRank<1> >
{
public:
  typedef CrsMatrix<MatrixValue,Device,Layout> matrix_type;
  typedef std::vector<InputViewType> input_multi_vector_type;
  typedef std::vector<OutputViewType> output_multi_vector_type;

  typedef Device execution_space;
  typedef typename execution_space::size_type size_type;
  typedef typename OutputViewType::value_type scalar_type;

  const matrix_type m_A;
  const input_multi_vector_type m_x;
  output_multi_vector_type m_y;
  const size_type m_num_vecs;

  Multiply( const matrix_type& A,
            const input_multi_vector_type& x,
            output_multi_vector_type& y )
  : m_A( A )
  , m_x( x )
  , m_y( y )
  , m_num_vecs( x.size() )
  {
  }

  //--------------------------------------------------------------------------

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type iRow ) const
  {
    const size_type iEntryBegin = m_A.graph.row_map[iRow];
    const size_type iEntryEnd   = m_A.graph.row_map[iRow+1];

    for (size_type iCol=0; iCol<m_num_vecs; iCol++) {

      scalar_type sum = 0.0;

      for ( size_type iEntry = iEntryBegin ; iEntry < iEntryEnd ; ++iEntry ) {
        sum += m_A.values(iEntry) * m_x[iCol](  m_A.graph.entries(iEntry) );
      }

      m_y[iCol]( iRow) = sum;

    }

  }

  static void apply( const matrix_type & A,
                     const input_multi_vector_type& x,
                     output_multi_vector_type& y )
  {
    const size_t n = A.graph.row_map.extent(0) - 1 ;
    Kokkos::parallel_for( n , Multiply(A,x,y) );

    // const size_t block_size = 20;
    // const size_t num_vecs = x.size();
    // input_multi_vector_type xx;
    // output_multi_vector_type yy;
    // xx.reserve(block_size);
    // yy.reserve(block_size);
    // for (size_t block=0; block<num_vecs; block+=block_size) {
    //   const size_t bs =
    //     block+block_size <= num_vecs ? block_size : num_vecs-block;
    //   xx.resize(bs);
    //   yy.resize(bs);
    //   for (size_t i=0; i<bs; ++i) {
    //     xx[i] = x[block+i];
    //     yy[i] = y[block+i];
    //   }
    //   Kokkos::parallel_for( n , Multiply(A,xx,yy) );
    // }
  }
};
#endif

// Matrix multivector multiply specializations for one column at a time
class SingleColumnMultivectorMultiply {};
template <typename MatrixValue,
          typename Layout,
          typename Device,
          typename InputMultiVectorType,
          typename OutputMultiVectorType,
          typename OrdinalType>
void multiply(const CrsMatrix<MatrixValue,Device,Layout>& A,
              const InputMultiVectorType& x,
              OutputMultiVectorType& y,
              const std::vector<OrdinalType>& col_indices,
              SingleColumnMultivectorMultiply)
{
  typedef CrsMatrix<MatrixValue,Device,Layout> MatrixType;

  typedef Kokkos::View<typename InputMultiVectorType::value_type*, typename InputMultiVectorType::array_layout, Device, Kokkos::MemoryUnmanaged> InputVectorType;
  typedef Kokkos::View<typename OutputMultiVectorType::value_type*, typename OutputMultiVectorType::array_layout, Device, Kokkos::MemoryUnmanaged> OutputVectorType;
  typedef Multiply<MatrixType,InputVectorType,OutputVectorType> multiply_type;
  for (size_t i=0; i<col_indices.size(); ++i) {
    InputVectorType x_view =
      Kokkos::subview( x , Kokkos::ALL() , col_indices[i] );
    OutputVectorType y_view =
      Kokkos::subview( y , Kokkos::ALL() , col_indices[i] );
    multiply_type::apply( A , x_view , y_view );
  }
}

template <typename MatrixValue,
          typename Layout,
          typename Device,
          typename InputVectorType,
          typename OutputVectorType>
void multiply(const CrsMatrix<MatrixValue,Device,Layout>& A,
              const std::vector<InputVectorType>& x,
              std::vector<OutputVectorType>& y,
              SingleColumnMultivectorMultiply)
{
  typedef CrsMatrix<MatrixValue,Device,Layout> MatrixType;
  typedef Multiply<MatrixType,InputVectorType,OutputVectorType> multiply_type;
  for (size_t i=0; i<x.size(); ++i) {
    multiply_type::apply( A , x[i] , y[i] );
  }
}

} // namespace Stokhos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

template <typename ValueType, typename Layout, typename Device>
typename Stokhos::CrsMatrix<ValueType,Device,Layout>::HostMirror
create_mirror(const Stokhos::CrsMatrix<ValueType,Device,Layout>& A) {
  typename Stokhos::CrsMatrix<ValueType,Device,Layout>::HostMirror mirror_A;
  mirror_A.values = Kokkos::create_mirror(A.values);
  mirror_A.graph = Kokkos::create_mirror(A.graph); // this deep copies
  mirror_A.dev_config = A.dev_config;
  return mirror_A;
}

template <typename ValueType, typename Layout, typename Device>
typename Stokhos::CrsMatrix<ValueType,Device,Layout>::HostMirror
create_mirror_view(const Stokhos::CrsMatrix<ValueType,Device,Layout>& A) {
  typename Stokhos::CrsMatrix<ValueType,Device,Layout>::HostMirror mirror_A;
  mirror_A.values = Kokkos::create_mirror_view(A.values);
  mirror_A.graph = Kokkos::create_mirror(A.graph); // this deep copies
  mirror_A.dev_config = A.dev_config;
  return mirror_A;
}

template <typename ValueType, typename Layout, typename DstDevice,
          typename SrcDevice>
void
deep_copy(const Stokhos::CrsMatrix<ValueType,DstDevice,Layout>& dst,
          const Stokhos::CrsMatrix<ValueType,SrcDevice,Layout>& src) {
  Kokkos::deep_copy(dst.values, src.values);
}

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Stokhos {

// MatrixMarket writer for CrsMatrix
template < typename MatrixValue, typename Layout, typename Device >
class MatrixMarketWriter< CrsMatrix<MatrixValue,Device,Layout> >
{
public:
  typedef CrsMatrix<MatrixValue,Device,Layout> matrix_type ;
  typedef Device execution_space ;
  typedef typename execution_space::size_type size_type ;

  static void write(const matrix_type& A, const std::string& filename) {
    std::ofstream file(filename.c_str());
    file.precision(16);
    file.setf(std::ios::scientific);

    typename matrix_type::HostMirror hA = Kokkos::create_mirror_view(A);
    Kokkos::deep_copy(hA, A);

    const size_type nRow = hA.graph.row_map.extent(0) - 1 ;

    // Write banner
    file << "%%MatrixMarket matrix coordinate real general" << std::endl;
    file << nRow << " " << nRow << " " << hA.values.extent(0) << std::endl;

    for (size_type row=0; row<nRow; ++row) {
      size_type entryBegin = hA.graph.row_map(row);
      size_type entryEnd = hA.graph.row_map(row+1);
      for (size_type entry=entryBegin; entry<entryEnd; ++entry) {
        file << row+1 << " " << hA.graph.entries(entry)+1 << " "
             << std::setw(22) << hA.values(entry) << std::endl;
      }
    }

    file.close();
  }
};

} // namespace Stokhos

#endif /* #ifndef STOKHOS_CRSMATRIX_HPP */
