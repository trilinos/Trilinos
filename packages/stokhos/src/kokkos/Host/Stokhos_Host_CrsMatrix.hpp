// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_HOST_CRSMATRIX_HPP
#define STOKHOS_HOST_CRSMATRIX_HPP

#include <fstream>
#include <iomanip>

#include "KokkosArray_Host.hpp"
#include "KokkosArray_ParallelFor.hpp"

#include "Stokhos_Multiply.hpp"
#include "Stokhos_CrsMatrix.hpp"

#include "Stokhos_ConfigDefs.h"
#ifdef HAVE_STOKHOS_MKL
#include "mkl.h"
#endif

namespace Stokhos {

template< typename MatrixValue , typename VectorValue >
class Multiply<
  CrsMatrix< MatrixValue , KokkosArray::Host > ,
  KokkosArray::View< VectorValue[] , KokkosArray::Host > ,
  KokkosArray::View< VectorValue[] , KokkosArray::Host > ,
  DefaultSparseMatOps >
{
public:
  typedef KokkosArray::Host                         device_type ;
  typedef device_type::size_type                    size_type ;
  typedef KokkosArray::View< VectorValue[] , device_type >  vector_type ;
  typedef CrsMatrix< MatrixValue , device_type >    matrix_type ;

  const matrix_type  m_A ;
  const vector_type  m_x ;
  const vector_type  m_y ;

  Multiply( const matrix_type & A ,
            const vector_type & x ,
            const vector_type & y )
  : m_A( A )
  , m_x( x )
  , m_y( y )
  {}

  //--------------------------------------------------------------------------

  inline
  void operator()( const size_type iRow ) const
  {
    const size_type iEntryBegin = m_A.graph.row_map[iRow];
    const size_type iEntryEnd   = m_A.graph.row_map[iRow+1];

    double sum = 0 ;

    for ( size_type iEntry = iEntryBegin ; iEntry < iEntryEnd ; ++iEntry ) {
      sum += m_A.values(iEntry) * m_x( m_A.graph.entries(iEntry) );
    }

    m_y(iRow) = sum ;
  }

  static void apply( const matrix_type & A ,
                     const vector_type & x ,
                     const vector_type & y )
  {
    const size_t row_count = A.graph.row_map.dimension(0) - 1 ;
    KokkosArray::parallel_for( row_count , Multiply(A,x,y) );
  }
};

template< typename MatrixValue , typename VectorValue >
class MMultiply<
  CrsMatrix< MatrixValue , KokkosArray::Host > ,
  KokkosArray::View< VectorValue** , KokkosArray::LayoutLeft, KokkosArray::Host > ,
  KokkosArray::View< VectorValue** , KokkosArray::LayoutLeft, KokkosArray::Host > ,
  DefaultSparseMatOps >
{
public:
  typedef KokkosArray::Host                                device_type ;
  typedef device_type::size_type                           size_type ;
  typedef KokkosArray::View< VectorValue** , KokkosArray::LayoutLeft, device_type >  multi_vector_type ;
  typedef CrsMatrix< MatrixValue , device_type >           matrix_type ;
  typedef VectorValue                                      value_type ;
  typedef int                                              Ordinal ;

  const matrix_type  m_A ;
  const multi_vector_type m_x ;
  const multi_vector_type m_y ;
  const std::vector<Ordinal> m_col_indices ;
  const size_t num_vecs ;

  MMultiply( const matrix_type & A ,
	     const multi_vector_type & x ,
	     const multi_vector_type & y ,
	     const std::vector<Ordinal> & col_indices )
  : m_A( A )
  , m_x( x )
  , m_y( y )
  , m_col_indices( col_indices )
  , num_vecs( col_indices.size() )
  {
    //std::cout << "num_vecs = " << num_vecs << std::endl;
  }

  //--------------------------------------------------------------------------

  inline
  void operator()( const size_type iRow ) const
  {
    const size_type iEntryBegin = m_A.graph.row_map[iRow];
    const size_type iEntryEnd   = m_A.graph.row_map[iRow+1];
    const size_t n = m_A.graph.row_map.dimension(0) - 1 ;

    for (size_t j=0; j<num_vecs; j++) {
      Ordinal col = m_col_indices[j];
      
      value_type y_tmp = 0.0;

      for ( size_type iEntry = iEntryBegin ; iEntry < iEntryEnd ; ++iEntry ) {
	y_tmp += m_A.values(iEntry) * m_x(  m_A.graph.entries(iEntry), col );
      }

      m_y( iRow, col ) = y_tmp;
    
    }

  }

  static void apply( const matrix_type & A ,
                     const multi_vector_type & x ,
                     const multi_vector_type & y ,
		     const std::vector<Ordinal> & col)
  {
    const size_t n = A.graph.row_map.dimension(0) - 1 ;
    const size_t block_size = 20;
    const size_t num_vecs = col.size();
    const size_t num_blocks = num_vecs / block_size;
    const size_t rem = num_vecs - num_blocks*block_size;
    const size_t bs = block_size;
    std::vector<Ordinal> block_col(block_size);
    for (size_t block=0; block<num_blocks; ++block) {
      for (size_t i=0; i<bs; ++i)
	block_col[i] = col[block*block_size+i];
      KokkosArray::parallel_for( n , MMultiply(A,x,y,block_col) );
    }
    if (rem > 0) {
      block_col.resize(rem);
      for (size_t i=0; i<block_size; ++i)
    	block_col[i] = col[num_blocks*block_size+i];
      KokkosArray::parallel_for( n , MMultiply(A,x,y,block_col) );
    }
  }
};

template< typename MatrixValue , typename VectorValue >
class MMultiply<
  CrsMatrix< MatrixValue , KokkosArray::Host > ,
  KokkosArray::View< VectorValue[] , KokkosArray::Host > ,
  KokkosArray::View< VectorValue[] , KokkosArray::Host > ,
  DefaultSparseMatOps >
{
public:
  typedef KokkosArray::Host                                device_type ;
  typedef device_type::size_type                           size_type ;
  typedef KokkosArray::View< VectorValue[] , device_type > vector_type ;
  typedef CrsMatrix< MatrixValue , device_type >           matrix_type ;
  typedef VectorValue                                      value_type ;

  const matrix_type  m_A ;
  const std::vector<vector_type> m_x ;
  const std::vector<vector_type> m_y ;

  MMultiply( const matrix_type & A ,
	     const std::vector<vector_type> & x ,
	     const std::vector<vector_type> & y )
  : m_A( A )
  , m_x( x )
  , m_y( y )
  {
  }

  //--------------------------------------------------------------------------

  inline
  void operator()( const size_type iRow ) const
  {
    const size_type iEntryBegin = m_A.graph.row_map[iRow];
    const size_type iEntryEnd   = m_A.graph.row_map[iRow+1];
    //const size_t n = m_A.graph.row_map.dimension(0) - 1 ;
    const size_t num_vecs = m_x.size();

    for (size_t j=0; j<num_vecs; j++) {
      
      value_type y_tmp = 0.0;

      for ( size_type iEntry = iEntryBegin ; iEntry < iEntryEnd ; ++iEntry ) {
	y_tmp += m_A.values(iEntry) * m_x[j](  m_A.graph.entries(iEntry) );
      }

      m_y[j]( iRow) = y_tmp;
    
    }

  }

  static void apply( const matrix_type & A ,
                     const std::vector<vector_type> & x ,
                     const std::vector<vector_type> & y )
  {
    const size_t n = A.graph.row_map.dimension(0) - 1 ;
    KokkosArray::parallel_for( n , MMultiply(A,x,y) );
  }
};

#ifdef HAVE_STOKHOS_MKL

class MKLSparseMatOps {}; 

template<>
class Multiply<
  CrsMatrix< double , KokkosArray::Host > ,
  KokkosArray::View< double[] , KokkosArray::Host > ,
  KokkosArray::View< double[] , KokkosArray::Host > ,
  MKLSparseMatOps >
{
public:
  typedef KokkosArray::Host                         device_type ;
  typedef device_type::size_type                    size_type ;
  typedef KokkosArray::View< double[] , device_type >  vector_type ;
  typedef CrsMatrix< double , device_type >    matrix_type ;

  static void apply( const matrix_type & A ,
                     const vector_type & x ,
                     const vector_type & y )
  {
    MKL_INT n = A.graph.row_map.dimension(0) - 1 ;
    double *A_values = A.values.ptr_on_device() ;
    MKL_INT *col_indices = A.graph.entries.ptr_on_device() ;
    MKL_INT *row_beg = const_cast<MKL_INT*>(A.graph.row_map.ptr_on_device()) ;
    MKL_INT *row_end = row_beg+1;
    char matdescra[6] = { 'G', 'x', 'N', 'C', 'x', 'x' };
    char trans = 'N';
    double alpha = 1.0;
    double beta = 0.0;

    double *x_values = x.ptr_on_device() ;
    double *y_values = y.ptr_on_device() ;
    
    mkl_dcsrmv(&trans, &n, &n, &alpha, matdescra, A_values, col_indices,
	       row_beg, row_end, x_values, &beta, y_values);
  }
};

namespace Impl {

  template <typename ValueType, typename OrdinalType, typename DeviceType>
  struct GatherTranspose {
    typedef ValueType value_type;
    typedef OrdinalType ordinal_type;
    typedef DeviceType device_type;
    typedef typename device_type::size_type size_type;
    typedef KokkosArray::View< value_type** , KokkosArray::LayoutLeft, device_type >  multi_vector_type ;
    typedef KokkosArray::View< value_type** , KokkosArray::LayoutLeft, device_type >  trans_multi_vector_type ;

    const multi_vector_type m_x;
    const trans_multi_vector_type m_xt;
    const std::vector<ordinal_type> m_indices;
    const size_t ncol;
    GatherTranspose(const multi_vector_type & x,
		    const trans_multi_vector_type& xt,
		    const std::vector<ordinal_type> & indices) :
      m_x(x), m_xt(xt), m_indices(indices), ncol(indices.size()) {}

    inline void operator()( const size_type row ) const {
      for (size_t col=0; col<ncol; ++col)
	m_xt(col,row) = m_x(row,m_indices[col]);
    }
    
    static void apply(const multi_vector_type & x,
		      const trans_multi_vector_type& xt,
		      const std::vector<ordinal_type> & indices) {
      const size_t n = xt.dimension(1);
      KokkosArray::parallel_for( n, GatherTranspose(x,xt,indices) );
    }
  };

  template <typename ValueType, typename OrdinalType, typename DeviceType>
  struct ScatterTranspose {
    typedef ValueType value_type;
    typedef OrdinalType ordinal_type;
    typedef DeviceType device_type;
    typedef typename device_type::size_type size_type;
    typedef KokkosArray::View< value_type** , KokkosArray::LayoutLeft, device_type >  multi_vector_type ;
    typedef KokkosArray::View< value_type** , KokkosArray::LayoutLeft, device_type >  trans_multi_vector_type ;

    const multi_vector_type m_x;
    const trans_multi_vector_type m_xt;
    const std::vector<ordinal_type> m_indices;
    const size_t ncol;
    ScatterTranspose(const multi_vector_type & x,
		     const trans_multi_vector_type& xt,
		     const std::vector<ordinal_type> & indices) :
      m_x(x), m_xt(xt), m_indices(indices), ncol(indices.size()) {}

    inline void operator()( const size_type row ) const {
      for (size_t col=0; col<ncol; ++col)
	m_x(row,m_indices[col]) = m_xt(col,row);
    }
    
    static void apply(const multi_vector_type & x,
		      const trans_multi_vector_type& xt,
		      const std::vector<ordinal_type> & indices) {
      const size_t n = xt.dimension(1);
      KokkosArray::parallel_for( n, ScatterTranspose(x,xt,indices) );
    }
  };

  template <typename ValueType, typename DeviceType>
  struct GatherVecTranspose {
    typedef ValueType value_type;
    typedef DeviceType device_type;
    typedef typename device_type::size_type size_type;
    typedef KokkosArray::View< value_type[] , device_type > vector_type ;
    typedef KokkosArray::View< value_type** , KokkosArray::LayoutLeft, device_type >  trans_multi_vector_type ;

    const std::vector<vector_type> m_x;
    const trans_multi_vector_type m_xt;
    const size_t ncol;
    GatherVecTranspose(const std::vector<vector_type> & x,
		       const trans_multi_vector_type& xt) :
      m_x(x), m_xt(xt), ncol(x.size()) {}

    inline void operator()( const size_type row ) const {
      for (size_t col=0; col<ncol; ++col)
	m_xt(col,row) = m_x[col](row);
    }
    
    static void apply(const std::vector<vector_type> & x,
		      const trans_multi_vector_type& xt) {
      const size_t n = xt.dimension(1);
      KokkosArray::parallel_for( n, GatherVecTranspose(x,xt) );
    }
  };

  template <typename ValueType, typename DeviceType>
  struct ScatterVecTranspose {
    typedef ValueType value_type;
    typedef DeviceType device_type;
    typedef typename device_type::size_type size_type;
    typedef KokkosArray::View< value_type[] , device_type > vector_type ;
    typedef KokkosArray::View< value_type** , KokkosArray::LayoutLeft, device_type >  trans_multi_vector_type ;

    const std::vector<vector_type> m_x;
    const trans_multi_vector_type m_xt;
    const size_t ncol;
    ScatterVecTranspose(const std::vector<vector_type> & x,
			const trans_multi_vector_type& xt) :
      m_x(x), m_xt(xt), ncol(x.size()) {}

    inline void operator()( const size_type row ) const {
      for (size_t col=0; col<ncol; ++col)
	m_x[col](row) = m_xt(col,row);
    }
    
    static void apply(const std::vector<vector_type> & x,
		      const trans_multi_vector_type& xt) {
      const size_t n = xt.dimension(1);
      KokkosArray::parallel_for( n, ScatterVecTranspose(x,xt) );
    }
  };

} // namespace Impl

template<>
class MMultiply<
  CrsMatrix< double , KokkosArray::Host > ,
  KokkosArray::View< double** , KokkosArray::LayoutLeft, KokkosArray::Host > ,
  KokkosArray::View< double** , KokkosArray::LayoutLeft, KokkosArray::Host > ,
  MKLSparseMatOps >
{
public:
  typedef KokkosArray::Host device_type ;
  typedef device_type::size_type size_type ;
  typedef KokkosArray::View< double** , KokkosArray::LayoutLeft, device_type >  multi_vector_type ;
  typedef KokkosArray::View< double** , KokkosArray::LayoutLeft, device_type >  trans_multi_vector_type ;
  typedef KokkosArray::View< double[] , device_type > vector_type ;
  typedef CrsMatrix< double , device_type > matrix_type ;
  typedef double value_type ;
  typedef int ordinal_type ;

public:

  static void apply( const matrix_type & A ,
                     const multi_vector_type & x ,
                     const multi_vector_type & y ,
		     const std::vector<ordinal_type> & indices)
  {
    MKL_INT n = A.graph.row_map.dimension(0) - 1 ;
    double *A_values = A.values.ptr_on_device() ;
    MKL_INT *col_indices = A.graph.entries.ptr_on_device() ;
    MKL_INT *row_beg = const_cast<MKL_INT*>(A.graph.row_map.ptr_on_device()) ;
    MKL_INT *row_end = row_beg+1;
    char matdescra[6] = { 'G', 'x', 'N', 'C', 'x', 'x' };
    char trans = 'N';
    double alpha = 1.0;
    double beta = 0.0;

    // Copy columns of x into a contiguous vector
    MKL_INT ncol = indices.size();
    trans_multi_vector_type xx( "xx" , ncol , n );
    trans_multi_vector_type yy( "yy" , ncol , n );
    Impl::GatherTranspose<value_type,ordinal_type,device_type>::apply(x,xx,indices);
    double *x_values = xx.ptr_on_device() ;
    double *y_values = yy.ptr_on_device() ;
    
    // Call MKLs CSR x multi-vector (row-based) multiply
    mkl_dcsrmm(&trans, &n, &ncol, &n, &alpha, matdescra, A_values, col_indices,
	       row_beg, row_end, x_values, &ncol, &beta, y_values, &ncol);

    // Copy columns out of continguous multivector
    Impl::ScatterTranspose<value_type,ordinal_type,device_type>::apply(y,yy,indices);

  }
};

template<>
class MMultiply<
  CrsMatrix< double , KokkosArray::Host > ,
  KokkosArray::View< double[] , KokkosArray::Host > ,
  KokkosArray::View< double[] , KokkosArray::Host > ,
  MKLSparseMatOps >
{
public:
  typedef KokkosArray::Host                                device_type ;
  typedef device_type::size_type                           size_type ;
  typedef KokkosArray::View< double[] , device_type > vector_type ;
  typedef CrsMatrix< double , device_type >           matrix_type ;
  typedef double                                      value_type ;
  typedef KokkosArray::View< double** , KokkosArray::LayoutLeft, device_type >  trans_multi_vector_type ;

public:

  static void apply( const matrix_type & A ,
                     const std::vector<vector_type> & x ,
                     const std::vector<vector_type> & y )
  {
    MKL_INT n = A.graph.row_map.dimension(0) - 1 ;
    double *A_values = A.values.ptr_on_device() ;
    MKL_INT *col_indices = A.graph.entries.ptr_on_device() ;
    MKL_INT *row_beg = const_cast<MKL_INT*>(A.graph.row_map.ptr_on_device()) ;
    MKL_INT *row_end = row_beg+1;
    char matdescra[6] = { 'G', 'x', 'N', 'C', 'x', 'x' };
    char trans = 'N';
    double alpha = 1.0;
    double beta = 0.0;

    // Copy columns of x into a contiguous vector
    MKL_INT ncol = x.size();
    trans_multi_vector_type xx( "xx" , ncol , n );
    trans_multi_vector_type yy( "yy" , ncol , n );
    Impl::GatherVecTranspose<value_type,device_type>::apply(x,xx);
    double *x_values = xx.ptr_on_device() ;
    double *y_values = yy.ptr_on_device() ;
    
    // Call MKLs CSR x multi-vector (row-based) multiply
    mkl_dcsrmm(&trans, &n, &ncol, &n, &alpha, matdescra, A_values, col_indices,
	       row_beg, row_end, x_values, &ncol, &beta, y_values, &ncol);

    // Copy columns out of continguous multivector
    Impl::ScatterVecTranspose<value_type,device_type>::apply(y,yy);

  }

};

template<>
class Multiply<
  CrsMatrix< float , KokkosArray::Host > ,
  KokkosArray::View< float[] , KokkosArray::Host > ,
  KokkosArray::View< float[] , KokkosArray::Host > ,
  MKLSparseMatOps >
{
public:
  typedef KokkosArray::Host                         device_type ;
  typedef device_type::size_type                    size_type ;
  typedef KokkosArray::View< float[] , device_type >  vector_type ;
  typedef CrsMatrix< float , device_type >    matrix_type ;

  static void apply( const matrix_type & A ,
                     const vector_type & x ,
                     const vector_type & y )
  {
    MKL_INT n = A.graph.row_map.dimension(0) - 1 ;
    float *A_values = A.values.ptr_on_device() ;
    MKL_INT *col_indices = A.graph.entries.ptr_on_device() ;
    MKL_INT *row_beg = const_cast<MKL_INT*>(A.graph.row_map.ptr_on_device()) ;
    MKL_INT *row_end = row_beg+1;
    char matdescra[6] = { 'G', 'x', 'N', 'C', 'x', 'x' };
    char trans = 'N';
    float alpha = 1.0;
    float beta = 0.0;

    float *x_values = x.ptr_on_device() ;
    float *y_values = y.ptr_on_device() ;
    
    mkl_scsrmv(&trans, &n, &n, &alpha, matdescra, A_values, col_indices,
	       row_beg, row_end, x_values, &beta, y_values);
  }
};

template<>
class MMultiply<
  CrsMatrix< float , KokkosArray::Host > ,
  KokkosArray::View< float** , KokkosArray::LayoutLeft, KokkosArray::Host > ,
  KokkosArray::View< float** , KokkosArray::LayoutLeft, KokkosArray::Host > ,
  MKLSparseMatOps >
{
public:
  typedef KokkosArray::Host device_type ;
  typedef device_type::size_type size_type ;
  typedef KokkosArray::View< float** , KokkosArray::LayoutLeft, device_type >  multi_vector_type ;
  typedef KokkosArray::View< float** , KokkosArray::LayoutLeft, device_type >  trans_multi_vector_type ;
  typedef KokkosArray::View< float[] , device_type > vector_type ;
  typedef CrsMatrix< float , device_type > matrix_type ;
  typedef float value_type ;
  typedef int ordinal_type ;

  static void apply( const matrix_type & A ,
                     const multi_vector_type & x ,
                     const multi_vector_type & y ,
		     const std::vector<ordinal_type> & indices)
  {
    MKL_INT n = A.graph.row_map.dimension(0) - 1 ;
    float *A_values = A.values.ptr_on_device() ;
    MKL_INT *col_indices = A.graph.entries.ptr_on_device() ;
    MKL_INT *row_beg = const_cast<MKL_INT*>(A.graph.row_map.ptr_on_device()) ;
    MKL_INT *row_end = row_beg+1;
    char matdescra[6] = { 'G', 'x', 'N', 'C', 'x', 'x' };
    char trans = 'N';
    float alpha = 1.0;
    float beta = 0.0;

    // Copy columns of x into a contiguous vector
    MKL_INT ncol = indices.size();
    trans_multi_vector_type xx( "xx" , ncol , n );
    trans_multi_vector_type yy( "yy" , ncol , n );
    Impl::GatherTranspose<value_type,ordinal_type,device_type>::apply(x,xx,indices);
    float *x_values = xx.ptr_on_device() ;
    float *y_values = yy.ptr_on_device() ;
    
    // Call MKLs CSR x multi-vector (row-based) multiply
    mkl_scsrmm(&trans, &n, &ncol, &n, &alpha, matdescra, A_values, col_indices,
	       row_beg, row_end, x_values, &ncol, &beta, y_values, &ncol);

    // Copy columns out of continguous multivector
    Impl::ScatterTranspose<value_type,ordinal_type,device_type>::apply(y,yy,indices);

  }
};

template<>
class MMultiply<
  CrsMatrix< float , KokkosArray::Host > ,
  KokkosArray::View< float[] , KokkosArray::Host > ,
  KokkosArray::View< float[] , KokkosArray::Host > ,
  MKLSparseMatOps >
{
public:
  typedef KokkosArray::Host                                device_type ;
  typedef device_type::size_type                           size_type ;
  typedef KokkosArray::View< float[] , device_type > vector_type ;
  typedef CrsMatrix< float , device_type >           matrix_type ;
  typedef float                                      value_type ;
  typedef KokkosArray::View< float** , KokkosArray::LayoutLeft, device_type >  trans_multi_vector_type ;

  static void apply( const matrix_type & A ,
                     const std::vector<vector_type> & x ,
                     const std::vector<vector_type> & y )
  {
    MKL_INT n = A.graph.row_map.dimension(0) - 1 ;
    float *A_values = A.values.ptr_on_device() ;
    MKL_INT *col_indices = A.graph.entries.ptr_on_device() ;
    MKL_INT *row_beg = const_cast<MKL_INT*>(A.graph.row_map.ptr_on_device()) ;
    MKL_INT *row_end = row_beg+1;
    char matdescra[6] = { 'G', 'x', 'N', 'C', 'x', 'x' };
    char trans = 'N';
    float alpha = 1.0;
    float beta = 0.0;

    // Copy columns of x into a contiguous vector
    MKL_INT ncol = x.size();
    trans_multi_vector_type xx( "xx" , ncol , n );
    trans_multi_vector_type yy( "yy" , ncol , n );
    Impl::GatherVecTranspose<value_type,device_type>::apply(x,xx);
    float *x_values = xx.ptr_on_device() ;
    float *y_values = yy.ptr_on_device() ;
    
    // Call MKLs CSR x multi-vector (row-based) multiply
    mkl_scsrmm(&trans, &n, &ncol, &n, &alpha, matdescra, A_values, col_indices,
	       row_beg, row_end, x_values, &ncol, &beta, y_values, &ncol);

    // Copy columns out of continguous multivector
    Impl::ScatterVecTranspose<value_type,device_type>::apply(y,yy);

  }

};

#endif

template< typename MatrixValue>
class MatrixMarketWriter<MatrixValue, KokkosArray::Host>
{
public:
  typedef KokkosArray::Host                         device_type ;
  typedef device_type::size_type                    size_type ;
  typedef CrsMatrix< MatrixValue , device_type >    matrix_type ;

  static void write(const matrix_type & A ,
		    const std::string& filename) {
    std::ofstream file(filename.c_str());
    file.precision(16);
    file.setf(std::ios::scientific);

    const size_type nRow = A.graph.row_count();
    
    // Write banner
    file << "%%MatrixMarket matrix coordinate real general" << std::endl;
    file << nRow << " " << nRow << " " << A.graph.entry_count() << std::endl;

    for (size_type row=0; row<nRow; ++row) {
      size_type entryBegin = A.graph.row_entry_begin(row);
      size_type entryEnd = A.graph.row_entry_end(row);
      for (size_type entry=entryBegin; entry<entryEnd; ++entry) {
	file << row+1 << " " << A.graph.column(entry)+1 << " " 
	     << std::setw(22) << A.values(entry) << std::endl;
      }
    }

    file.close();
  }
};

//----------------------------------------------------------------------------

} // namespace Stokhos

#endif /* #ifndef STOKHOS_HOST_CRSMATRIX_HPP */

