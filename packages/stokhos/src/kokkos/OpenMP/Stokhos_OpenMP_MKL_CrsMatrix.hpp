// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_OPENMP_MKL_CRSMATRIX_HPP
#define STOKHOS_OPENMP_MKL_CRSMATRIX_HPP

#include "Kokkos_Macros.hpp"
#include "Stokhos_ConfigDefs.h"
#if defined(KOKKOS_ENABLE_OPENMP) && defined(HAVE_STOKHOS_MKL)

#include "Kokkos_Core.hpp"
#include "Stokhos_CrsMatrix.hpp"
#include "mkl.h"

// Implementation of CrsMatrix multiply using MKL and OpenMP

namespace Stokhos {

namespace Impl {

  template <typename ValueType, typename OrdinalType, typename ExecutionSpace>
  struct GatherTranspose {
    typedef ValueType value_type;
    typedef OrdinalType ordinal_type;
    typedef ExecutionSpace execution_space;
    typedef typename execution_space::size_type size_type;
    typedef Kokkos::View< value_type** , Kokkos::LayoutLeft, execution_space >  multi_vector_type ;
    typedef Kokkos::View< value_type** , Kokkos::LayoutLeft, execution_space >  trans_multi_vector_type ;

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
      const size_t n = xt.extent(1);
      Kokkos::parallel_for( n, GatherTranspose(x,xt,indices) );
    }
  };

  template <typename ValueType, typename OrdinalType, typename ExecutionSpace>
  struct ScatterTranspose {
    typedef ValueType value_type;
    typedef OrdinalType ordinal_type;
    typedef ExecutionSpace execution_space;
    typedef typename execution_space::size_type size_type;
    typedef Kokkos::View< value_type** , Kokkos::LayoutLeft, execution_space >  multi_vector_type ;
    typedef Kokkos::View< value_type** , Kokkos::LayoutLeft, execution_space >  trans_multi_vector_type ;

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
      const size_t n = xt.extent(1);
      Kokkos::parallel_for( n, ScatterTranspose(x,xt,indices) );
    }
  };

  template <typename ValueType, typename ExecutionSpace>
  struct GatherVecTranspose {
    typedef ValueType value_type;
    typedef ExecutionSpace execution_space;
    typedef typename execution_space::size_type size_type;
    typedef Kokkos::View< value_type* , execution_space > vector_type ;
    typedef Kokkos::View< value_type** , Kokkos::LayoutLeft, execution_space >  trans_multi_vector_type ;

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
      const size_t n = xt.extent(1);
      Kokkos::parallel_for( n, GatherVecTranspose(x,xt) );
    }
  };

  template <typename ValueType, typename ExecutionSpace>
  struct ScatterVecTranspose {
    typedef ValueType value_type;
    typedef ExecutionSpace execution_space;
    typedef typename execution_space::size_type size_type;
    typedef Kokkos::View< value_type* , execution_space > vector_type ;
    typedef Kokkos::View< value_type** , Kokkos::LayoutLeft, execution_space >  trans_multi_vector_type ;

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
      const size_t n = xt.extent(1);
      Kokkos::parallel_for( n, ScatterVecTranspose(x,xt) );
    }
  };

} // namespace Impl

// Specializations of multiply using Intel MKL
class MKLMultiply {};

void multiply(const CrsMatrix< double , Kokkos::OpenMP >& A,
              const Kokkos::View< double* , Kokkos::OpenMP >& x,
              Kokkos::View< double* , Kokkos::OpenMP >& y,
              MKLMultiply tag)
{
  MKL_INT n = A.graph.row_map.extent(0) - 1 ;
  double *A_values = A.values.data() ;
  MKL_INT *col_indices = A.graph.entries.data() ;
  MKL_INT *row_beg = const_cast<MKL_INT*>(A.graph.row_map.data()) ;
  MKL_INT *row_end = row_beg+1;
  char matdescra[6] = { 'G', 'x', 'N', 'C', 'x', 'x' };
  char trans = 'N';
  double alpha = 1.0;
  double beta = 0.0;

  double *x_values = x.data() ;
  double *y_values = y.data() ;

  mkl_dcsrmv(&trans, &n, &n, &alpha, matdescra, A_values, col_indices,
             row_beg, row_end, x_values, &beta, y_values);
}

void multiply(const CrsMatrix< float , Kokkos::OpenMP >& A,
              const Kokkos::View< float* , Kokkos::OpenMP >& x,
              Kokkos::View< float* , Kokkos::OpenMP >& y,
              MKLMultiply tag)
{
  MKL_INT n = A.graph.row_map.extent(0) - 1 ;
  float *A_values = A.values.data() ;
  MKL_INT *col_indices = A.graph.entries.data() ;
  MKL_INT *row_beg = const_cast<MKL_INT*>(A.graph.row_map.data()) ;
  MKL_INT *row_end = row_beg+1;
  char matdescra[6] = { 'G', 'x', 'N', 'C', 'x', 'x' };
  char trans = 'N';
  float alpha = 1.0;
  float beta = 0.0;

  float *x_values = x.data() ;
  float *y_values = y.data() ;

  mkl_scsrmv(&trans, &n, &n, &alpha, matdescra, A_values, col_indices,
             row_beg, row_end, x_values, &beta, y_values);
}

void multiply(const CrsMatrix< double , Kokkos::OpenMP >& A,
              const std::vector< Kokkos::View< double* , Kokkos::OpenMP > >& x,
              std::vector< Kokkos::View< double* , Kokkos::OpenMP > >& y,
              MKLMultiply tag)
{
  typedef Kokkos::OpenMP execution_space ;
  typedef double value_type ;
  typedef Kokkos::View< double** , Kokkos::LayoutLeft, execution_space >  trans_multi_vector_type ;

  MKL_INT n = A.graph.row_map.extent(0) - 1 ;
  double *A_values = A.values.data() ;
  MKL_INT *col_indices = A.graph.entries.data() ;
  MKL_INT *row_beg = const_cast<MKL_INT*>(A.graph.row_map.data()) ;
  MKL_INT *row_end = row_beg+1;
  char matdescra[6] = { 'G', 'x', 'N', 'C', 'x', 'x' };
  char trans = 'N';
  double alpha = 1.0;
  double beta = 0.0;

  // Copy columns of x into a contiguous vector
  MKL_INT ncol = x.size();
  trans_multi_vector_type xx( "xx" , ncol , n );
  trans_multi_vector_type yy( "yy" , ncol , n );
  Impl::GatherVecTranspose<value_type,execution_space>::apply(x,xx);
  double *x_values = xx.data() ;
  double *y_values = yy.data() ;

  // Call MKLs CSR x multi-vector (row-based) multiply
  mkl_dcsrmm(&trans, &n, &ncol, &n, &alpha, matdescra, A_values, col_indices,
             row_beg, row_end, x_values, &ncol, &beta, y_values, &ncol);

  // Copy columns out of continguous multivector
  Impl::ScatterVecTranspose<value_type,execution_space>::apply(y,yy);
}

void multiply(const CrsMatrix< float , Kokkos::OpenMP >& A,
              const std::vector< Kokkos::View< float* , Kokkos::OpenMP > >& x,
              std::vector< Kokkos::View< float* , Kokkos::OpenMP > >& y,
              MKLMultiply tag)
{
  typedef Kokkos::OpenMP execution_space ;
  typedef float value_type ;
  typedef Kokkos::View< float** , Kokkos::LayoutLeft, execution_space >  trans_multi_vector_type ;

  MKL_INT n = A.graph.row_map.extent(0) - 1 ;
  float *A_values = A.values.data() ;
  MKL_INT *col_indices = A.graph.entries.data() ;
  MKL_INT *row_beg = const_cast<MKL_INT*>(A.graph.row_map.data()) ;
  MKL_INT *row_end = row_beg+1;
  char matdescra[6] = { 'G', 'x', 'N', 'C', 'x', 'x' };
  char trans = 'N';
  float alpha = 1.0;
  float beta = 0.0;

  // Copy columns of x into a contiguous vector
  MKL_INT ncol = x.size();
  trans_multi_vector_type xx( "xx" , ncol , n );
  trans_multi_vector_type yy( "yy" , ncol , n );
  Impl::GatherVecTranspose<value_type,execution_space>::apply(x,xx);
  float *x_values = xx.data() ;
  float *y_values = yy.data() ;

  // Call MKLs CSR x multi-vector (row-based) multiply
  mkl_scsrmm(&trans, &n, &ncol, &n, &alpha, matdescra, A_values, col_indices,
             row_beg, row_end, x_values, &ncol, &beta, y_values, &ncol);

  // Copy columns out of continguous multivector
  Impl::ScatterVecTranspose<value_type,execution_space>::apply(y,yy);
}

template <typename ordinal_type>
void multiply(
  const CrsMatrix< double , Kokkos::OpenMP >& A,
  const Kokkos::View< double** , Kokkos::LayoutLeft, Kokkos::OpenMP >& x,
  Kokkos::View< double** , Kokkos::LayoutLeft, Kokkos::OpenMP >& y,
  const std::vector<ordinal_type>& indices,
  MKLMultiply tag)
{
  typedef Kokkos::OpenMP execution_space ;
  typedef double value_type ;
  typedef Kokkos::View< double** , Kokkos::LayoutLeft, execution_space >  trans_multi_vector_type ;

  MKL_INT n = A.graph.row_map.extent(0) - 1 ;
  double *A_values = A.values.data() ;
  MKL_INT *col_indices = A.graph.entries.data() ;
  MKL_INT *row_beg = const_cast<MKL_INT*>(A.graph.row_map.data()) ;
  MKL_INT *row_end = row_beg+1;
  char matdescra[6] = { 'G', 'x', 'N', 'C', 'x', 'x' };
  char trans = 'N';
  double alpha = 1.0;
  double beta = 0.0;

  // Copy columns of x into a contiguous vector
  MKL_INT ncol = indices.size();
  trans_multi_vector_type xx( "xx" , ncol , n );
  trans_multi_vector_type yy( "yy" , ncol , n );
  Impl::GatherTranspose<value_type,ordinal_type,execution_space>::apply(x,xx,indices);
  double *x_values = xx.data() ;
  double *y_values = yy.data() ;

  // Call MKLs CSR x multi-vector (row-based) multiply
  mkl_dcsrmm(&trans, &n, &ncol, &n, &alpha, matdescra, A_values, col_indices,
             row_beg, row_end, x_values, &ncol, &beta, y_values, &ncol);

  // Copy columns out of continguous multivector
  Impl::ScatterTranspose<value_type,ordinal_type,execution_space>::apply(y,yy,indices);
}

template <typename ordinal_type>
void multiply(
  const CrsMatrix< float , Kokkos::OpenMP >& A,
  const Kokkos::View< float** , Kokkos::LayoutLeft, Kokkos::OpenMP >& x,
  Kokkos::View< float** , Kokkos::LayoutLeft, Kokkos::OpenMP >& y,
  const std::vector<ordinal_type>& indices,
  MKLMultiply tag)
{
  typedef Kokkos::OpenMP execution_space ;
  typedef float value_type ;
  typedef Kokkos::View< float** , Kokkos::LayoutLeft, execution_space >  trans_multi_vector_type ;

  MKL_INT n = A.graph.row_map.extent(0) - 1 ;
  float *A_values = A.values.data() ;
  MKL_INT *col_indices = A.graph.entries.data() ;
  MKL_INT *row_beg = const_cast<MKL_INT*>(A.graph.row_map.data()) ;
  MKL_INT *row_end = row_beg+1;
  char matdescra[6] = { 'G', 'x', 'N', 'C', 'x', 'x' };
  char trans = 'N';
  float alpha = 1.0;
  float beta = 0.0;

  // Copy columns of x into a contiguous vector
  MKL_INT ncol = indices.size();
  trans_multi_vector_type xx( "xx" , ncol , n );
  trans_multi_vector_type yy( "yy" , ncol , n );
  Impl::GatherTranspose<value_type,ordinal_type,execution_space>::apply(x,xx,indices);
  float *x_values = xx.data() ;
  float *y_values = yy.data() ;

  // Call MKLs CSR x multi-vector (row-based) multiply
  mkl_scsrmm(&trans, &n, &ncol, &n, &alpha, matdescra, A_values, col_indices,
             row_beg, row_end, x_values, &ncol, &beta, y_values, &ncol);

  // Copy columns out of continguous multivector
  Impl::ScatterTranspose<value_type,ordinal_type,execution_space>::apply(y,yy,indices);
}

//----------------------------------------------------------------------------

} // namespace Stokhos

#endif

#endif /* #ifndef STOKHOS_OPENMP_MKL_CRSMATRIX_HPP */
