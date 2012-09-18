/*
//@HEADER
// ************************************************************************
// 
//   KokkosArray: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
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
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOSARRAY_HOST_CRSMATRIX_HPP
#define KOKKOSARRAY_HOST_CRSMATRIX_HPP

#include <fstream>
#include <iomanip>

namespace KokkosArray {
namespace Impl {

template< typename MatrixValue , typename VectorValue >
class Multiply<
  CrsMatrix< MatrixValue , Host > ,
  KokkosArray::View< VectorValue[] , Host > ,
  KokkosArray::View< VectorValue[] , Host > >
{
public:
  typedef Host                                      device_type ;
  typedef device_type::size_type                    size_type ;
  typedef View< VectorValue[] , device_type >  vector_type ;
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
    parallel_for( row_count , Multiply(A,x,y) );
  }
};

template< typename MatrixValue , typename VectorValue >
class MMultiply<
  CrsMatrix< MatrixValue , Host > ,
  KokkosArray::View< VectorValue[] , Host > ,
  KokkosArray::View< VectorValue[] , Host > >
{
public:
  typedef Host                                      device_type ;
  typedef device_type::size_type                    size_type ;
  typedef View< VectorValue[] , device_type >       vector_type ;
  typedef CrsMatrix< MatrixValue , device_type >    matrix_type ;
  typedef VectorValue                               value_type ;

  const matrix_type  m_A ;
  const std::vector<vector_type>  m_x ;
  const std::vector<vector_type>  m_y ;

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
    const size_t num_vecs = m_x.size();
    
    std::vector<double> sum(num_vecs, 0);

    for ( size_type iEntry = iEntryBegin ; iEntry < iEntryEnd ; ++iEntry ) {
      double a = m_A.values(iEntry);
      for (size_t col=0; col<num_vecs; col++) 
	sum[col] += a * m_x[col]( m_A.graph.entries(iEntry) );
    }
    
    for (size_t col=0; col<num_vecs; col++)
      m_y[col](iRow) = sum[col] ;
  }

  static void apply( const matrix_type & A ,
                     const std::vector<vector_type> & x ,
                     const std::vector<vector_type> & y )
  {
    parallel_for( A.graph.row_map.length() , MMultiply(A,x,y) );
  }
};

template< typename MatrixValue>
class MatrixMarketWriter<MatrixValue,Host>
{
public:
  typedef Host                                      device_type ;
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

} // namespace Impl
} // namespace KokkosArray

#endif /* #ifndef KOKKOSARRAY_HOST_CRSMATRIX_HPP */

