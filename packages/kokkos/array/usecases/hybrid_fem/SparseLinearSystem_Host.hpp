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

#ifndef SPARSELINEARSYSTEM_HOST_HPP
#define SPARSELINEARSYSTEM_HOST_HPP

#include <SparseLinearSystem.hpp>
#include <KokkosArray_Host.hpp>

namespace KokkosArray {
namespace Impl {

template< typename AScalarType , typename VScalarType >
struct Multiply< CrsMatrix<AScalarType,Host> ,
                 View<VScalarType*,Host > ,
                 View<VScalarType*,Host > >
{
  typedef Host                                      device_type ;
  typedef device_type::size_type                    size_type ;
  typedef View< VScalarType* , device_type >  vector_type ;
  typedef CrsMatrix< AScalarType , device_type >    matrix_type ;

private:

  matrix_type  m_A ;
  vector_type  m_x ;
  vector_type  m_y ;

public:

  //--------------------------------------------------------------------------

  inline
  void operator()( const size_type iRow ) const
  {
    const size_type iEntryBegin = m_A.graph.row_map[iRow];
    const size_type iEntryEnd   = m_A.graph.row_map[iRow+1];

    double sum = 0 ;

    for ( size_type iEntry = iEntryBegin ; iEntry < iEntryEnd ; ++iEntry ) {
      sum += m_A.coefficients(iEntry) * m_x( m_A.graph.entries(iEntry) );
    }

    m_y(iRow) = sum ;
  }

  static void apply( const matrix_type & A ,
                     const size_type nrow ,
                     const size_type , // ncol ,
                     const vector_type & x ,
                     const vector_type & y )
  {
    Multiply op ;
    op.m_A = A ;
    op.m_x = x ;
    op.m_y = y ;
    parallel_for( nrow , op );
  }
};

} /* namespace Impl */
} /* namespace KokkosArray */


#endif /* #ifndef SPARSELINEARSYSTEM_HOST_HPP */

