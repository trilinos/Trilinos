/*
// @HEADER
// ************************************************************************
//
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                    Copyright 2008 Sandia Corporation
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
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
// @HEADER

*/

#ifndef util_ArrayVector_hpp
#define util_ArrayVector_hpp

#include <vector>
#include <Array.hpp>

namespace phdmesh {

//----------------------------------------------------------------------
/** \brief Multidimensional array access wrapper to std::vector.
 */
template< class Dimension , typename DataType >
class Array< Dimension , std::vector<DataType> > {
private:
  typedef std::vector<DataType> VType ;
  Dimension             m_dim ;
  std::vector<DataType> m_vec ;
public:

  enum { DEEP_COPY = true };

  typedef Dimension dimension_type ;
  typedef DataType  data_type ;

  const dimension_type & dimension() const { return m_dim ; }

  data_type * data() const { return & const_cast<VType&>(m_vec)[0] ; }

  data_type & operator()( unsigned i1 , unsigned i2 ,
                          unsigned i3 , unsigned i4 ,
                          unsigned i5 , unsigned i6 ,
                          unsigned i7 , unsigned i8 ) const
    { return const_cast<VType&>(m_vec)[ m_dim(i1,i2,i3,i4,i5,i6,i7,i8) ]; }

  data_type & operator()( unsigned i1 , unsigned i2 ,
                          unsigned i3 , unsigned i4 ,
                          unsigned i5 , unsigned i6 ,
                          unsigned i7 ) const
    { return const_cast<VType&>(m_vec)[ m_dim(i1,i2,i3,i4,i5,i6,i7) ]; }

  data_type & operator()( unsigned i1 , unsigned i2 ,
                          unsigned i3 , unsigned i4 ,
                          unsigned i5 , unsigned i6 ) const
    { return const_cast<VType&>(m_vec)[ m_dim(i1,i2,i3,i4,i5,i6) ]; }

  data_type & operator()( unsigned i1 , unsigned i2 ,
                          unsigned i3 , unsigned i4 ,
                          unsigned i5 ) const
    { return const_cast<VType&>(m_vec)[ m_dim(i1,i2,i3,i4,i5) ]; }

  data_type & operator()( unsigned i1 , unsigned i2 ,
                          unsigned i3 , unsigned i4 ) const
    { return const_cast<VType&>(m_vec)[ m_dim(i1,i2,i3,i4) ]; }

  data_type & operator()( unsigned i1 , unsigned i2 ,
                          unsigned i3 ) const
    { return const_cast<VType&>(m_vec)[ m_dim(i1,i2,i3) ]; }

  data_type & operator()( unsigned i1 , unsigned i2 ) const
    { return const_cast<VType&>(m_vec)[ m_dim(i1,i2) ]; }

  data_type & operator()( unsigned i1 ) const
    { return const_cast<VType&>(m_vec)[ m_dim(i1) ]; }

  Array() : m_dim(), m_vec() {}

  Array( const Array & rhs )
    : m_dim( rhs.m_dim ), m_vec( rhs.m_vec ) {}

  Array & operator = ( const Array & rhs )
    { m_dim = rhs.m_dim ; m_vec = rhs.m_vec ; return *this ; }

  //--------------------------------

  template< class OtherDimension , class OtherStorage >
  Array( const Array< OtherDimension , OtherStorage > & rhs )
    : m_dim( rhs.dimension() ),
      m_vec( rhs.data() , rhs.data() + rhs.dimension().size() ) {}

  template< class OtherDimension , class OtherStorage >
  Array & operator = ( const Array< OtherDimension , OtherStorage > & rhs )
    {
      m_dim = rhs.dimension();
      m_vec.assign( rhs.data() , rhs.data() + rhs.dimension().size() );
      return *this ;
    }

  //--------------------------------
  // Class specific constructors:

  explicit Array( const dimension_type & arg_dim )
    : m_dim( arg_dim ), m_vec( arg_dim.size() ) {}

  explicit Array( unsigned n1 , unsigned n2 ,
                  unsigned n3 , unsigned n4 ,
                  unsigned n5 , unsigned n6 ,
                  unsigned n7 , unsigned n8 )
    : m_dim( n1 , n2 , n3 , n4 , n5 , n6 , n7 , n8 ),
      m_vec( n1 * n2 * n3 * n4 * n5 * n6 * n7 * n8 ) {}

  explicit Array( unsigned n1 , unsigned n2 ,
                  unsigned n3 , unsigned n4 ,
                  unsigned n5 , unsigned n6 ,
                  unsigned n7 )
    : m_dim( n1 , n2 , n3 , n4 , n5 , n6 , n7 ),
      m_vec( n1 * n2 * n3 * n4 * n5 * n6 * n7 ) {}

  explicit Array( unsigned n1 , unsigned n2 ,
                  unsigned n3 , unsigned n4 ,
                  unsigned n5 , unsigned n6 )
    : m_dim( n1 , n2 , n3 , n4 , n5 , n6 ),
      m_vec( n1 * n2 * n3 * n4 * n5 * n6 ) {}

  explicit Array( unsigned n1 , unsigned n2 ,
                  unsigned n3 , unsigned n4 ,
                  unsigned n5 )
    : m_dim( n1 , n2 , n3 , n4 , n5 ),
      m_vec( n1 * n2 * n3 * n4 * n5 ) {}

  explicit Array( unsigned n1 , unsigned n2 ,
                  unsigned n3 , unsigned n4 )
    : m_dim( n1 , n2 , n3 , n4 ),
      m_vec( n1 * n2 * n3 * n4 ) {}

  explicit Array( unsigned n1 , unsigned n2 , unsigned n3 )
    : m_dim( n1 , n2 , n3 ),
      m_vec( n1 * n2 * n3 ) {}

  explicit Array( unsigned n1 , unsigned n2 )
    : m_dim( n1 , n2 ),
      m_vec( n1 * n2 ) {}

  explicit Array( unsigned n1 )
    : m_dim( n1 ),
      m_vec( n1 ) {}

  //----------------------------------
  /** \brief Truncated view of the array */
  Array< typename Dimension::TruncatedType , DataType * >
    operator[]( unsigned i ) const
      {
        typedef typename Dimension::TruncatedType DimTrunc ;
        return Array< DimTrunc , DataType * >
                ( (const DimTrunc &) m_dim , data() + m_dim[i] );
      }
  //----------------------------------
};

//----------------------------------------------------------------------

}

#endif

