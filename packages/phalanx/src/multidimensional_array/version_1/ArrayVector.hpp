/*------------------------------------------------------------------------*/
/*      phdMesh : Parallel Heterogneous Dynamic unstructured Mesh         */
/*                Copyright (2008) Sandia Corporation                     */
/*                                                                        */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*                                                                        */
/*  This library is free software; you can redistribute it and/or modify  */
/*  it under the terms of the GNU Lesser General Public License as        */
/*  published by the Free Software Foundation; either version 2.1 of the  */
/*  License, or (at your option) any later version.                       */
/*                                                                        */
/*  This library is distributed in the hope that it will be useful,       */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU     */
/*  Lesser General Public License for more details.                       */
/*                                                                        */
/*  You should have received a copy of the GNU Lesser General Public      */
/*  License along with this library; if not, write to the Free Software   */
/*  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307   */
/*  USA                                                                   */
/*------------------------------------------------------------------------*/
/**
 * @author H. Carter Edwards  <hcedwar@sandia.gov>
 * @date   June 2008
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

