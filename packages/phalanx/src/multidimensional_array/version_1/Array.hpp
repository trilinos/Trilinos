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

#ifndef util_Array_hpp
#define util_Array_hpp

#include <Dimension.hpp>

namespace phdmesh {

//----------------------------------------------------------------------
/** \class Array
 *  \brief A multidimensional array is the composition of a
 *         contigous block of member data and a compatible
 *         multi-index mapping.
 *
 *  \param Dimension
 *  Template parameter for a multi-index mapping class.
 *  This class is required to provide:
 *  -  enum { Rank = ... };
 *     which is the rank of the dimension and
 *  -  unsigned operator()( unsigned i1 , ... ) const ;
 *     which is the multi-index to offset map.
 *
 * \param StorageType
 *  Template parameter for the type of contiguous block of member data.
 *  It is expected for the Array template class to be specialized
 *  for each StorageType.
 *  
 *  An Array specialization must provide the following
 *  type definitions and methods:
 *  -#  typedef ... data_type ;
 *  -#  typedef ... dimension_type ;
 *  -#  const dimension_type & dimension() const ;
 *  -#  data_type * data() const ;
 *  -#  data_type & operator()( unsigned i1 , ... ) const ;
 *
 *  Constructors and assignment operators may be provided;
 *  however, these are specific to the storage type.
 */
template< class Dimension , typename StorageType > class Array ;

//----------------------------------------------------------------------
/** \brief  Multidimensional array access to an anonymous block of members.
 *
 *  This specialization for the multidimensional rray with a
 *  simple pointer storage type has a constructor that accepts
 *  a dimension and a data_type pointer.
 *  This constructor trusts the user to input a pointer that references
 *  a contiguous block of memory which is equal to or greather than
 *  the size of the dimension.  This specialization also has
 *  a copy constructor and assignment operator that accepts any
 *  other multidimensional array type with compatible
 *  'dimension_type' and 'data_type'.
 */
template< class Dimension , typename DataType >
class Array< Dimension , DataType * > {
private:
  Dimension  m_dim ;
  DataType * m_ptr ;
public:

  //----------------------------------
  /** \brief Does this array type perform a deep copy */
  enum { DEEP_COPY = false };

  typedef DataType  data_type ;
  typedef Dimension dimension_type ;

  //----------------------------------
  /** \brief Query the dimension */
  const dimension_type & dimension() const { return m_dim ; }

  /** \brief Pointer to contiguous block of memory */
  data_type * data() const { return m_ptr ; }

  //----------------------------------
  /** \brief Access member via Rank 8 multi-index */
  data_type & operator()( unsigned i1 , unsigned i2 ,
                          unsigned i3 , unsigned i4 ,
                          unsigned i5 , unsigned i6 ,
                          unsigned i7 , unsigned i8 ) const
    {
      DimHelper< dimension_type::Rank , 8 >::good();
      return m_ptr[ m_dim(i1,i2,i3,i4,i5,i6,i7,i8) ];
    }

  /** \brief Access member via Rank 7 multi-index */
  data_type & operator()( unsigned i1 , unsigned i2 ,
                          unsigned i3 , unsigned i4 ,
                          unsigned i5 , unsigned i6 ,
                          unsigned i7 ) const
    { return m_ptr[ m_dim(i1,i2,i3,i4,i5,i6,i7) ]; }

  /** \brief Access member via Rank 6 multi-index */
  data_type & operator()( unsigned i1 , unsigned i2 ,
                          unsigned i3 , unsigned i4 ,
                          unsigned i5 , unsigned i6 ) const
    { return m_ptr[ m_dim(i1,i2,i3,i4,i5,i6) ]; }

  /** \brief Access member via Rank 5 multi-index */
  data_type & operator()( unsigned i1 , unsigned i2 ,
                          unsigned i3 , unsigned i4 ,
                          unsigned i5 ) const
    { return m_ptr[ m_dim(i1,i2,i3,i4,i5) ]; }

  /** \brief Access member via Rank 4 multi-index */
  data_type & operator()( unsigned i1 , unsigned i2 ,
                          unsigned i3 , unsigned i4 ) const
    { return m_ptr[ m_dim(i1,i2,i3,i4) ]; }

  /** \brief Access member via Rank 3 multi-index */
  data_type & operator()( unsigned i1 , unsigned i2 ,
                          unsigned i3 ) const
    { return m_ptr[ m_dim(i1,i2,i3) ]; }

  /** \brief Access member via Rank 2 multi-index */
  data_type & operator()( unsigned i1 , unsigned i2 ) const
    { return m_ptr[ m_dim(i1,i2) ]; }

  /** \brief Access member via Rank 1 multi-index */
  data_type & operator()( unsigned i1 ) const
    { return m_ptr[ m_dim(i1) ]; }

  //----------------------------------
  // Required constructors and assignment operators:

  Array() : m_dim(), m_ptr(NULL) {}

  Array( const Array & rhs ) : m_dim( rhs.m_dim ), m_ptr( rhs.m_ptr ) {}

  Array & operator = ( const Array & rhs )
    { m_dim = rhs.m_dim ; m_ptr = rhs.m_ptr ; return *this ; }

  //----------------------------------
  // Take pointer to other compatible array types.

  template< class OtherDimension , class OtherStorage >
  Array( const Array< OtherDimension , OtherStorage > & rhs )
    : m_dim( rhs.dimension() ), // Dimension must be compatible
      m_ptr( rhs.data() ) {}    // Data type must be compatible

  template< class OtherDimension , class OtherStorage >
  Array & operator = ( const Array< OtherDimension , OtherStorage > & rhs )
    {
      m_dim = rhs.dimension(); // Dimension must be compatible
      m_ptr = rhs.data();      // Data type must be compatible
      return *this ;
    }

  //----------------------------------
  // Class specific constructors:

  Array( const dimension_type & arg_dim , data_type * arg_ptr )
    : m_dim( arg_dim ), m_ptr( arg_ptr ) {}

  Array( unsigned n1 , unsigned n2 , unsigned n3 , unsigned n4 ,
         unsigned n5 , unsigned n6 , unsigned n7 , unsigned n8 ,
         data_type * arg_ptr )
    : m_dim( n1 , n2 , n3 , n4 , n5 , n6 , n7 , n8 ), m_ptr( arg_ptr ) {}

  Array( unsigned n1 , unsigned n2 , unsigned n3 , unsigned n4 ,
         unsigned n5 , unsigned n6 , unsigned n7 ,
         data_type * arg_ptr )
    : m_dim( n1 , n2 , n3 , n4 , n5 , n6 , n7 ), m_ptr( arg_ptr ) {}

  Array( unsigned n1 , unsigned n2 , unsigned n3 , unsigned n4 ,
         unsigned n5 , unsigned n6 ,
         data_type * arg_ptr )
    : m_dim( n1 , n2 , n3 , n4 , n5 , n6 ), m_ptr( arg_ptr ) {}

  Array( unsigned n1 , unsigned n2 , unsigned n3 , unsigned n4 ,
         unsigned n5 ,
         data_type * arg_ptr )
    : m_dim( n1 , n2 , n3 , n4 , n5 ), m_ptr( arg_ptr ) {}

  Array( unsigned n1 , unsigned n2 , unsigned n3 , unsigned n4 ,
         data_type * arg_ptr )
    : m_dim( n1 , n2 , n3 , n4 ), m_ptr( arg_ptr ) {}

  Array( unsigned n1 , unsigned n2 , unsigned n3 ,
         data_type * arg_ptr )
    : m_dim( n1 , n2 , n3 ), m_ptr( arg_ptr ) {}

  Array( unsigned n1 , unsigned n2 ,
         data_type * arg_ptr )
    : m_dim( n1 , n2 ), m_ptr( arg_ptr ) {}

  Array( unsigned n1 ,
         data_type * arg_ptr )
    : m_dim( n1 ), m_ptr( arg_ptr ) {}

  //----------------------------------
  /** \brief Truncated view of the array */
  Array< typename Dimension::TruncatedType , DataType * >
    operator[]( unsigned i ) const
      {
        typedef typename Dimension::TruncatedType DimTrunc ;
        return Array< DimTrunc , DataType * >
                ( (const DimTrunc &) m_dim , m_ptr + m_dim[i] );
      }
  //----------------------------------
};

//----------------------------------------------------------------------

}

#endif

