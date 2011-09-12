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

//----------------------------------------------------------------------

namespace phdmesh {

/** \class ArrayNatural
 *  \brief A multidimensional array index mapping into a contiguous storage.
 *
 *  The Scalar template parameter defines the member type of the array.
 *  The Tag? template parameters define the ordinates.
 *
 *  Template parameter for a multi-index mapping class.
 *  This class is required to provide:
 *  -  enum { Rank = ... };
 *     which is the rank of the dimension and
 *  -  unsigned operator()( unsigned i1 , ... ) const ;
 *     which is the multi-index to offset map.
 *
 */
template< typename Scalar ,
          class Tag1 = void , class Tag2 = void ,
          class Tag3 = void , class Tag4 = void ,
          class Tag5 = void , class Tag6 = void ,
          class Tag7 = void , class Tag8 = void >
class ArrayNatural ;

template< typename Scalar ,
          class Tag1 = void , class Tag2 = void ,
          class Tag3 = void , class Tag4 = void ,
          class Tag5 = void , class Tag6 = void ,
          class Tag7 = void , class Tag8 = void >
class ArrayFortran ;

template< class Tag1 = void , class Tag2 = void ,
          class Tag3 = void , class Tag4 = void ,
          class Tag5 = void , class Tag6 = void ,
          class Tag7 = void , class Tag8 = void >
class ArrayDimension ;

class ArrayDimTag ;

}

//----------------------------------------------------------------------

#include <ArrayHelper.hpp>

//----------------------------------------------------------------------

namespace phdmesh {

struct ArrayDimTag {
  virtual const char * name() const = 0 ;

  virtual std::string to_string( size_t , int ) const ;
  virtual int          to_index( size_t , const std::string & ) const ;

  // A derived type must provide the following static method:
  //
  // static const Derived_ArrayDimTag_Type & descriptor();

protected:
  virtual ~ArrayDimTag();
  ArrayDimTag() {}
private:
  ArrayDimTag( const ArrayDimTag & );
  ArrayDimTag & operator = ( const ArrayDimTag & );
};

//----------------------------------------------------------------------

/** \cond */

template<>
class ArrayDimension<void,void,void,void,void,void,void,void> {
public:
  size_t              rank ;
  size_t              dimension[8];
  const ArrayDimTag * tags[8];
};

template< class Tag1 >
class ArrayDimension<Tag1,void,void,void,void,void,void,void> {
public:
  enum { Rank = 1 };
  size_t dimension[ Rank ];
};

template< class Tag1 , class Tag2 >
class ArrayDimension<Tag1,Tag2,void,void,void,void,void,void> {
public:
  enum { Rank = 2 };
  size_t dimension[ Rank ];
};

template< class Tag1 , class Tag2 , class Tag3 >
class ArrayDimension<Tag1,Tag2,Tag3,void,void,void,void,void> {
public:
  enum { Rank = 3 };
  size_t dimension[ Rank ];
};

template< class Tag1 , class Tag2 , class Tag3 , class Tag4 >
class ArrayDimension<Tag1,Tag2,Tag3,Tag4,void,void,void,void> {
public:
  enum { Rank = 4 };
  size_t dimension[ Rank ];
};

template< class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 >
class ArrayDimension<Tag1,Tag2,Tag3,Tag4,Tag5,void,void,void> {
public:
  enum { Rank = 5 };
  size_t dimension[ Rank ];
};

template< class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 >
class ArrayDimension<Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,void,void> {
public:
  enum { Rank = 6 };
  size_t dimension[ Rank ];
};

template< class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 >
class ArrayDimension<Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,void> {
public:
  enum { Rank = 7 };
  size_t dimension[ Rank ];
};

template< class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
class ArrayDimension {
public:
  enum { Rank = 8 };
  size_t dimension[ Rank ];
};

//----------------------------------------------------------------------

template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
class ArrayNatural
{
public:
  // Required by all arrays:
  typedef Scalar  value_type ;
  typedef int     index_type ;
  typedef size_t  size_type ;
  typedef const ArrayDimTag * tag_type ;

  // Required by compile-time knowledge arrays:

  typedef ArrayDimTagList<Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8> TagList ;

  enum { Rank       = TagList::Rank };
  enum { Natural    = true };
  enum { Reverse    = false };
  enum { Contiguous = true };

  template < unsigned Ord >
  struct Tag {
    typedef typename ArrayDimTagListAt<TagList,Ord>::type type ;
  };

private:

  typedef ArrayHelper< TagList::Rank > Helper ;

  Scalar * m_ptr ;
  size_t   m_stride[ TagList::Rank ];

  template< typename , class , class , class , class ,
                       class , class , class , class >
  friend class ArrayNatural ;

  template< typename , class , class , class , class ,
                       class , class , class , class >
  friend class ArrayFortran ;

public:

  //----------------------------------

  unsigned rank() const { return Rank ; }
  bool natural() const { return Natural ; }
  bool reverse() const { return Reverse ; }
  bool contiguous() const { return Contiguous ; }

  tag_type tag( unsigned ord ) const
    {
      array_bounds_checking( Rank , ord );
      return array_tags<TagList>()[ord] ;
    }

  //----------------------------------

  size_type size() const { return m_stride[ Rank - 1 ]; }

  template < unsigned Ord > size_type dimension() const
    {
      array_bounds_check_ordinal_is_less<Rank,Ord>();
      enum { I = ( Rank - 1 ) - Ord };
      return I ? m_stride[I] / m_stride[I-1] : m_stride[I] ;
    }

  size_type dimension( index_type ord ) const
    {
      array_bounds_checking( Rank , ord );
      const int i = ( Rank - 1 ) - ord ;
      return i ? m_stride[i] / m_stride[i-1] : m_stride[i] ;
    }

  void dimensions( std::vector<size_type> & n )
    {
      n.resize( Rank );
      n[Rank-1] = m_stride[0] ;
      for ( int i = 1 ; i < Rank ; ++i ) {
        n[ ( Rank - 1 ) - i ] = m_stride[i] / m_stride[i-1] ;
      }
    }
  
  //----------------------------------
  /** \brief Access member data */
  value_type * contiguous_data() const { return m_ptr ; }

  /** \brief Access member via full ordering of members. */
  value_type & operator[]( size_type i ) const
    {
      ARRAY_BOUNDS_CHECKING_ORDINAL( size() , i );
      return m_ptr[ i ];
    }

  /** \brief Access member via Rank 8 multi-index */
  value_type & operator()( index_type i1 , index_type i2 ,
                           index_type i3 , index_type i4 ,
                           index_type i5 , index_type i6 ,
                           index_type i7 , index_type i8 ) const
    {
      array_bounds_check_rank_is_equal<Rank,8>();
      ARRAY_BOUNDS_CHECKING_8(i1,i2,i3,i4,i5,i6,i7,i8);
      return m_ptr[ ARRAY_INDEX_NATURAL_8(m_stride,i1,i2,i3,i4,i5,i6,i7,i8) ];
    }

  /** \brief Access member via Rank 7 multi-index */
  value_type & operator()( index_type i1 , index_type i2 ,
                           index_type i3 , index_type i4 ,
                           index_type i5 , index_type i6 ,
                           index_type i7 ) const
    {
      array_bounds_check_rank_is_equal<Rank,7>();
      ARRAY_BOUNDS_CHECKING_7(i1,i2,i3,i4,i5,i6,i7);
      return m_ptr[ ARRAY_INDEX_NATURAL_7(m_stride,i1,i2,i3,i4,i5,i6,i7) ];
    }

  /** \brief Access member via Rank 6 multi-index */
  value_type & operator()( index_type i1 , index_type i2 ,
                           index_type i3 , index_type i4 ,
                           index_type i5 , index_type i6 ) const
    {
      array_bounds_check_rank_is_equal<Rank,6>();
      ARRAY_BOUNDS_CHECKING_6(i1,i2,i3,i4,i5,i6);
      return m_ptr[ ARRAY_INDEX_NATURAL_6(m_stride,i1,i2,i3,i4,i5,i6) ];
    }

  /** \brief Access member via Rank 5 multi-index */
  value_type & operator()( index_type i1 , index_type i2 ,
                           index_type i3 , index_type i4 ,
                           index_type i5 ) const
    {
      array_bounds_check_rank_is_equal<Rank,5>();
      ARRAY_BOUNDS_CHECKING_5(i1,i2,i3,i4,i5);
      return m_ptr[ ARRAY_INDEX_NATURAL_5(m_stride,i1,i2,i3,i4,i5) ];
    }

  /** \brief Access member via Rank 4 multi-index */
  value_type & operator()( index_type i1 , index_type i2 ,
                           index_type i3 , index_type i4 ) const
    {
      array_bounds_check_rank_is_equal<Rank,4>();
      ARRAY_BOUNDS_CHECKING_4(i1,i2,i3,i4);
      return m_ptr[ ARRAY_INDEX_NATURAL_4(m_stride,i1,i2,i3,i4) ];
    }

  /** \brief Access member via Rank 3 multi-index */
  value_type & operator()( index_type i1 , index_type i2 ,
                           index_type i3 ) const
    {
      array_bounds_check_rank_is_equal<Rank,3>();
      ARRAY_BOUNDS_CHECKING_3(i1,i2,i3);
      return m_ptr[ ARRAY_INDEX_NATURAL_3(m_stride,i1,i2,i3) ];
    }

  /** \brief Access member via Rank 2 multi-index */
  value_type & operator()( index_type i1 , index_type i2 ) const
    {
      array_bounds_check_rank_is_equal<Rank,2>();
      ARRAY_BOUNDS_CHECKING_2(i1,i2);
      return m_ptr[ ARRAY_INDEX_NATURAL_2(m_stride,i1,i2) ];
    }

  /** \brief Access member via Rank 1 multi-index */
  value_type & operator()( index_type i1 ) const
    {
      array_bounds_check_rank_is_equal<Rank,1>();
      ARRAY_BOUNDS_CHECKING_1(i1);
      return m_ptr[ i1 ];
    }

  //----------------------------------
  // Required constructors and assignment operators:

  ArrayNatural() : m_ptr(NULL) { Helper::zero( m_stride ); }

  ArrayNatural( const ArrayNatural & rhs )
    : m_ptr( rhs.m_ptr ) { Helper::copy( m_stride , rhs.m_stride ); }

  ArrayNatural & operator = ( const ArrayNatural & rhs )
    {
      m_ptr = rhs.m_ptr ;
      Helper::copy( m_stride , rhs.m_stride );
      return *this ;
    }

  //----------------------------------
  // Copy a reverse-map array

  typedef typename
    ArrayReverse<Scalar,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>::fortran_type
      ReverseType ;

  ArrayNatural( const ReverseType & rhs )
    : m_ptr( rhs.m_ptr ) { Helper::copy( m_stride , rhs.m_stride ); }

  ArrayNatural & operator = ( const ReverseType & rhs )
    {
      m_ptr = rhs.m_ptr ;
      Helper::copy( m_stride , rhs.m_stride );
      return *this ;
    }

  //----------------------------------
  // Truncated view

  typedef typename
    ArrayTruncate<Scalar,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>::natural_type
      TruncatedViewType ;

  TruncatedViewType truncated_view( index_type i ) const
    {
      return TruncatedViewType( m_ptr + m_stride[ Rank - 2 ] * i ,
                                ArrayStride(), m_stride );
    }

  // "Power user" constructor, promising a properly strided input
  ArrayNatural( value_type * arg_ptr ,
                const ArrayStride & , const size_type * arg_stride )
    : m_ptr( arg_ptr ) { Helper::copy( m_stride , arg_stride ); }

  // "Power user" constructor, promising a properly strided input
  ArrayNatural( value_type * arg_ptr ,
                const ArrayStride & ,
                const size_type * arg_stride ,
                const size_type   arg_append )
    : m_ptr( arg_ptr )
    { Helper::stride_append( m_stride , arg_stride , arg_append ); }

  //----------------------------------
  // Take pointer to other compatible array types.

  template< class OtherArrayType >
  explicit ArrayNatural( const OtherArrayType & rhs )
    : m_ptr( rhs.contiguous_data() )
    { Helper::stride( m_stride , rhs ); }

  template< class OtherArrayType >
  ArrayNatural & operator = ( const OtherArrayType & rhs )
    {
      m_ptr = rhs.contiguous_data();
      Helper::stride( m_stride , rhs );
      return *this ;
    }

  //----------------------------------
  // Class specific constructors:

  ArrayNatural( value_type * arg_ptr ,
         size_type n1 , size_type n2 , size_type n3 , size_type n4 ,
         size_type n5 , size_type n6 , size_type n7 , size_type n8 )
    : m_ptr( arg_ptr )
    {
      array_bounds_check_rank_is_equal<Rank,8>();
      ARRAY_STRIDE_NATURAL_8( m_stride, n1, n2, n3, n4, n5, n6, n7, n8 );
    }

  ArrayNatural( value_type * arg_ptr ,
         size_type n1 , size_type n2 , size_type n3 , size_type n4 ,
         size_type n5 , size_type n6 , size_type n7 )
    : m_ptr( arg_ptr )
    {
      array_bounds_check_rank_is_equal<Rank,7>();
      ARRAY_STRIDE_NATURAL_7( m_stride, n1, n2, n3, n4, n5, n6, n7 );
    }

  ArrayNatural( value_type * arg_ptr ,
         size_type n1 , size_type n2 , size_type n3 , size_type n4 ,
         size_type n5 , size_type n6 )
    : m_ptr( arg_ptr )
    {
      array_bounds_check_rank_is_equal<Rank,6>();
      ARRAY_STRIDE_NATURAL_6( m_stride, n1, n2, n3, n4, n5, n6 );
    }

  ArrayNatural( value_type * arg_ptr ,
         size_type n1 , size_type n2 , size_type n3 , size_type n4 ,
         size_type n5 )
    : m_ptr( arg_ptr )
    {
      array_bounds_check_rank_is_equal<Rank,5>();
      ARRAY_STRIDE_NATURAL_5( m_stride, n1, n2, n3, n4, n5 );
    }

  ArrayNatural( value_type * arg_ptr ,
         size_type n1 , size_type n2 , size_type n3 , size_type n4 )
    : m_ptr( arg_ptr )
    {
      array_bounds_check_rank_is_equal<Rank,4>();
      ARRAY_STRIDE_NATURAL_4( m_stride, n1, n2, n3, n4 );
    }

  ArrayNatural( value_type * arg_ptr ,
         size_type n1 , size_type n2 , size_type n3 )
    : m_ptr( arg_ptr )
    {
      array_bounds_check_rank_is_equal<Rank,3>();
      ARRAY_STRIDE_NATURAL_3( m_stride, n1, n2, n3 );
    }

  ArrayNatural( value_type * arg_ptr ,
         size_type n1 , size_type n2 )
    : m_ptr( arg_ptr )
    {
      array_bounds_check_rank_is_equal<Rank,2>();
      ARRAY_STRIDE_NATURAL_2( m_stride, n1, n2 );
    }

  ArrayNatural( value_type * arg_ptr , size_type n1 )
    : m_ptr( arg_ptr )
    {
      array_bounds_check_rank_is_equal<Rank,1>();
      m_stride[0] = n1 ;
    }
};

//----------------------------------------------------------------------

template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
class ArrayFortran
{
public:
  // Required by all arrays:
  typedef Scalar           value_type ;
  typedef int              index_type ;
  typedef size_t           size_type ;
  typedef const ArrayDimTag * tag_type ;

  // Required by compile-time knowledge arrays:

  typedef ArrayDimTagList<Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8> TagList ;

  enum { Rank       = TagList::Rank };
  enum { Natural    = false };
  enum { Reverse    = true };
  enum { Contiguous = true };

  template < unsigned Ord >
  struct Tag {
    typedef typename ArrayDimTagListAt<TagList,Ord>::type type ;
  };

private:

  typedef ArrayHelper< TagList::Rank > Helper ;

  Scalar * m_ptr ;
  size_t m_stride[ TagList::Rank ];

  template< typename ,
            class , class , class , class ,
            class , class , class , class >
  friend class ArrayNatural ;

  template< typename ,
            class , class , class , class ,
            class , class , class , class >
  friend class ArrayFortran ;

public:

  //----------------------------------

  unsigned rank() const { return Rank ; }
  bool natural() const { return Natural ; }
  bool reverse() const { return Reverse ; }
  bool contiguous() const { return Contiguous ; }

  tag_type tag( unsigned ord ) const
    {
      array_bounds_checking( Rank , ord );
      return array_tags<TagList>()[ord] ;
    }

  //----------------------------------

  size_type size() const { return m_stride[ Rank - 1 ]; }

  template < unsigned Ord > size_type dimension() const
    {
      array_bounds_check_ordinal_is_less<Rank,Ord>();
      return Ord ? m_stride[Ord] / m_stride[Ord-1] : m_stride[Ord] ;
    }

  size_type dimension( index_type ord ) const
    {
      array_bounds_checking( Rank , ord );
      return ord ? m_stride[ord] / m_stride[ord-1] : m_stride[ord] ;
    }

  void dimensions( std::vector<size_type> & n )
    {
      n.resize( Rank );
      n[0] = m_stride[0] ;
      for ( int i = 1 ; i < Rank ; ++i ) {
        n[i] = m_stride[i] / m_stride[i-1] ;
      }
    }
  
  //----------------------------------
  /** \brief Access member data */
  value_type * contiguous_data() const { return m_ptr ; }

  /** \brief Access member via full ordering of members. */
  value_type & operator[]( size_type i ) const
    {
      ARRAY_BOUNDS_CHECKING_ORDINAL( size() , i );
      return m_ptr[ i ];
    }

  /** \brief Access member via Rank 8 multi-index */
  value_type & operator()( index_type i1 , index_type i2 ,
                           index_type i3 , index_type i4 ,
                           index_type i5 , index_type i6 ,
                           index_type i7 , index_type i8 ) const
    {
      array_bounds_check_rank_is_equal<Rank,8>();
      ARRAY_BOUNDS_CHECKING_8(i1,i2,i3,i4,i5,i6,i7,i8);
      return m_ptr[ ARRAY_INDEX_FORTRAN_8(m_stride,i1,i2,i3,i4,i5,i6,i7,i8) ];
    }

  /** \brief Access member via Rank 7 multi-index */
  value_type & operator()( index_type i1 , index_type i2 ,
                           index_type i3 , index_type i4 ,
                           index_type i5 , index_type i6 ,
                           index_type i7 ) const
    {
      array_bounds_check_rank_is_equal<Rank,7>();
      ARRAY_BOUNDS_CHECKING_7(i1,i2,i3,i4,i5,i6,i7);
      return m_ptr[ ARRAY_INDEX_FORTRAN_7(m_stride,i1,i2,i3,i4,i5,i6,i7) ];
    }

  /** \brief Access member via Rank 6 multi-index */
  value_type & operator()( index_type i1 , index_type i2 ,
                           index_type i3 , index_type i4 ,
                           index_type i5 , index_type i6 ) const
    {
      array_bounds_check_rank_is_equal<Rank,6>();
      ARRAY_BOUNDS_CHECKING_6(i1,i2,i3,i4,i5,i6);
      return m_ptr[ ARRAY_INDEX_FORTRAN_6(m_stride,i1,i2,i3,i4,i5,i6) ];
    }

  /** \brief Access member via Rank 5 multi-index */
  value_type & operator()( index_type i1 , index_type i2 ,
                           index_type i3 , index_type i4 ,
                           index_type i5 ) const
    {
      array_bounds_check_rank_is_equal<Rank,5>();
      ARRAY_BOUNDS_CHECKING_5(i1,i2,i3,i4,i5);
      return m_ptr[ ARRAY_INDEX_FORTRAN_5(m_stride,i1,i2,i3,i4,i5) ];
    }

  /** \brief Access member via Rank 4 multi-index */
  value_type & operator()( index_type i1 , index_type i2 ,
                           index_type i3 , index_type i4 ) const
    {
      array_bounds_check_rank_is_equal<Rank,4>();
      ARRAY_BOUNDS_CHECKING_4(i1,i2,i3,i4);
      return m_ptr[ ARRAY_INDEX_FORTRAN_4(m_stride,i1,i2,i3,i4) ];
    }

  /** \brief Access member via Rank 3 multi-index */
  value_type & operator()( index_type i1 , index_type i2 ,
                           index_type i3 ) const
    {
      array_bounds_check_rank_is_equal<Rank,3>();
      ARRAY_BOUNDS_CHECKING_3(i1,i2,i3);
      return m_ptr[ ARRAY_INDEX_FORTRAN_3(m_stride,i1,i2,i3) ];
    }

  /** \brief Access member via Rank 2 multi-index */
  value_type & operator()( index_type i1 , index_type i2 ) const
    {
      array_bounds_check_rank_is_equal<Rank,2>();
      ARRAY_BOUNDS_CHECKING_2(i1,i2);
      return m_ptr[ ARRAY_INDEX_FORTRAN_2(m_stride,i1,i2) ];
    }

  /** \brief Access member via Rank 1 multi-index */
  value_type & operator()( index_type i1 ) const
    {
      array_bounds_check_rank_is_equal<Rank,1>();
      ARRAY_BOUNDS_CHECKING_1(i1);
      return m_ptr[ i1 ];
    }

  //----------------------------------
  // Required constructors and assignment operators:

  ArrayFortran() : m_ptr(NULL) { Helper::zero( m_stride ); }

  ArrayFortran( const ArrayFortran & rhs )
    : m_ptr( rhs.m_ptr ) { Helper::copy( m_stride , rhs.m_stride ); }

  ArrayFortran & operator = ( const ArrayFortran & rhs )
    {
      m_ptr = rhs.m_ptr ;
      Helper::copy( m_stride , rhs.m_stride );
      return *this ;
    }

  //----------------------------------
  // Copy a reverse-map array

  typedef typename
    ArrayReverse<Scalar,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>::natural_type
      ReverseType ;

  ArrayFortran( const ReverseType & rhs )
    : m_ptr( rhs.m_ptr ) { Helper::copy( m_stride , rhs.m_stride ); }

  ArrayFortran & operator = ( const ReverseType & rhs )
    {
      m_ptr = rhs.m_ptr ;
      Helper::copy( m_stride , rhs.m_stride );
      return *this ;
    }

  //----------------------------------
  // Truncated view

  typedef typename
    ArrayTruncate<Scalar,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>::fortran_type
      TruncatedViewType ;

  TruncatedViewType truncated_view( index_type i ) const
    {
      return TruncatedViewType( m_ptr + m_stride[ Rank - 2 ] * i ,
                                ArrayStride() , m_stride );
    }

  // "Power user" constructor, promising a properly strided input
  ArrayFortran( value_type * arg_ptr ,
                const ArrayStride & , const size_type * arg_stride )
    : m_ptr( arg_ptr ) { Helper::copy( m_stride , arg_stride ); }

  // "Power user" constructor, promising a properly strided input
  ArrayFortran( value_type * arg_ptr ,
                const ArrayStride & ,
                const size_type * arg_stride ,
                const size_type   arg_append )
    : m_ptr( arg_ptr )
    { Helper::stride_append( m_stride , arg_stride , arg_append ); }

  //----------------------------------
  // Take pointer to other compatible array types.

  template< class OtherArrayType >
  explicit ArrayFortran( const OtherArrayType & rhs )
    : m_ptr( rhs.contiguous_data() )
    { Helper::stride( m_stride , rhs ); }

  template< class OtherArrayType >
  ArrayFortran & operator = ( const OtherArrayType & rhs )
    {
      m_ptr = rhs.contiguous_data();
      Helper::stride( m_stride , rhs );
      return *this ;
    }

  //----------------------------------
  // Class specific constructors:

  ArrayFortran( value_type * arg_ptr ,
         size_type n1 , size_type n2 , size_type n3 , size_type n4 ,
         size_type n5 , size_type n6 , size_type n7 , size_type n8 )
    : m_ptr( arg_ptr )
    {
      array_bounds_check_rank_is_equal<Rank,8>();
      ARRAY_STRIDE_FORTRAN_8( m_stride, n1, n2, n3, n4, n5, n6, n7, n8 );
    }

  ArrayFortran( value_type * arg_ptr ,
         size_type n1 , size_type n2 , size_type n3 , size_type n4 ,
         size_type n5 , size_type n6 , size_type n7 )
    : m_ptr( arg_ptr )
    {
      array_bounds_check_rank_is_equal<Rank,7>();
      ARRAY_STRIDE_FORTRAN_7( m_stride, n1, n2, n3, n4, n5, n6, n7 );
    }

  ArrayFortran( value_type * arg_ptr ,
         size_type n1 , size_type n2 , size_type n3 , size_type n4 ,
         size_type n5 , size_type n6 )
    : m_ptr( arg_ptr )
    {
      array_bounds_check_rank_is_equal<Rank,6>();
      ARRAY_STRIDE_FORTRAN_6( m_stride, n1, n2, n3, n4, n5, n6 );
    }

  ArrayFortran( value_type * arg_ptr ,
         size_type n1 , size_type n2 , size_type n3 , size_type n4 ,
         size_type n5 )
    : m_ptr( arg_ptr )
    {
      array_bounds_check_rank_is_equal<Rank,5>();
      ARRAY_STRIDE_FORTRAN_5( m_stride, n1, n2, n3, n4, n5 );
    }

  ArrayFortran( value_type * arg_ptr ,
         size_type n1 , size_type n2 , size_type n3 , size_type n4 )
    : m_ptr( arg_ptr )
    {
      array_bounds_check_rank_is_equal<Rank,4>();
      ARRAY_STRIDE_FORTRAN_4( m_stride, n1, n2, n3, n4 );
    }

  ArrayFortran( value_type * arg_ptr ,
         size_type n1 , size_type n2 , size_type n3 )
    : m_ptr( arg_ptr )
    {
      array_bounds_check_rank_is_equal<Rank,3>();
      ARRAY_STRIDE_FORTRAN_3( m_stride, n1, n2, n3 );
    }

  ArrayFortran( value_type * arg_ptr ,
         size_type n1 , size_type n2 )
    : m_ptr( arg_ptr )
    {
      array_bounds_check_rank_is_equal<Rank,2>();
      ARRAY_STRIDE_FORTRAN_2( m_stride, n1, n2 );
    }

  ArrayFortran( value_type * arg_ptr ,
         size_type n1 )
    : m_ptr( arg_ptr )
    {
      array_bounds_check_rank_is_equal<Rank,1>();
      m_stride[0] = n1 ;
    }
};

//----------------------------------------------------------------------
// Specialization for arrays with rank and dimension tags as runtime info.

template< typename Scalar >
class ArrayNatural< Scalar, void,void,void,void,void,void,void,void>
{
public:

  typedef Scalar           value_type ;
  typedef int              index_type ;
  typedef size_t           size_type ;
  typedef const ArrayDimTag * tag_type ;

  enum { Natural    = true };
  enum { Reverse    = false };
  enum { Contiguous = true };

private:

  Scalar  * m_ptr ;
  size_type m_rank ;
  size_type m_stride[8];
  tag_type  m_tag[8] ;

  typedef ArrayHelper<8> Helper ;

  template< typename ,
            class , class , class , class ,
            class , class , class , class >
  friend class ArrayNatural ;

  template< typename ,
            class , class , class , class ,
            class , class , class , class >
  friend class ArrayFortran ;

public:

  //----------------------------------

  unsigned rank() const { return m_rank ; }
  bool natural() const { return Natural ; }
  bool reverse() const { return Reverse ; }
  bool contiguous() const { return Contiguous ; }

  tag_type tag( unsigned ord ) const
    {
      array_bounds_checking( m_rank , ord );
      return m_tag[ ( m_rank - 1 ) - ord ];
    }

  //----------------------------------

  size_type size() const { return m_stride[ m_rank - 1 ]; }

  size_type dimension( index_type ord ) const
    {
      array_bounds_checking( m_rank , ord );
      const int i = ( m_rank - 1 ) - ord ;
      return i ? m_stride[i] / m_stride[i-1] : m_stride[i] ;
    }

  void dimensions( std::vector<size_type> & n )
    {
      n.resize( m_rank );
      n[m_rank-1] = m_stride[0] ;
      for ( int i = 1 ; i < m_rank ; ++i ) {
        n[ ( m_rank - 1 ) - i ] = m_stride[i] / m_stride[i-1] ;
      }
    }
  
  //----------------------------------
  /** \brief Access member data */
  value_type * contiguous_data() const { return m_ptr ; }

  /** \brief Access member via full ordering of members. */
  value_type & operator[]( size_type i ) const
    {
      ARRAY_BOUNDS_CHECKING_ORDINAL( size() , i );
      return m_ptr[ i ];
    }

  /** \brief Access member via Rank 8 multi-index */
  value_type & operator()( index_type i1 , index_type i2 ,
                           index_type i3 , index_type i4 ,
                           index_type i5 , index_type i6 ,
                           index_type i7 , index_type i8 ) const
    {
      ARRAY_BOUNDS_CHECKING_8(i1,i2,i3,i4,i5,i6,i7,i8);
      return m_ptr[ ARRAY_INDEX_NATURAL_8(m_stride,i1,i2,i3,i4,i5,i6,i7,i8) ];
    }

  /** \brief Access member via Rank 7 multi-index */
  value_type & operator()( index_type i1 , index_type i2 ,
                           index_type i3 , index_type i4 ,
                           index_type i5 , index_type i6 ,
                           index_type i7 ) const
    {
      ARRAY_BOUNDS_CHECKING_7(i1,i2,i3,i4,i5,i6,i7);
      return m_ptr[ ARRAY_INDEX_NATURAL_7(m_stride,i1,i2,i3,i4,i5,i6,i7) ];
    }

  /** \brief Access member via Rank 6 multi-index */
  value_type & operator()( index_type i1 , index_type i2 ,
                           index_type i3 , index_type i4 ,
                           index_type i5 , index_type i6 ) const
    {
      ARRAY_BOUNDS_CHECKING_6(i1,i2,i3,i4,i5,i6);
      return m_ptr[ ARRAY_INDEX_NATURAL_6(m_stride,i1,i2,i3,i4,i5,i6) ];
    }

  /** \brief Access member via Rank 5 multi-index */
  value_type & operator()( index_type i1 , index_type i2 ,
                           index_type i3 , index_type i4 ,
                           index_type i5 ) const
    {
      ARRAY_BOUNDS_CHECKING_5(i1,i2,i3,i4,i5);
      return m_ptr[ ARRAY_INDEX_NATURAL_5(m_stride,i1,i2,i3,i4,i5) ];
    }

  /** \brief Access member via Rank 4 multi-index */
  value_type & operator()( index_type i1 , index_type i2 ,
                           index_type i3 , index_type i4 ) const
    {
      ARRAY_BOUNDS_CHECKING_4(i1,i2,i3,i4);
      return m_ptr[ ARRAY_INDEX_NATURAL_4(m_stride,i1,i2,i3,i4) ];
    }

  /** \brief Access member via Rank 3 multi-index */
  value_type & operator()( index_type i1 , index_type i2 ,
                           index_type i3 ) const
    {
      ARRAY_BOUNDS_CHECKING_3(i1,i2,i3);
      return m_ptr[ ARRAY_INDEX_NATURAL_3(m_stride,i1,i2,i3) ];
    }

  /** \brief Access member via Rank 2 multi-index */
  value_type & operator()( index_type i1 , index_type i2 ) const
    {
      ARRAY_BOUNDS_CHECKING_2(i1,i2);
      return m_ptr[ ARRAY_INDEX_NATURAL_2(m_stride,i1,i2) ];
    }

  /** \brief Access member via Rank 1 multi-index */
  value_type & operator()( index_type i1 ) const
    {
      ARRAY_BOUNDS_CHECKING_1(i1);
      return m_ptr[ i1 ];
    }

  //----------------------------------
  /** \brief Truncate view of the array */

  typedef ArrayNatural<Scalar> TruncateViewType ;

  TruncateViewType truncated_view( index_type i ) const
    {
      return TruncateViewType( m_ptr + m_stride[ m_rank - 2 ] * i ,
                               m_rank - 1 ,
                               ArrayStride() ,
                               m_stride ,
                               m_tag );
    }

  ArrayNatural( value_type      * arg_ptr ,
                size_type         arg_rank ,
                const ArrayStride & ,
                const size_type * arg_stride ,
                tag_type        * arg_tag )
    : m_ptr( arg_ptr ), m_rank( arg_rank )
    {
      Helper::copy( m_stride , arg_stride );
      Helper::copy( m_tag ,    arg_tag );
    }

  //----------------------------------
  // Required constructors and assignment operators:

  ArrayNatural()
    : m_ptr(NULL), m_rank(0)
    {
      Helper::zero( m_stride );
      Helper::zero( m_tag );
    }

  ArrayNatural( const ArrayNatural & rhs )
    : m_ptr( rhs.m_ptr ), m_rank( rhs.m_rank )
    {
      Helper::copy( m_stride , rhs.m_stride );
      Helper::copy( m_tag ,    rhs.m_tag );
    }

  ArrayNatural & operator = ( const ArrayNatural & rhs )
    {
      m_ptr  = rhs.m_ptr ;
      m_rank = rhs.m_rank ;
      Helper::copy( m_stride , rhs.m_stride );
      Helper::copy( m_tag ,    rhs.m_tag );
      return *this ;
    }

  //----------------------------------

  ArrayNatural( const ArrayFortran<Scalar> & rhs )
    : m_ptr( rhs.m_ptr ), m_rank( rhs.m_rank )
    {
      Helper::copy( m_stride , rhs.m_stride );
      Helper::copy( m_tag ,    rhs.m_tag );
    }

  ArrayNatural & operator = ( const ArrayFortran<Scalar> & rhs )
    {
      m_ptr  = rhs.m_ptr ;
      m_rank = rhs.m_rank ;
      Helper::copy( m_stride , rhs.m_stride );
      Helper::copy( m_tag ,    rhs.m_tag );
      return *this ;
    }

  //----------------------------------
  // Take pointer to other compatible array types.

  template< class OtherArrayType >
  ArrayNatural( const OtherArrayType & rhs )
    : m_ptr( rhs.contiguous_data() ), m_rank( rhs.rank() )
    {
      Helper::stride( m_stride , rhs );
      Helper::tags( m_tag , rhs );
    }

  template< class OtherArrayType >
  ArrayNatural & operator = ( const OtherArrayType & rhs )
    {
      m_ptr  = rhs.contiguous_data();
      m_rank = rhs.rank();
      Helper::stride( m_stride , rhs );
      Helper::tags(   m_tag ,    rhs );
      return *this ;
    }

  //----------------------------------
  // Class specific constructors:

  ArrayNatural( value_type * arg_ptr ,
                const std::vector<size_type> & arg_size ,
                const std::vector<tag_type>  & arg_tag )
    : m_ptr( arg_ptr ),
      m_rank( arg_size.size() == arg_tag.size() ? arg_size.size() : 0 )
    {
      array_stride_from_natural_sizes( m_rank , m_stride , & arg_size[0] );
      array_stride_natural_tag( m_rank , m_tag , & arg_tag[0] );
    }
};

//----------------------------------------------------------------------
// Specialization for arrays with rank and dimension tags as runtime info.

template< typename Scalar >
class ArrayFortran< Scalar, void,void,void,void,void,void,void,void>
{
public:

  typedef Scalar              value_type ;
  typedef int                 index_type ;
  typedef size_t              size_type ;
  typedef const ArrayDimTag * tag_type ;

  enum { Natural    = false };
  enum { Reverse    = true };
  enum { Contiguous = true };

private:

  Scalar  * m_ptr ;
  size_type m_rank ;
  size_type m_stride[8];
  tag_type  m_tag[8] ;

  typedef ArrayHelper<8> Helper ;

  template< typename ,
            class , class , class , class ,
            class , class , class , class >
  friend class ArrayNatural ;

  template< typename ,
            class , class , class , class ,
            class , class , class , class >
  friend class ArrayFortran ;

public:

  //----------------------------------

  unsigned rank() const { return m_rank ; }
  bool natural() const { return Natural ; }
  bool reverse() const { return Reverse ; }
  bool contiguous() const { return Contiguous ; }

  tag_type tag( unsigned ord ) const
    {
      array_bounds_checking( m_rank , ord );
      return m_tag[ ord ];
    }

  //----------------------------------

  size_type size() const { return m_stride[ m_rank - 1 ]; }

  size_type dimension( index_type ord ) const
    {
      array_bounds_checking( m_rank , ord );
      return ord ? m_stride[ord] / m_stride[ord-1] : m_stride[ord] ;
    }

  void dimensions( std::vector<size_type> & n )
    {
      n.resize( m_rank );
      n[0] = m_stride[0] ;
      for ( int i = 1 ; i < m_rank ; ++i ) {
        n[i] = m_stride[i] / m_stride[i-1] ;
      }
    }
  
  //----------------------------------
  /** \brief Access member data */
  value_type * contiguous_data() const { return m_ptr ; }

  /** \brief Access member via full ordering of members. */
  value_type & operator[]( size_type i ) const
    {
      ARRAY_BOUNDS_CHECKING_ORDINAL( size() , i );
      return m_ptr[ i ];
    }

  /** \brief Access member via Rank 8 multi-index */
  value_type & operator()( index_type i1 , index_type i2 ,
                           index_type i3 , index_type i4 ,
                           index_type i5 , index_type i6 ,
                           index_type i7 , index_type i8 ) const
    {
      ARRAY_BOUNDS_CHECKING_8(i1,i2,i3,i4,i5,i6,i7,i8);
      return m_ptr[ ARRAY_INDEX_FORTRAN_8(m_stride,i1,i2,i3,i4,i5,i6,i7,i8) ];
    }

  /** \brief Access member via Rank 7 multi-index */
  value_type & operator()( index_type i1 , index_type i2 ,
                           index_type i3 , index_type i4 ,
                           index_type i5 , index_type i6 ,
                           index_type i7 ) const
    {
      ARRAY_BOUNDS_CHECKING_7(i1,i2,i3,i4,i5,i6,i7);
      return m_ptr[ ARRAY_INDEX_FORTRAN_7(m_stride,i1,i2,i3,i4,i5,i6,i7) ];
    }

  /** \brief Access member via Rank 6 multi-index */
  value_type & operator()( index_type i1 , index_type i2 ,
                           index_type i3 , index_type i4 ,
                           index_type i5 , index_type i6 ) const
    {
      ARRAY_BOUNDS_CHECKING_6(i1,i2,i3,i4,i5,i6);
      return m_ptr[ ARRAY_INDEX_FORTRAN_6(m_stride,i1,i2,i3,i4,i5,i6) ];
    }

  /** \brief Access member via Rank 5 multi-index */
  value_type & operator()( index_type i1 , index_type i2 ,
                           index_type i3 , index_type i4 ,
                           index_type i5 ) const
    {
      ARRAY_BOUNDS_CHECKING_5(i1,i2,i3,i4,i5);
      return m_ptr[ ARRAY_INDEX_FORTRAN_5(m_stride,i1,i2,i3,i4,i5) ];
    }

  /** \brief Access member via Rank 4 multi-index */
  value_type & operator()( index_type i1 , index_type i2 ,
                           index_type i3 , index_type i4 ) const
    {
      ARRAY_BOUNDS_CHECKING_4(i1,i2,i3,i4);
      return m_ptr[ ARRAY_INDEX_FORTRAN_4(m_stride,i1,i2,i3,i4) ];
    }

  /** \brief Access member via Rank 3 multi-index */
  value_type & operator()( index_type i1 , index_type i2 ,
                           index_type i3 ) const
    {
      ARRAY_BOUNDS_CHECKING_3(i1,i2,i3);
      return m_ptr[ ARRAY_INDEX_FORTRAN_3(m_stride,i1,i2,i3) ];
    }

  /** \brief Access member via Rank 2 multi-index */
  value_type & operator()( index_type i1 , index_type i2 ) const
    {
      ARRAY_BOUNDS_CHECKING_2(i1,i2);
      return m_ptr[ ARRAY_INDEX_FORTRAN_2(m_stride,i1,i2) ];
    }

  /** \brief Access member via Rank 1 multi-index */
  value_type & operator()( index_type i1 ) const
    {
      ARRAY_BOUNDS_CHECKING_1(i1);
      return m_ptr[ i1 ];
    }

  //----------------------------------
  /** \brief Truncate view of the array */

  typedef ArrayFortran<Scalar> TruncateViewType ;

  TruncateViewType truncated_view( index_type i ) const
    {
      return TruncateViewType( m_ptr + m_stride[ m_rank - 2 ] * i ,
                               m_rank - 1 ,
                               ArrayStride() ,
                               m_stride ,
                               m_tag );
    }

  ArrayFortran( value_type      * arg_ptr ,
                size_type         arg_rank ,
                const ArrayStride & ,
                const size_type * arg_stride ,
                tag_type        * arg_tag )
    : m_ptr( arg_ptr ), m_rank( arg_rank )
    {
      Helper::copy( m_stride , arg_stride );
      Helper::copy( m_tag ,    arg_tag );
    }

  //----------------------------------
  // Required constructors and assignment operators:

  ArrayFortran()
    : m_ptr(NULL), m_rank(0)
    {
      Helper::zero( m_stride );
      Helper::zero( m_tag );
    }

  ArrayFortran( const ArrayFortran & rhs )
    : m_ptr( rhs.m_ptr ), m_rank( rhs.m_rank )
    {
      Helper::copy( m_stride , rhs.m_stride );
      Helper::copy( m_tag ,    rhs.m_tag );
    }

  ArrayFortran & operator = ( const ArrayFortran & rhs )
    {
      m_ptr  = rhs.m_ptr ;
      m_rank = rhs.m_rank ;
      Helper::copy( m_stride , rhs.m_stride );
      Helper::copy( m_tag ,    rhs.m_tag );
      return *this ;
    }

  ArrayFortran( const ArrayNatural<Scalar> & rhs )
    : m_ptr( rhs.m_ptr ), m_rank( rhs.m_rank )
    {
      Helper::copy( m_stride , rhs.m_stride );
      Helper::copy( m_tag ,    rhs.m_tag );
    }

  ArrayFortran & operator = ( const ArrayNatural<Scalar> & rhs )
    {
      m_ptr  = rhs.m_ptr ;
      m_rank = rhs.m_rank ;
      Helper::copy( m_stride , rhs.m_stride );
      Helper::copy( m_tag ,    rhs.m_tag );
      return *this ;
    }

  //----------------------------------
  // Take pointer to other compatible array types.

  template< class OtherArrayType >
  ArrayFortran( const OtherArrayType & rhs )
    : m_ptr( rhs.contiguous_data() ), m_rank( rhs.rank() )
    {
      Helper::stride( m_stride , rhs );
      Helper::tags( m_tag , rhs );
    }

  template< class OtherArrayType >
  ArrayFortran & operator = ( const OtherArrayType & rhs )
    {
      m_ptr = rhs.contiguous_data();
      m_rank = rhs.m_rank();
      Helper::stride( m_stride , rhs );
      Helper::tags( m_tag , rhs );
      return *this ;
    }

  //----------------------------------
  // Class specific constructors:

  ArrayFortran( value_type * arg_ptr ,
                const std::vector<size_type> & arg_size ,
                const std::vector<tag_type>  & arg_tag )
    : m_ptr( arg_ptr ),
      m_rank( arg_size.size() == arg_tag.size() ? arg_size.size() : 0 )
    {
      array_stride_from_fortran_sizes( m_rank , m_stride , & arg_size[0] );
      array_stride_fortran_tag( m_rank , m_tag , & arg_tag[0] );
    }
};

//----------------------------------------------------------------------

/** \endcond */

}

#endif

