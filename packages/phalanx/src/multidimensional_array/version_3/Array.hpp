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

#ifdef ARRAY_BOUNDS_CHECKING
#define ARRAY_CHECK( X ) X
#else
#define ARRAY_CHECK( X )
#endif

#include <vector>
#include <string>
#include <SimpleArrayOps.hpp>

namespace phdmesh {

/**
 * \defgroup mdarray_module
 * \author H. Carter Edwards  <hcedwar@sandia.gov>
 * \date   June 2008
 */

//----------------------------------------------------------------------
/** \brief  Define <b> Natural </b> (C-language) or
 *          <b> Fortran </b> ordering of array dimensions.
 *          A RankZero array does not have an ordering.
 *  \ingroup mdarray_module
 */
enum ArrayOrder {
  /** \brief  Use the Natural or C-language ordering for multi-dimensions
   *          where the right-most dimension is stride-one.
   */
  NaturalOrder ,

  /** \brief  Use the Reverse or Fortran-language ordering for multi-dimensions
   *          where the left-most dimension is stride-one.
   */
  FortranOrder ,

  /** \brief  Special tag to indicate that an array specification has
   *          degenerated to rank-zero, i.e. is no longer an array.
   */
  RankZero
};

template< typename Scalar , ArrayOrder Order , 
          class Tag1 = void , class Tag2 = void ,
          class Tag3 = void , class Tag4 = void ,
          class Tag5 = void , class Tag6 = void ,
          class Tag7 = void , class Tag8 = void >
class Array ;

//----------------------------------------------------------------------

template< class ArrayType , class Tag > struct ArrayAppend ;

//----------------------------------------------------------------------
/** \class  ArrayDimTag
 *  \brief  Virtual base class for array dimension tags supplied to
 *          the Array template class.
 *  \ingroup mdarray_module 
 *  \sa Array
 *
 *  A derived array dimension tag class must provide the
 *  <b> name </b> method and <b> tag </b> singleton method
 *  as in the following example.
 *  <PRE>
 *  struct MyTag : public phdmesh::ArrayDimTag {
 *    const char * name() const ;
 *    static const MyTag & tag();
 *  };
 *  </PRE>
 *  An example implementation of these methods is as follows.
 *  <PRE>
 *  const char * MyTag::name() const
 *  { static const char my_name[] = "MyTag" ; return my_name ; }
 *
 *  const MyTag & MyTag::tag()
 *  { static const MyTag my_tag ; return my_tag ; }
 *  </PRE>
 */
struct ArrayDimTag {

  /** \brief Name of the tag, typically the name of the derived class. */
  virtual const char * name() const = 0 ;

  /** \brief  Given a dimension and index produce a string for output.
   *          Default to converting <b> index </b> to a string.
   */
  virtual std::string to_string( unsigned dimension ,
                                 unsigned index ) const ;

  /** \brief Given a dimension and input strige produce an index.
   *          Default to converting <b> label </b> to an integer.
   */
  virtual unsigned to_index( unsigned dimension ,
                             const std::string & label ) const ; 
 
protected:
  virtual ~ArrayDimTag();
  ArrayDimTag() {}
private:
  ArrayDimTag( const ArrayDimTag & );
  ArrayDimTag & operator = ( const ArrayDimTag & );
};

/** \class  ArrayDimension
 *  \brief  An anonymous array dimension tag,
 *          which is NOT the recommended usage.
 *  \ingroup mdarray_module
 */
struct ArrayDimension : public ArrayDimTag {

  const char * name() const ;

  static const ArrayDimension & tag();

private:
  ~ArrayDimension();
  ArrayDimension() {}
  ArrayDimension( const ArrayDimension & );
  ArrayDimension & operator = ( const ArrayDimension & );
};

} // namespace phdmesh

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#include <ArrayPrivate.hpp>

namespace phdmesh {

//----------------------------------------------------------------------
/** \class  Array
 *  \brief  The <b> preferred </b> multi-dimensional Array interface
 *          with <b> compile-time </b> user-defined dimension ordinates.
 *  \ingroup mdarray_module
 *  \nosubgrouping
 *
 *  \param Scalar  The "plain old data" type of the array's member data.
 *  \param array_order An <b> ArrayOrder </b> value that specifies whether to
 *                     use Natural (a.k.a. C-language) or Fortran ordering 
 *                     for the multi-dimensions and multi-indices.
 *  \param Tag#  The <b> Tag# </b> template parameters document the
 *               user-defiend purpose of each dimension ordinate.
 *               The <b> Rank </b> of the array (i.e. the number of dimensions)
 *               is the number of user-defined dimension tags, up to eight.
 *               A user-defined dimension <b> Tag# </b> must be derived from
 *               the <b> ArrayDimTag </b> template class.
 *
 *  \sa ArrayDimTag ArrayOrder
 */
template< typename Scalar , ArrayOrder array_order ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
class Array
{
private:
  typedef
    Array<void,array_order,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>
      md_type ;
public:
  /** \name Array Attributes
   *  \{
   */

  /** \brief  Type of member data. */
  typedef Scalar  value_type ;

  /** \brief  Type for sizes. */
  typedef unsigned size_type ;

  /** \brief  Type of runtime dimension tags. */
  typedef const ArrayDimTag * tag_type ;

  //----------------------------------

  /** \brief  Rank of the array is the number of non-void dimension tags. */
  enum { Rank = md_type::Rank };

  /** \brief  If the multidimension follows the natural ordering */
  enum { Natural = NaturalOrder == array_order };

  /** \brief  If the multidimension follows the reverse (Fortran) ordering */
  enum { Reverse = FortranOrder == array_order };

  /** \brief  If the member data storage is contiguous */
  enum { Contiguous = true };

  /** \brief  Rank of the array is the number of non-void dimension tags. */
  unsigned rank()   const { return Rank ; }

  /** \brief  If the multidimension follows the natural ordering */
  bool natural()    const { return Natural ; }

  /** \brief  If the multidimension follows the reverse (Fortran) ordering */
  bool reverse()    const { return Reverse ; }

  /** \brief  If the member data storage is contiguous */
  bool contiguous() const { return Contiguous ; }

  //----------------------------------

#ifndef DOXYGEN_COMPILE
  /** \brief  Access the dimension tag-type for a given ordinate. */
  template < unsigned ordinate >
  struct Tag { typedef typename ArrayTagAt<Array,ordinate>::type type ; };
#endif

  /** \brief  Access the dimension tag-singleton for a given ordinate. */
  tag_type tag( const unsigned ordinate ) const
    {
      array_check_ordinal( Rank , ordinate );
      return
        array_dim_tags<Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>()[ordinate];
    }

  //----------------------------------
  /** \brief  Dimension of the given ordinate. */
  template < unsigned ordinate > unsigned dimension() const
    {
      array_check_ordinal_is_less<ordinate,Rank>();
      return ArrayStrideDim<array_order,Rank,ordinate>::dimension(m_stride);
    }

  /** \brief  Dimension of the given ordinate. */
  unsigned dimension( const unsigned ordinate ) const
    {
      array_check_ordinal( Rank , ordinate );
      return ArrayStrideDim<array_order,Rank>::dimension(m_stride,ordinate);
    }

  /** \brief  Dimensions of all ordinates. */
  void dimensions( std::vector<unsigned> & n )
    {
      n.resize( Rank );
      for ( unsigned i = 0 ; i < Rank ; ++i ) { n[i] = dimension(i); }
    }

  /** \brief  Total number of member data items. */
  size_type size() const { return m_stride[ Rank - 1 ]; }

  /** \} */
  //----------------------------------
  /** \name Member data access operators
   *  \{
   */

  /** \brief  Subarray type that removes the slowest striding dimension
   *          (first natural or last fortran ordinate).
   */
  typedef typename ArrayTruncate<Array>::type TruncateType ;

  /** \brief  Generate a subarray view of the array with the
   *          slowest striding ordinate offset by <b> i </b>
   *          and removed.
   */
  TruncateType truncate( const unsigned i ) const
    {
      TruncateType tmp ;
      tmp.m_ptr = m_ptr + i * ( 1 < Rank ? m_stride[ Rank - 2 ] : 1 );
      Copy<Rank-1>( tmp.m_stride , m_stride );
      return tmp ;
    }

  //----------------------------------
  /** \brief Pointer to contiguous block of member data. */
  value_type * contiguous_data() const { return m_ptr ; }

  /** \brief Access member via offset into contiguous block. */
  value_type & operator[]( size_type i ) const
    {
      ARRAY_CHECK( array_check_offset(size(),i) );
      return m_ptr[ i ];
    }

  /** \brief Access member of a Rank 8 array */
  value_type & operator()( const unsigned i1 , const unsigned i2 ,
                           const unsigned i3 , const unsigned i4 ,
                           const unsigned i5 , const unsigned i6 ,
                           const unsigned i7 , const unsigned i8 ) const
    { return m_ptr[
        array_offset<array_order,Rank>(m_stride,i1,i2,i3,i4,i5,i6,i7,i8) ];
    }

  /** \brief Access member of a Rank 7 array */
  value_type & operator()( const unsigned i1 , const unsigned i2 ,
                           const unsigned i3 , const unsigned i4 ,
                           const unsigned i5 , const unsigned i6 ,
                           const unsigned i7 ) const
    { return m_ptr[
        array_offset<array_order,Rank>(m_stride,i1,i2,i3,i4,i5,i6,i7) ];
    }

  /** \brief Access member of a Rank 6 array */
  value_type & operator()( const unsigned i1 , const unsigned i2 ,
                           const unsigned i3 , const unsigned i4 ,
                           const unsigned i5 , const unsigned i6 ) const
    { return m_ptr[
        array_offset<array_order,Rank>(m_stride,i1,i2,i3,i4,i5,i6) ];
    }

  /** \brief Access member of a Rank 5 array */
  value_type & operator()( const unsigned i1 , const unsigned i2 ,
                           const unsigned i3 , const unsigned i4 ,
                           const unsigned i5 ) const
    { return m_ptr[ array_offset<array_order,Rank>(m_stride,i1,i2,i3,i4,i5) ]; }

  /** \brief Access member of a Rank 4 array */
  value_type & operator()( const unsigned i1 , const unsigned i2 ,
                           const unsigned i3 , const unsigned i4 ) const
    { return m_ptr[ array_offset<array_order,Rank>(m_stride,i1,i2,i3,i4) ]; }

  /** \brief Access member of a Rank 3 array */
  value_type & operator()( const unsigned i1 , const unsigned i2 ,
                           const unsigned i3 ) const
    { return m_ptr[ array_offset<array_order,Rank>(m_stride,i1,i2,i3) ]; }

  /** \brief Access member of a Rank 2 array */
  value_type & operator()( const unsigned i1 , const unsigned i2 ) const
    { return m_ptr[ array_offset<array_order,Rank>(m_stride,i1,i2) ]; }

  /** \brief Access member of a Rank 1 array */
  value_type & operator()( const unsigned i1 ) const
    { return m_ptr[ array_offset<array_order,Rank>(m_stride,i1) ]; }

  /** \} */
  //----------------------------------
  /** \name Constructors and Assignment Operators
   * \{
   */

  /** \brief  The compatible multidimensional array with
   *          reversed multi-index ordering and dimension tags.
   */
  typedef typename ArrayReverse< Array >::type ReverseType ;

  /** \brief Default constructor */
  Array() : m_ptr(NULL) { Copy<Rank>( m_stride , (size_type) 0 ); }

  /** \brief Copy constructor */
  Array( const Array & rhs )
    : m_ptr( rhs.m_ptr ) { Copy<Rank>( m_stride , rhs.m_stride ); }

  /** \brief Assignment operator */
  Array & operator = ( const Array & rhs )
    {
      m_ptr = rhs.m_ptr ;
      Copy<Rank>( m_stride , rhs.m_stride );
      return *this ;
    }

  /** \brief Copy constructor for compatible reverse type. */
  Array( const ReverseType & rhs )
    : m_ptr( rhs.m_ptr ) { Copy<Rank>( m_stride , rhs.m_stride ); }

  /** \brief Assignment operator for compatible reverse type. */
  Array & operator = ( const ReverseType & rhs )
    {
      m_ptr = rhs.m_ptr ;
      Copy<Rank>( m_stride , rhs.m_stride );
      return *this ;
    }

  /** \brief Construct with array of dimensions. */
  Array( value_type * arg_ptr , const unsigned * const dims )
    : m_ptr( arg_ptr ) { md_type::assign( m_stride , dims ); }

  /** \brief  Construct a Rank 8 array */
  Array( value_type * arg_ptr ,
         const unsigned n1 , const unsigned n2 ,
         const unsigned n3 , const unsigned n4 ,
         const unsigned n5 , const unsigned n6 ,
         const unsigned n7 , const unsigned n8 )
    : m_ptr( arg_ptr )
    { md_type::assign( m_stride , n1 , n2 , n3 , n4 , n5 , n6 , n7 , n8 ); }

  /** \brief  Construct a Rank 7..8 array; use Tag#::Size for defaults.
   *          The input dimensions are the 7 slowest strides.
   */
  Array( value_type * arg_ptr ,
         const unsigned n1 , const unsigned n2 ,
         const unsigned n3 , const unsigned n4 ,
         const unsigned n5 , const unsigned n6 ,
         const unsigned n7 )
    : m_ptr( arg_ptr )
    { md_type::assign( m_stride , n1 , n2 , n3 , n4 , n5 , n6 , n7 ); }

  /** \brief  Construct a Rank 6..8 array; use Tag#::Size for defaults.
   *          The input dimensions are the 6 slowest strides.
   */
  Array( value_type * arg_ptr ,
         const unsigned n1 , const unsigned n2 ,
         const unsigned n3 , const unsigned n4 ,
         const unsigned n5 , const unsigned n6 )
    : m_ptr( arg_ptr )
    { md_type::assign( m_stride , n1 , n2 , n3 , n4 , n5 , n6 ); }

  /** \brief  Construct a Rank 5..8 array; use Tag#::Size for defaults.
   *          The input dimensions are the 5 slowest strides.
   */
  Array( value_type * arg_ptr ,
         const unsigned n1 , const unsigned n2 ,
         const unsigned n3 , const unsigned n4 ,
         const unsigned n5 )
    : m_ptr( arg_ptr )
    { md_type::assign( m_stride , n1 , n2 , n3 , n4 , n5 ); }

  /** \brief  Construct a Rank 4..8 array; use Tag#::Size for defaults.
   *          The input dimensions are the 4 slowest strides.
   */
  Array( value_type * arg_ptr ,
         const unsigned n1 , const unsigned n2 ,
         const unsigned n3 , const unsigned n4 )
    : m_ptr( arg_ptr )
    { md_type::assign( m_stride , n1 , n2 , n3 , n4 ); }

  /** \brief  Construct a Rank 3..8 array; use Tag#::Size for defaults.
   *          The input dimensions are the 3 slowest strides.
   */
  Array( value_type * arg_ptr ,
         const unsigned n1 , const unsigned n2 ,
         const unsigned n3 )
    : m_ptr( arg_ptr )
    { md_type::assign( m_stride , n1 , n2 , n3 ); }

  /** \brief  Construct a Rank 2..8 array; use Tag#::Size for defaults.
   *          The input dimensions are the 2 slowest strides.
   */
  Array( value_type * arg_ptr , const unsigned n1 , const unsigned n2 )
    : m_ptr( arg_ptr ) { md_type::assign( m_stride , n1 , n2 ); }

  /** \brief  Construct a Rank 1..8 array; use Tag#::Size for defaults.
   *          The input dimension is the slowest stride.
   */
  Array( value_type * arg_ptr , const unsigned n1 )
    : m_ptr( arg_ptr ) { md_type::assign( m_stride , n1 ); }

  /** \brief  Construct a Rank 1..8 array; use Tag#::Size for defaults. */
  Array( value_type * arg_ptr )
    : m_ptr( arg_ptr ) { md_type::assign( m_stride ); }

  /** \} */
protected:

  value_type * m_ptr ;
  size_type    m_stride[ Rank ];

  template< typename , ArrayOrder ,
            class , class , class , class ,
            class , class , class , class >
  friend class phdmesh::Array ;
};

//----------------------------------------------------------------------

#ifndef DOXYGEN_COMPILE

/** \brief  Specialization for an array with Rank = 0.
 *  \ingroup mdarray_module
 */
template< typename Scalar >
class Array<Scalar,RankZero,void,void,void,void,void,void,void,void>
{
public:

  typedef Scalar              value_type ;
  typedef unsigned            size_type ;
  typedef const ArrayDimTag * tag_type ;

  //----------------------------------

  enum { Rank       = 0 };
  enum { Natural    = false };
  enum { Reverse    = false };
  enum { Contiguous = true };

  unsigned rank()   const { return Rank ; }
  bool natural()    const { return Natural ; }
  bool reverse()    const { return Reverse ; }
  bool contiguous() const { return Contiguous ; }

  //----------------------------------

  /** \brief  Total number of member data items. */
  size_type size() const { return m_ptr ? 1 : 0 ; }

  //----------------------------------
  /** \brief Pointer to contiguous block of member data. */
  value_type * contiguous_data() const { return m_ptr ; }

  /** \brief Access member via Rank 0 multi-index */
  value_type & operator()() const { return *m_ptr ; }

  //----------------------------------
  // Required constructors and assignment operators:

  Array() : m_ptr(NULL) {}

  Array( const Array & rhs ) : m_ptr( rhs.m_ptr ) {}

  Array & operator = ( const Array & rhs )
    { m_ptr = rhs.m_ptr ; return *this ; }

  //----------------------------------
  // Class specific constructors:

  Array( value_type * arg_ptr ) : m_ptr( arg_ptr ) {}

protected:

  value_type * m_ptr ;

  template< typename , ArrayOrder ,
            class , class , class , class ,
            class , class , class , class >
  friend class phdmesh::Array ;
};

#endif /* DOXYGEN_COMPILE */
//----------------------------------------------------------------------
//----------------------------------------------------------------------
/** \brief  The <b> not-preferred </b> multi-dimensional Array interface
 *          with <b> runtime </b> user-defined dimension ordinates.
 *          Typically used when runtime-polymorphic arrays are passed to
 *          functions.
 *          
 *  \ingroup mdarray_module
 *  \nosubgrouping
 *
 *  \param Scalar  The "plain old data" type of the array's member data.
 *  \param array_order An <b> ArrayOrder </b> value that specifies whether to
 *                     use Natural (a.k.a. C-language) or Fortran ordering 
 *                     for the multi-dimensions and multi-indices.
 */
template< typename Scalar , ArrayOrder array_order >
class Array<Scalar,array_order,void,void,void,void,void,void,void,void>
{
public:
  /** \name Array Attributes
   *  \{
   */

  /** \brief  Type of member data. */
  typedef Scalar  value_type ;

  /** \brief  Type for sizes. */
  typedef unsigned size_type ;

  /** \brief  Type of runtime dimension tags. */
  typedef const ArrayDimTag * tag_type ;

  //----------------------------------

  /** \brief  If the multidimension follows the natural ordering */
  enum { Natural = NaturalOrder == array_order };

  /** \brief  If the multidimension follows the reverse (Fortran) ordering */
  enum { Reverse = FortranOrder == array_order };

  /** \brief  If the member data storage is contiguous */
  enum { Contiguous = true };

  /** \brief  Rank of the array is the number of non-void dimension tags. */
  unsigned rank()   const { return m_rank ; }

  /** \brief  If the multidimension follows the natural ordering */
  bool natural()    const { return Natural ; }

  /** \brief  If the multidimension follows the reverse (Fortran) ordering */
  bool reverse()    const { return Reverse ; }

  /** \brief  If the member data storage is contiguous */
  bool contiguous() const { return Contiguous ; }

  //----------------------------------

  /** \brief  Access the dimension tag-singleton for a given ordinate. */
  tag_type tag( const unsigned ordinal ) const
    {
      array_check_ordinal( m_rank , ordinal );
      const int i = Natural ? ( m_rank - 1 ) - ordinal : ordinal ;
      return m_tag[i];
    }

  //----------------------------------

  /** \brief  Dimension of the given ordinate. */
  unsigned dimension( const unsigned ordinal ) const
    {
      array_check_ordinal( m_rank , ordinal );
      const int i = Natural ? ( m_rank - 1 ) - ordinal : ordinal ;
      return i ? m_stride[i] / m_stride[i-1] : m_stride[i] ;
    }

  /** \brief  Dimension of all ordinate. */
  void dimensions( std::vector<unsigned> & n )
    {
      n.resize( m_rank );
      for ( unsigned i = 0 ; i < m_rank ; ++i ) { n[i] = dimension(i); }
    }

  /** \brief  Total number of data items. */
  size_type size() const { return m_stride[ m_rank - 1 ]; }

  /** \} */
  //----------------------------------
  /** \name Member data access operators
   *  \{
   */

  /** \brief  Generate a subarray view of the array with the slowest
   *          striding ordinate offset by <b> i </b> and removed.
   */
  Array truncate( const unsigned i ) const
    {
      Array tmp ;
      if ( 1 < m_rank ) {
        tmp.m_ptr  = m_ptr + m_stride[ m_rank - 2 ] * i ;
        tmp.m_rank = m_rank - 1 ;
        unsigned k ;
        for ( k = 0 ; k < m_rank - 1 ; ++k ) { tmp.m_stride[i] = m_stride[i] ; }
        for (       ; k < 8          ; ++k ) { tmp.m_stride[i] = 0 ; }
        for ( k = 0 ; k < m_rank - 1 ; ++k ) { tmp.m_tag[i] = m_tag[i] ; }
        for (       ; k < 8          ; ++k ) { tmp.m_tag[i] = NULL ; }
      }
      return tmp ;
    }

  /** \brief Pointer to contiguous block of member data. */
  value_type * contiguous_data() const { return m_ptr ; }

  /** \brief Access member via full ordering of members. */
  value_type & operator[]( size_type i ) const
    {
      ARRAY_CHECK( array_check_offset(size(),i) );
      return m_ptr[ i ];
    }

  //----------------------------------
  /** \brief Access member via Rank 8 multi-index */
  value_type & operator()( const unsigned i1 , const unsigned i2 ,
                           const unsigned i3 , const unsigned i4 ,
                           const unsigned i5 , const unsigned i6 ,
                           const unsigned i7 , const unsigned i8 ) const
    {
      ARRAY_CHECK( array_check_rank( m_rank , 8 ) );
      return m_ptr[
        array_offset<array_order,8>(m_stride,i1,i2,i3,i4,i5,i6,i7,i8) ];
    }

  value_type & operator()( const unsigned i1 , const unsigned i2 ,
                           const unsigned i3 , const unsigned i4 ,
                           const unsigned i5 , const unsigned i6 ,
                           const unsigned i7 ) const
    {
      ARRAY_CHECK( array_check_rank( m_rank , 7 ) );
      return m_ptr[
        array_offset<array_order,7>(m_stride,i1,i2,i3,i4,i5,i6,i7) ];
    }

  value_type & operator()( const unsigned i1 , const unsigned i2 ,
                           const unsigned i3 , const unsigned i4 ,
                           const unsigned i5 , const unsigned i6 ) const
    {
      ARRAY_CHECK( array_check_rank( m_rank , 6 ) );
      return m_ptr[ array_offset<array_order,6>(m_stride,i1,i2,i3,i4,i5,i6) ];
    }

  value_type & operator()( const unsigned i1 , const unsigned i2 ,
                           const unsigned i3 , const unsigned i4 ,
                           const unsigned i5 ) const
    {
      ARRAY_CHECK( array_check_rank( m_rank , 5 ) );
      return m_ptr[ array_offset<array_order,5>(m_stride,i1,i2,i3,i4,i5) ];
    }

  value_type & operator()( const unsigned i1 , const unsigned i2 ,
                           const unsigned i3 , const unsigned i4 ) const
    {
      ARRAY_CHECK( array_check_rank( m_rank , 4 ) );
      return m_ptr[ array_offset<array_order,4>(m_stride,i1,i2,i3,i4) ];
    }

  value_type & operator()( const unsigned i1 , const unsigned i2 ,
                           const unsigned i3 ) const
    {
      ARRAY_CHECK( array_check_rank( m_rank , 3 ) );
      return m_ptr[ array_offset<array_order,3>(m_stride,i1,i2,i3) ];
    }

  value_type & operator()( const unsigned i1 , const unsigned i2 ) const
    {
      ARRAY_CHECK( array_check_rank( m_rank , 2 ) );
      return m_ptr[ array_offset<array_order,2>(m_stride,i1,i2) ];
    }

  value_type & operator()( const unsigned i1 ) const
    {
      ARRAY_CHECK( array_check_rank( m_rank , 1 ) );
      return m_ptr[ array_offset<array_order,1>(m_stride,i1) ];
    }

  /** \} */
  //----------------------------------
  /** \name Constructors and Assignment Operators
   * \{
   */

  typedef typename ArrayReverse< Array >::type ReverseType ;

  //----------------------------------

  Array()
    : m_ptr(NULL), m_rank(0)
    {
      Copy<8>( m_stride , (size_type) 0 );
      Copy<8>( m_tag , (tag_type) NULL );
    }

  Array( const Array & rhs )
    : m_ptr( rhs.m_ptr ), m_rank( rhs.m_rank )
    {
      Copy<8>( m_stride , rhs.m_stride );
      Copy<8>( m_tag , rhs.m_tag );
    }

  Array & operator = ( const Array & rhs )
    {
      m_ptr = rhs.m_ptr ;
      m_rank = rhs.m_rank ;
      Copy<8>( m_stride , rhs.m_stride );
      Copy<8>( m_tag , rhs.m_tag );
      return *this ;
    }

  /** \brief Copy constructor for reverse type. */
  Array( const ReverseType & rhs )
    : m_ptr( rhs.m_ptr ), m_rank( rhs.m_rank )
    {
      Copy<8>( m_stride , rhs.m_stride );
      Copy<8>( m_tag , rhs.m_tag );
    }

  /** \brief Assignment operator for reverse type. */
  Array & operator = ( const ReverseType & rhs )
    {
      m_ptr = rhs.m_ptr ;
      m_rank = rhs.m_rank ;
      Copy<8>( m_stride , rhs.m_stride );
      Copy<8>( m_tag , rhs.m_tag );
      return *this ;
    }

  //----------------------------------

  /** \brief  Copy constructor from an Array with compile-time
   *          defined rank and dimension tags.
   */
  template< ArrayOrder order ,
            class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
            class Tag5 , class Tag6 , class Tag7 , class Tag8 >
  Array(
    const Array<value_type,order,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8> & rhs )
  : m_ptr( rhs.m_ptr ), m_rank( 0 )
  {
    typedef Array<value_type,order,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8> a_t ;
    enum { inRank    = a_t::Rank };
    enum { inNatural = a_t::Natural };
    m_rank = inRank ;
    Copy< inRank >(     m_stride , rhs.m_stride );
    Copy< 8 - inRank >( m_stride + inRank , (size_type) 0 );
    unsigned i = 0 ;
    if ( inNatural ) {
      for ( ; i < inRank ; ++i ) { m_tag[i] = rhs.tag((inRank-1)-i); }
    }
    else {
      for ( ; i < inRank ; ++i ) { m_tag[i] = rhs.tag(i); }
    }
    for ( ; i < 8 ; ++i ) { m_tag[i] = NULL ; }
  }

  //----------------------------------
  // Class specific constructors:

  Array( value_type * ptr ,
         const unsigned rank ,
         const unsigned * const dims ,
         const tag_type  * const tags )
    : m_ptr( ptr ), m_rank( rank )
    {
      if ( Natural ) {
        size_type n = 1 ;
        unsigned i ;
        for ( i = 0 ; i < rank ; ++i ) { m_stride[i] = n *= dims[(rank-1)-i]; }
        for (       ; i < 8    ; ++i ) { m_stride[i] = 0 ; }
        for ( i = 0 ; i < rank ; ++i ) { m_tag[i] = tags[(rank-1)-i]; }
        for (       ; i < 8    ; ++i ) { m_tag[i] = NULL ; }
      }
      else {
        size_type n = 1 ;
        unsigned i ;
        for ( i = 0 ; i < rank ; ++i ) { m_stride[i] = n *= dims[i] ; }
        for (       ; i < 8    ; ++i ) { m_stride[i] = 0 ; }
        for ( i = 0 ; i < rank ; ++i ) { m_tag[i] = tags[i]; }
        for (       ; i < 8    ; ++i ) { m_tag[i] = NULL ; }
      }
    }

  /** \} */
protected:

  value_type * m_ptr ;
  unsigned     m_rank ;
  size_type    m_stride[8];
  tag_type     m_tag[8] ;

  template< typename , ArrayOrder ,
            class , class , class , class ,
            class , class , class , class >
  friend class phdmesh::Array ;
};

//----------------------------------------------------------------------
//----------------------------------------------------------------------

} // namespace phdmesh

#undef ARRAY_CHECK

#endif

