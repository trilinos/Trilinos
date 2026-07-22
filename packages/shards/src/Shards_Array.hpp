// @HEADER
// *****************************************************************************
//                Shards : Shared Discretization Tools
//
// Copyright 2008-2011 NTESS and the Shards contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef Shards_Array_hpp
#define Shards_Array_hpp

//----------------------------------------------------------------------

#include <vector>
#include <string>
#include <Shards_SimpleArrayOps.hpp>

//----------------------------------------------------------------------
// Macro to compile in array bounds checking:

#ifdef SHARDS_ARRAY_BOUNDS_CHECKING
#define SHARDS_ARRAY_CHECK( X ) X
#else
#define SHARDS_ARRAY_CHECK( X )
#endif

//----------------------------------------------------------------------

namespace shards {

/** \addtogroup  shards_package_array
 *  \{
 */

namespace array_traits {
typedef int int_t ;
} // namespace array_traits

//----------------------------------------------------------------------
/** \brief  Define <b> Natural </b> (C-language) or
 *          <b> Fortran </b> ordering of array dimensions.
 *          A RankZero array does not have an ordering.
 *
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

//----------------------------------------------------------------------

template< typename Scalar , ArrayOrder Order , 
          class Tag1 = void , class Tag2 = void ,
          class Tag3 = void , class Tag4 = void ,
          class Tag5 = void , class Tag6 = void ,
          class Tag7 = void , class Tag8 = void >
class Array ;

//----------------------------------------------------------------------

/** \class  ArrayDimTag
 *  \brief  Abstract base class for array dimension tags supplied to
 *          the Array template class.
 *  \sa Array
 *
 *  A derived array dimension tag class must provide the
 *  <b> name </b> method and <b> tag </b> singleton method
 *  as in the following example.
 *  <PRE>
 *  struct MyTag : public shards::ArrayDimTag {
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
class ArrayDimTag {
public:

  typedef array_traits::int_t size_type ;

  /** \brief Name of the tag, typically the name of the derived class. */
  virtual const char * name() const = 0 ;

  /** \brief  Given a dimension and index produce a string for output.
   *
   *          Default to converting <b> index </b> to a string.
   */
  virtual std::string to_string( size_type dimension ,
                                 size_type index ) const ;

  /** \brief Given a dimension and input strige produce an index.
   *
   *          Default to converting <b> label </b> to an integer.
   */
  virtual size_type to_index( size_type dimension ,
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
 */
class ArrayDimension : public ArrayDimTag {
public:

  const char * name() const ;

  /** \brief Singleton tag for ArrayDimension */
  static const ArrayDimension & tag();

private:
  ~ArrayDimension();
  ArrayDimension();
  ArrayDimension( const ArrayDimension & );
  ArrayDimension & operator = ( const ArrayDimension & );
};

/** \brief  Macro for declaration of a simple ArrayDimTag
 *  \param ADT  name of the tag.
 */
#define SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION( ADT ) \
  class ADT : public shards::ArrayDimTag { \
  public: \
    const char * name() const ; \
    static const ADT & tag(); \
  private: \
    ~ADT(); \
    ADT(); \
    ADT( const ADT & ); \
    ADT & operator = ( const ADT & ); \
  };

/** \brief  Macro for implementing the body of a simple ArrayDimTag
 *  \param ADT  name of the tag.
 */
#define SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION( ADT ) \
  ADT::ADT() {} \
  ADT::~ADT() {} \
  const char * ADT::name() const { static const char n[] = # ADT; return n; } \
  const ADT & ADT::tag() { static const ADT self ; return self ; }

//----------------------------------------------------------------------
//----------------------------------------------------------------------

/** \} */

} // namespace shards

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// Private implementation details for the array

#ifndef DOXYGEN_COMPILE

namespace shards {
namespace array_traits {

//----------------------------------------------------------------------
/** \brief  Return the total number of members from the array stride */
template< typename iType >
inline
iType stride_size(
  const iType & rank ,
  const iType * const stride )
{ return 0 < rank ? stride[ rank - 1 ] : 0 ; }

/** \brief  Generate natural dimension from array stride */
template< typename iType >
inline
void stride_to_natural_dimensions(
  const iType   rank ,
  const iType * const stride ,
        iType * const dim )
{
  iType n = 1 ;
  for ( iType i = 0 ; i < rank ; ++i )
    { dim[(rank-1)-i] = stride[i] / n ; n = stride[i] ; }
}

/** \brief  Generate natural indices from array stride */
template< typename iType >
inline
void stride_to_natural_indices(
  const iType   rank ,
  const iType * const stride ,
        iType   offset ,
        iType * const indices )
{
  iType * i = indices ;
  for ( const iType * s = stride + rank - 1 ; stride < s-- ; ++i ) {
    *i = offset / *s ;
    offset %= *s ;
  }
  *i = offset ;
}

/** \brief  Generate array stride from natural dimensions */
template< typename iType >
inline
void stride_from_natural_dimensions( 
  const iType rank ,
        iType * const stride ,
  const iType * const dim )
{
  iType n = 1 ;
  for ( iType i = 0 ; i < rank ; ++i ) { stride[i] = n *= dim[(rank-1)-i]; }
}

//----------------------------------------------------------------------

void throw_bad_conversion( const int_t lhs_rank ,
                           const ArrayDimTag * const lhs_tags[] ,
                           const int_t rhs_rank ,
                           const ArrayDimTag * const rhs_tags[] );

void check_rank( const int_t rank ,
                 const int_t test_rank );

void check_range( const int_t index , const int_t bound );

void check_indices( const bool ,
                    const int_t rank ,
                    const int_t * const stride ,
                    const int_t = 0 ,
                    const int_t = 0 ,
                    const int_t = 0 ,
                    const int_t = 0 ,
                    const int_t = 0 ,
                    const int_t = 0 ,
                    const int_t = 0 ,
                    const int_t = 0 );

void init_dim(
         int_t dst_stride[] ,     
   const int_t src_dimension[] ,  
   const int_t rank , const bool natural );
 
void init_tags(
   const ArrayDimTag *       dst_tag[] , 
   const ArrayDimTag * const src_tag[] , 
   const int_t rank , const bool natural );

//----------------------------------------------------------------------

template< int_t , int_t > struct CheckRank ;

template<> struct CheckRank<0,0> { static void ok(){} };
template<> struct CheckRank<1,1> { static void ok(){} };
template<> struct CheckRank<2,2> { static void ok(){} };
template<> struct CheckRank<3,3> { static void ok(){} };
template<> struct CheckRank<4,4> { static void ok(){} };
template<> struct CheckRank<5,5> { static void ok(){} };
template<> struct CheckRank<6,6> { static void ok(){} };
template<> struct CheckRank<7,7> { static void ok(){} };
template<> struct CheckRank<8,8> { static void ok(){} };

//----------------------------------------------------------------------

template< int_t Index , int_t Bound > struct CheckRange ;

template<> struct CheckRange<0,8> { static void ok(){} };
template<> struct CheckRange<1,8> { static void ok(){} };
template<> struct CheckRange<2,8> { static void ok(){} };
template<> struct CheckRange<3,8> { static void ok(){} };
template<> struct CheckRange<4,8> { static void ok(){} };
template<> struct CheckRange<5,8> { static void ok(){} };
template<> struct CheckRange<6,8> { static void ok(){} };
template<> struct CheckRange<7,8> { static void ok(){} };

template<> struct CheckRange<0,7> { static void ok(){} };
template<> struct CheckRange<1,7> { static void ok(){} };
template<> struct CheckRange<2,7> { static void ok(){} };
template<> struct CheckRange<3,7> { static void ok(){} };
template<> struct CheckRange<4,7> { static void ok(){} };
template<> struct CheckRange<5,7> { static void ok(){} };
template<> struct CheckRange<6,7> { static void ok(){} };

template<> struct CheckRange<0,6> { static void ok(){} };
template<> struct CheckRange<1,6> { static void ok(){} };
template<> struct CheckRange<2,6> { static void ok(){} };
template<> struct CheckRange<3,6> { static void ok(){} };
template<> struct CheckRange<4,6> { static void ok(){} };
template<> struct CheckRange<5,6> { static void ok(){} };

template<> struct CheckRange<0,5> { static void ok(){} };
template<> struct CheckRange<1,5> { static void ok(){} };
template<> struct CheckRange<2,5> { static void ok(){} };
template<> struct CheckRange<3,5> { static void ok(){} };
template<> struct CheckRange<4,5> { static void ok(){} };

template<> struct CheckRange<0,4> { static void ok(){} };
template<> struct CheckRange<1,4> { static void ok(){} };
template<> struct CheckRange<2,4> { static void ok(){} };
template<> struct CheckRange<3,4> { static void ok(){} };

template<> struct CheckRange<0,3> { static void ok(){} };
template<> struct CheckRange<1,3> { static void ok(){} };
template<> struct CheckRange<2,3> { static void ok(){} };

template<> struct CheckRange<0,2> { static void ok(){} };
template<> struct CheckRange<1,2> { static void ok(){} };

template<> struct CheckRange<0,1> { static void ok(){} };

//----------------------------------------------------------------------

template< class , int_t > struct TagAt ;

template< typename Scalar , ArrayOrder order ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
struct TagAt<Array<Scalar,order,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>,0>
{ typedef Tag1 type ; };

template< typename Scalar , ArrayOrder order ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
struct TagAt<Array<Scalar,order,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>,1>
{ typedef Tag2 type ; };

template< typename Scalar , ArrayOrder order ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
struct TagAt<Array<Scalar,order,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>,2>
{ typedef Tag3 type ; };

template< typename Scalar , ArrayOrder order ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
struct TagAt<Array<Scalar,order,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>,3>
{ typedef Tag4 type ; };

template< typename Scalar , ArrayOrder order ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
struct TagAt<Array<Scalar,order,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>,4>
{ typedef Tag5 type ; };

template< typename Scalar , ArrayOrder order ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
struct TagAt<Array<Scalar,order,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>,5>
{ typedef Tag6 type ; };

template< typename Scalar , ArrayOrder order ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
struct TagAt<Array<Scalar,order,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>,6>
{ typedef Tag7 type ; };

template< typename Scalar , ArrayOrder order ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
struct TagAt<Array<Scalar,order,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>,7>
{ typedef Tag8 type ; };

//----------------------------------------------------------------------
//----------------------------------------------------------------------

template< ArrayOrder , int_t Rank , int_t Ordinal = 0 > struct StrideDim ;

template< int_t Rank , int_t Ordinal >
struct StrideDim<RankZero,Rank,Ordinal> {

  template< typename iType >
  static iType dimension( const iType * )
    { return 0 ; }

  template< typename iType >
  static iType dimension( const iType * , iType )
    { return 0 ; }
};

template< int_t Rank >
struct StrideDim<FortranOrder,Rank,0> {

  template< typename iType >
  static iType dimension( const iType * stride )
    {
      array_traits::CheckRange<0,Rank>::ok();
      return stride[0];
    }

  template< typename iType >
  static iType dimension( const iType * stride , iType ordinal )
    {
      array_traits::check_range(ordinal,Rank);
      return ordinal ? stride[ordinal] / stride[ordinal-1] : stride[0] ;
    }
};

template< int_t Rank >
struct StrideDim<NaturalOrder,Rank,0> {

  template< typename iType >
  static iType dimension( const iType * stride )
    {
      array_traits::CheckRange<0,Rank>::ok();
      return stride[0];
    }

  template< typename iType >
  static iType dimension( const iType * stride , iType ordinal )
    {
      array_traits::check_range(ordinal,Rank);
      ordinal = ( Rank - 1 ) - ordinal ;
      return ordinal ? stride[ordinal] / stride[ordinal-1] : stride[0] ;
    }
};

template< int_t Rank , int_t Ordinal >
struct StrideDim<FortranOrder,Rank,Ordinal> {

  template< typename iType >
  static iType dimension( const iType * stride )
    {
      array_traits::CheckRange<Ordinal,Rank>::ok();
      return stride[Ordinal] / stride[Ordinal-1];
    }
};

template< int_t Rank , int_t Ordinal >
struct StrideDim<NaturalOrder,Rank,Ordinal> {

  template< typename iType >
  static iType dimension( const iType * stride )
    {
      enum { I = ( Rank - 1 ) - Ordinal };
      array_traits::CheckRange<Ordinal,Rank>::ok();
      return stride[I] / stride[I-1];
    }
};

//----------------------------------------------------------------------

template< ArrayOrder > struct Offset ;

template<>
struct Offset<FortranOrder> {

  template< typename isType , typename iType >
  static iType op( const isType * const stride ,
                   const iType & i1 , const iType & i2 ,
                   const iType & i3 , const iType & i4 ,
                   const iType & i5 , const iType & i6 ,
                   const iType & i7 , const iType & i8 )
  {
    SHARDS_ARRAY_CHECK(check_indices(false,8,stride,i1,i2,i3,i4,i5,i6,i7,i8));
    return i1             + i2 * stride[0] +
           i3 * stride[1] + i4 * stride[2] +
           i5 * stride[3] + i6 * stride[4] +
           i7 * stride[5] + i8 * stride[6] ;
  }

  template< typename isType , typename iType >
  static iType op( const isType * const stride ,
                   const iType & i1 , const iType & i2 ,
                   const iType & i3 , const iType & i4 ,
                   const iType & i5 , const iType & i6 ,
                   const iType & i7 )
  {
    SHARDS_ARRAY_CHECK(check_indices(false,7,stride,i1,i2,i3,i4,i5,i6,i7));
    return i1             + i2 * stride[0] +
           i3 * stride[1] + i4 * stride[2] +
           i5 * stride[3] + i6 * stride[4] +
           i7 * stride[5] ;
  }

  template< typename isType , typename iType >
  static iType op( const isType * const stride ,
                   const iType & i1 , const iType & i2 ,
                   const iType & i3 , const iType & i4 ,
                   const iType & i5 , const iType & i6 )
  {
    SHARDS_ARRAY_CHECK(check_indices(false,6,stride,i1,i2,i3,i4,i5,i6));
    return i1             + i2 * stride[0] +
           i3 * stride[1] + i4 * stride[2] +
           i5 * stride[3] + i6 * stride[4] ;
  }

  template< typename isType , typename iType >
  static iType op( const isType * const stride ,
                   const iType & i1 , const iType & i2 ,
                   const iType & i3 , const iType & i4 ,
                   const iType & i5 )
  {
    SHARDS_ARRAY_CHECK(check_indices(false,5,stride,i1,i2,i3,i4,i5));
    return i1             + i2 * stride[0] +
           i3 * stride[1] + i4 * stride[2] +
           i5 * stride[3] ;
  }

  template< typename isType , typename iType >
  static iType op( const isType * const stride ,
                   const iType & i1 , const iType & i2 ,
                   const iType & i3 , const iType & i4 )
  {
    SHARDS_ARRAY_CHECK(check_indices(false,4,stride,i1,i2,i3,i4));
    return i1             + i2 * stride[0] +
           i3 * stride[1] + i4 * stride[2] ;
  }

  template< typename isType , typename iType >
  static iType op( const isType * const stride ,
                   const iType & i1 , const iType & i2 ,
                   const iType & i3 )
  {
    SHARDS_ARRAY_CHECK(check_indices(false,3,stride,i1,i2,i3));
    return i1 + i2 * stride[0] + i3 * stride[1] ;
  }

  template< typename isType , typename iType >
  static iType op( const isType * const stride ,
                   const iType & i1 , const iType & i2 )
  {
    SHARDS_ARRAY_CHECK(check_indices(false,2,stride,i1,i2));
    return i1 + i2 * stride[0] ;
  }

  template< typename isType , typename iType >
  static iType op( const isType * const SHARDS_ARRAY_CHECK( stride ) ,
                   const iType & i1 )
  {
    SHARDS_ARRAY_CHECK(check_indices(false,1,stride,i1));
    return i1 ;
  }
};

//----------------------------------------------------------------------

template<>
struct Offset<NaturalOrder> {

  template< typename isType , typename iType >
  static iType op( const isType * const stride ,
                   const iType & i1 , const iType & i2 ,
                   const iType & i3 , const iType & i4 ,
                   const iType & i5 , const iType & i6 ,
                   const iType & i7 , const iType & i8 )
  {
    SHARDS_ARRAY_CHECK(check_indices(true,8,stride,i1,i2,i3,i4,i5,i6,i7,i8));
    return i8             + i7 * stride[0] +
           i6 * stride[1] + i5 * stride[2] +
           i4 * stride[3] + i3 * stride[4] +
           i2 * stride[5] + i1 * stride[6] ;
  }

  template< typename isType , typename iType >
  static iType op( const isType * const stride ,
                   const iType & i1 , const iType & i2 ,
                   const iType & i3 , const iType & i4 ,
                   const iType & i5 , const iType & i6 ,
                   const iType & i7 )
  {
    SHARDS_ARRAY_CHECK(check_indices(true,7,stride,i1,i2,i3,i4,i5,i6,i7));
    return i7             + i6 * stride[0] +
           i5 * stride[1] + i4 * stride[2] +
           i3 * stride[3] + i2 * stride[4] +
           i1 * stride[5] ;
  }

  template< typename isType , typename iType >
  static iType op( const isType * const stride ,
                   const iType & i1 , const iType & i2 ,
                   const iType & i3 , const iType & i4 ,
                   const iType & i5 , const iType & i6 )
  {
    SHARDS_ARRAY_CHECK(check_indices(true,6,stride,i1,i2,i3,i4,i5,i6));
    return i6             + i5 * stride[0] +
           i4 * stride[1] + i3 * stride[2] +
           i2 * stride[3] + i1 * stride[4] ;
  }

  template< typename isType , typename iType >
  static iType op( const isType * const stride ,
                   const iType & i1 , const iType & i2 ,
                   const iType & i3 , const iType & i4 ,
                   const iType & i5 )
  {
    SHARDS_ARRAY_CHECK(check_indices(true,5,stride,i1,i2,i3,i4,i5));
    return i5             + i4 * stride[0] +
           i3 * stride[1] + i2 * stride[2] +
           i1 * stride[3] ;
  }

  template< typename isType , typename iType >
  static iType op( const isType * const stride ,
                   const iType & i1 , const iType & i2 ,
                   const iType & i3 , const iType & i4 )
  {
    SHARDS_ARRAY_CHECK(check_indices(true,4,stride,i1,i2,i3,i4));
    return i4             + i3 * stride[0] +
           i2 * stride[1] + i1 * stride[2] ;
  }

  template< typename isType , typename iType >
  static iType op( const isType * const stride ,
                   const iType & i1 , const iType & i2 ,
                   const iType & i3 )
  {
    SHARDS_ARRAY_CHECK(check_indices(true,3,stride,i1,i2,i3));
    return i3 + i2 * stride[0] + i1 * stride[1] ;
  }

  template< typename isType , typename iType >
  static iType op( const isType * const stride ,
                   const iType & i1 , const iType & i2 )
  {
    SHARDS_ARRAY_CHECK(check_indices(true,2,stride,i1,i2));
    return i2 + i1 * stride[0] ;
  }

  template< typename isType , typename iType >
  static iType op( const isType * const SHARDS_ARRAY_CHECK( stride ) ,
                   const iType & i1 )
  {
    SHARDS_ARRAY_CHECK(check_indices(true,1,stride,i1));
    return i1 ;
  }
};

//----------------------------------------------------------------------
//----------------------------------------------------------------------

template< typename Scalar , ArrayOrder Order ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
struct Helper ;

//----------------------------------------------------------------------
/** \brief  Rank 8 array types */

template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
struct Helper<Scalar,NaturalOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>
{
  typedef
    Array<Scalar,FortranOrder,Tag8,Tag7,Tag6,Tag5,Tag4,Tag3,Tag2,Tag1>
      reverse ;

  typedef
    Array<Scalar,NaturalOrder,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8,void>
      truncate ;

  enum { Rank = 8 };

  static bool verify( const int_t rank , const ArrayDimTag * const tags[] )
    {
      return rank    == Rank &&
             tags[0] == & Tag8::tag() &&
             tags[1] == & Tag7::tag() &&
             tags[2] == & Tag6::tag() &&
             tags[3] == & Tag5::tag() &&
             tags[4] == & Tag4::tag() &&
             tags[5] == & Tag3::tag() &&
             tags[6] == & Tag2::tag() &&
             tags[7] == & Tag1::tag() ;
    }

  static void assign_tags( const ArrayDimTag * tags[] )
    {
      tags[0] = & Tag8::tag();
      tags[1] = & Tag7::tag();
      tags[2] = & Tag6::tag();
      tags[3] = & Tag5::tag();
      tags[4] = & Tag4::tag();
      tags[5] = & Tag3::tag();
      tags[6] = & Tag2::tag();
      tags[7] = & Tag1::tag();
    }

  template< typename iType >
  static void assign( iType * stride )
    {
        stride[7] = Tag1::Size * (
        stride[6] = Tag2::Size * (
        stride[5] = Tag3::Size * (
        stride[4] = Tag4::Size * (
        stride[3] = Tag5::Size * (
        stride[2] = Tag6::Size * (
        stride[1] = Tag7::Size * (
        stride[0] = Tag8::Size )))))));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n1 )
    {
        stride[7] = n1 * (
        stride[6] = Tag2::Size * (
        stride[5] = Tag3::Size * (
        stride[4] = Tag4::Size * (
        stride[3] = Tag5::Size * (
        stride[2] = Tag6::Size * (
        stride[1] = Tag7::Size * (
        stride[0] = Tag8::Size )))))));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n1 ,
                      const iType & n2 )
    {
        stride[7] = n1 * (
        stride[6] = n2 * (
        stride[5] = Tag3::Size * (
        stride[4] = Tag4::Size * (
        stride[3] = Tag5::Size * (
        stride[2] = Tag6::Size * (
        stride[1] = Tag7::Size * (
        stride[0] = Tag8::Size )))))));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n1 ,
                      const iType & n2 ,
                      const iType & n3 )
    {
        stride[7] = n1 * (
        stride[6] = n2 * (
        stride[5] = n3 * (
        stride[4] = Tag4::Size * (
        stride[3] = Tag5::Size * (
        stride[2] = Tag6::Size * (
        stride[1] = Tag7::Size * (
        stride[0] = Tag8::Size )))))));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n1 ,
                      const iType & n2 ,
                      const iType & n3 ,
                      const iType & n4 )
    {
        stride[7] = n1 * (
        stride[6] = n2 * (
        stride[5] = n3 * (
        stride[4] = n4 * (
        stride[3] = Tag5::Size * (
        stride[2] = Tag6::Size * (
        stride[1] = Tag7::Size * (
        stride[0] = Tag8::Size )))))));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n1 ,
                      const iType & n2 ,
                      const iType & n3 ,
                      const iType & n4 ,
                      const iType & n5 )
    {
        stride[7] = n1 * (
        stride[6] = n2 * (
        stride[5] = n3 * (
        stride[4] = n4 * (
        stride[3] = n5 * (
        stride[2] = Tag6::Size  * (
        stride[1] = Tag7::Size  * (
        stride[0] = Tag8::Size  )))))));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n1 ,
                      const iType & n2 ,
                      const iType & n3 ,
                      const iType & n4 ,
                      const iType & n5 ,
                      const iType & n6 )
    {
        stride[7] = n1 * (
        stride[6] = n2 * (
        stride[5] = n3 * (
        stride[4] = n4 * (
        stride[3] = n5 * (
        stride[2] = n6 * (
        stride[1] = Tag7::Size * (
        stride[0] = Tag8::Size )))))));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n1 ,
                      const iType & n2 ,
                      const iType & n3 ,
                      const iType & n4 ,
                      const iType & n5 ,
                      const iType & n6 ,
                      const iType & n7 )
    {
        stride[7] = n1 * (
        stride[6] = n2 * (
        stride[5] = n3 * (
        stride[4] = n4 * (
        stride[3] = n5 * (
        stride[2] = n6 * (
        stride[1] = n7 * (
        stride[0] = Tag8::Size )))))));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n1 ,
                      const iType & n2 ,
                      const iType & n3 ,
                      const iType & n4 ,
                      const iType & n5 ,
                      const iType & n6 ,
                      const iType & n7 ,
                      const iType & n8 )
    {
        stride[7] = n1 * (
        stride[6] = n2 * (
        stride[5] = n3 * (
        stride[4] = n4 * (
        stride[3] = n5 * (
        stride[2] = n6 * (
        stride[1] = n7 * (
        stride[0] = n8 )))))));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType * const dims )
    {
        stride[7] = dims[0] * (
        stride[6] = dims[1] * (
        stride[5] = dims[2] * (
        stride[4] = dims[3] * (
        stride[3] = dims[4] * (
        stride[2] = dims[5] * (
        stride[1] = dims[6] * (
        stride[0] = dims[7] )))))));
    }
};

template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
struct Helper<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>
{
  typedef
    Array<Scalar,NaturalOrder,Tag8,Tag7,Tag6,Tag5,Tag4,Tag3,Tag2,Tag1>
      reverse ;

  typedef
    Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,void>
      truncate ;

  enum { Rank = 8 };

  static bool verify( int_t rank , const ArrayDimTag * tags[] )
    {
      return rank == Rank &&
             tags[0] == & Tag1::tag() &&
             tags[1] == & Tag2::tag() &&
             tags[2] == & Tag3::tag() &&
             tags[3] == & Tag4::tag() &&
             tags[4] == & Tag5::tag() &&
             tags[5] == & Tag6::tag() &&
             tags[6] == & Tag7::tag() &&
             tags[7] == & Tag8::tag();
    }

  static void assign_tags( const ArrayDimTag * tags[] )
    {
      tags[0] = & Tag1::tag();
      tags[1] = & Tag2::tag();
      tags[2] = & Tag3::tag();
      tags[3] = & Tag4::tag();
      tags[4] = & Tag5::tag();
      tags[5] = & Tag6::tag();
      tags[6] = & Tag7::tag();
      tags[7] = & Tag8::tag();
    }

  template< typename iType >
  static void assign( iType * stride )
    {
        stride[7] = Tag8::Size * (
        stride[6] = Tag7::Size * (
        stride[5] = Tag6::Size * (
        stride[4] = Tag5::Size * (
        stride[3] = Tag4::Size * (
        stride[2] = Tag3::Size * (
        stride[1] = Tag2::Size * (
        stride[0] = Tag1::Size )))))));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n8 )
    {
        stride[7] = n8 * (
        stride[6] = Tag7::Size * (
        stride[5] = Tag6::Size * (
        stride[4] = Tag5::Size * (
        stride[3] = Tag4::Size * (
        stride[2] = Tag3::Size * (
        stride[1] = Tag2::Size * (
        stride[0] = Tag1::Size )))))));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n7 ,
                      const iType & n8 )
    {
        stride[7] = n8 * (
        stride[6] = n7 * (
        stride[5] = Tag6::Size * (
        stride[4] = Tag5::Size * (
        stride[3] = Tag4::Size * (
        stride[2] = Tag3::Size * (
        stride[1] = Tag2::Size * (
        stride[0] = Tag1::Size )))))));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n6 ,
                      const iType & n7 ,
                      const iType & n8 )
    {
        stride[7] = n8 * (
        stride[6] = n7 * (
        stride[5] = n6 * (
        stride[4] = Tag5::Size * (
        stride[3] = Tag4::Size * (
        stride[2] = Tag3::Size * (
        stride[1] = Tag2::Size * (
        stride[0] = Tag1::Size )))))));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n5 ,
                      const iType & n6 ,
                      const iType & n7 ,
                      const iType & n8 )
    {
        stride[7] = n8 * (
        stride[6] = n7 * (
        stride[5] = n6 * (
        stride[4] = n5 * (
        stride[3] = Tag4::Size * (
        stride[2] = Tag3::Size * (
        stride[1] = Tag2::Size * (
        stride[0] = Tag1::Size )))))));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n4 ,
                      const iType & n5 ,
                      const iType & n6 ,
                      const iType & n7 ,
                      const iType & n8 )
    {
        stride[7] = n8 * (
        stride[6] = n7 * (
        stride[5] = n6 * (
        stride[4] = n5 * (
        stride[3] = n4 * (
        stride[2] = Tag3::Size  * (
        stride[1] = Tag2::Size  * (
        stride[0] = Tag1::Size  )))))));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n3 ,
                      const iType & n4 ,
                      const iType & n5 ,
                      const iType & n6 ,
                      const iType & n7 ,
                      const iType & n8 )
    {
        stride[7] = n8 * (
        stride[6] = n7 * (
        stride[5] = n6 * (
        stride[4] = n5 * (
        stride[3] = n4 * (
        stride[2] = n3 * (
        stride[1] = Tag2::Size * (
        stride[0] = Tag1::Size )))))));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n2 ,
                      const iType & n3 ,
                      const iType & n4 ,
                      const iType & n5 ,
                      const iType & n6 ,
                      const iType & n7 ,
                      const iType & n8 )
    {
        stride[7] = n8 * (
        stride[6] = n7 * (
        stride[5] = n6 * (
        stride[4] = n5 * (
        stride[3] = n4 * (
        stride[2] = n3 * (
        stride[1] = n2 * (
        stride[0] = Tag1::Size )))))));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n1 ,
                      const iType & n2 ,
                      const iType & n3 ,
                      const iType & n4 ,
                      const iType & n5 ,
                      const iType & n6 ,
                      const iType & n7 ,
                      const iType & n8 )
    {
        stride[7] = n8 * (
        stride[6] = n7 * (
        stride[5] = n6 * (
        stride[4] = n5 * (
        stride[3] = n4 * (
        stride[2] = n3 * (
        stride[1] = n2 * (
        stride[0] = n1 )))))));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType * const dims )
    {
        stride[7] = dims[7] * (
        stride[6] = dims[6] * (
        stride[5] = dims[5] * (
        stride[4] = dims[4] * (
        stride[3] = dims[3] * (
        stride[2] = dims[2] * (
        stride[1] = dims[1] * (
        stride[0] = dims[0] )))))));
    }
};

//----------------------------------------------------------------------
/** \brief  Rank 7 array types */

template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 >
struct Helper<Scalar,NaturalOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,void>
{
  typedef
    Array<Scalar,FortranOrder,Tag7,Tag6,Tag5,Tag4,Tag3,Tag2,Tag1,void>
      reverse ;

  typedef
    Array<Scalar,NaturalOrder,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,void,void>
      truncate ;

  template< class TagA >
  struct append {
    typedef
      Array<Scalar,NaturalOrder,TagA,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>
        natural ;

    typedef
      Array<Scalar,FortranOrder,Tag7,Tag6,Tag5,Tag4,Tag3,Tag2,Tag1,TagA>
        fortran ;

    typedef  natural  type ;
    typedef  fortran  reverse ;
  };

  enum { Rank = 7 };

  static bool verify( int_t rank , const ArrayDimTag * tags[] )
    {
      return rank == Rank &&
             tags[0] == & Tag7::tag() &&
             tags[1] == & Tag6::tag() &&
             tags[2] == & Tag5::tag() &&
             tags[3] == & Tag4::tag() &&
             tags[4] == & Tag3::tag() &&
             tags[5] == & Tag2::tag() &&
             tags[6] == & Tag1::tag() ;
    }

  static void assign_tags( const ArrayDimTag * tags[] )
    {
      tags[0] = & Tag7::tag();
      tags[1] = & Tag6::tag();
      tags[2] = & Tag5::tag();
      tags[3] = & Tag4::tag();
      tags[4] = & Tag3::tag();
      tags[5] = & Tag2::tag();
      tags[6] = & Tag1::tag();
      tags[7] = NULL ;
    }

  template< typename iType >
  static void assign( iType * stride )
    {
      stride[7] = 0 ;
      stride[6] = Tag1::Size * (
      stride[5] = Tag2::Size * (
      stride[4] = Tag3::Size * (
      stride[3] = Tag4::Size * (
      stride[2] = Tag5::Size * (
      stride[1] = Tag6::Size * (
      stride[0] = Tag7::Size ))))));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n1 )
    {
      stride[7] = 0 ;
      stride[6] = n1 * (
      stride[5] = Tag2::Size * (
      stride[4] = Tag3::Size * (
      stride[3] = Tag4::Size * (
      stride[2] = Tag5::Size * (
      stride[1] = Tag6::Size * (
      stride[0] = Tag7::Size ))))));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n1 ,
                      const iType & n2 )
    {
      stride[7] = 0 ;
      stride[6] = n1 * (
      stride[5] = n2 * (
      stride[4] = Tag3::Size * (
      stride[3] = Tag4::Size * (
      stride[2] = Tag5::Size * (
      stride[1] = Tag6::Size * (
      stride[0] = Tag7::Size ))))));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n1 ,
                      const iType & n2 ,
                      const iType & n3 )
    {
      stride[7] = 0 ;
      stride[6] = n1 * (
      stride[5] = n2 * (
      stride[4] = n3 * (
      stride[3] = Tag4::Size * (
      stride[2] = Tag5::Size * (
      stride[1] = Tag6::Size * (
      stride[0] = Tag7::Size ))))));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n1 ,
                      const iType & n2 ,
                      const iType & n3 ,
                      const iType & n4 )
    {
      stride[7] = 0 ;
      stride[6] = n1 * (
      stride[5] = n2 * (
      stride[4] = n3 * (
      stride[3] = n4 * (
      stride[2] = Tag5::Size * (
      stride[1] = Tag6::Size * (
      stride[0] = Tag7::Size ))))));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n1 ,
                      const iType & n2 ,
                      const iType & n3 ,
                      const iType & n4 ,
                      const iType & n5 )
    {
      stride[7] = 0 ;
      stride[6] = n1 * (
      stride[5] = n2 * (
      stride[4] = n3 * (
      stride[3] = n4 * (
      stride[2] = n5 * (
      stride[1] = Tag6::Size * (
      stride[0] = Tag7::Size ))))));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n1 ,
                      const iType & n2 ,
                      const iType & n3 ,
                      const iType & n4 ,
                      const iType & n5 ,
                      const iType & n6 )
    {
      stride[7] = 0 ;
      stride[6] = n1 * (
      stride[5] = n2 * (
      stride[4] = n3 * (
      stride[3] = n4 * (
      stride[2] = n5 * (
      stride[1] = n6 * (
      stride[0] = Tag7::Size ))))));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n1 ,
                      const iType & n2 ,
                      const iType & n3 ,
                      const iType & n4 ,
                      const iType & n5 ,
                      const iType & n6 ,
                      const iType & n7 )
    {
      stride[7] = 0 ;
      stride[6] = n1 * (
      stride[5] = n2 * (
      stride[4] = n3 * (
      stride[3] = n4 * (
      stride[2] = n5 * (
      stride[1] = n6 * (
      stride[0] = n7 ))))));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType * const dims )
    {
      stride[7] = 0 ;
      stride[6] = dims[0] * (
      stride[5] = dims[1] * (
      stride[4] = dims[2] * (
      stride[3] = dims[3] * (
      stride[2] = dims[4] * (
      stride[1] = dims[5] * (
      stride[0] = dims[6] ))))));
    }
};

template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 >
struct Helper<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,void>
{
  typedef
    Array<Scalar,NaturalOrder,Tag7,Tag6,Tag5,Tag4,Tag3,Tag2,Tag1,void>
      reverse ;

  typedef
    Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,void,void>
      truncate ;

  template< class TagA >
  struct append {
    typedef
      Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,TagA>
        fortran ;

    typedef
      Array<Scalar,NaturalOrder,TagA,Tag7,Tag6,Tag5,Tag4,Tag3,Tag2,Tag1>
        natural ;

    typedef  fortran  type ;
    typedef  natural  reverse ;
  };

  enum { Rank = 7 };

  static bool verify( int_t rank , const ArrayDimTag * tags[] )
    {
      return rank == Rank &&
             tags[0] == & Tag1::tag() &&
             tags[1] == & Tag2::tag() &&
             tags[2] == & Tag3::tag() &&
             tags[3] == & Tag4::tag() &&
             tags[4] == & Tag5::tag() &&
             tags[5] == & Tag6::tag() &&
             tags[6] == & Tag7::tag();
    }

  static void assign_tags( const ArrayDimTag * tags[] )
    {
      tags[0] = & Tag1::tag();
      tags[1] = & Tag2::tag();
      tags[2] = & Tag3::tag();
      tags[3] = & Tag4::tag();
      tags[4] = & Tag5::tag();
      tags[5] = & Tag6::tag();
      tags[6] = & Tag7::tag();
      tags[7] = NULL ;
    }

  template< typename iType >
  static void assign( iType * stride )
    {
      stride[7] = 0 ;
      stride[6] = Tag7::Size * (
      stride[5] = Tag6::Size * (
      stride[4] = Tag5::Size * (
      stride[3] = Tag4::Size * (
      stride[2] = Tag3::Size * (
      stride[1] = Tag2::Size * (
      stride[0] = Tag1::Size ))))));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n7 )
    {
      stride[7] = 0 ;
      stride[6] = n7 * (
      stride[5] = Tag6::Size * (
      stride[4] = Tag5::Size * (
      stride[3] = Tag4::Size * (
      stride[2] = Tag3::Size * (
      stride[1] = Tag2::Size * (
      stride[0] = Tag1::Size ))))));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n6 ,
                      const iType & n7 )
    {
      stride[7] = 0 ;
      stride[6] = n7 * (
      stride[5] = n6 * (
      stride[4] = Tag5::Size * (
      stride[3] = Tag4::Size * (
      stride[2] = Tag3::Size * (
      stride[1] = Tag2::Size * (
      stride[0] = Tag1::Size ))))));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n5 ,
                      const iType & n6 ,
                      const iType & n7 )
    {
      stride[7] = 0 ;
      stride[6] = n7 * (
      stride[5] = n6 * (
      stride[4] = n5 * (
      stride[3] = Tag4::Size * (
      stride[2] = Tag3::Size * (
      stride[1] = Tag2::Size * (
      stride[0] = Tag1::Size ))))));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n4 ,
                      const iType & n5 ,
                      const iType & n6 ,
                      const iType & n7 )
    {
      stride[7] = 0 ;
      stride[6] = n7 * (
      stride[5] = n6 * (
      stride[4] = n5 * (
      stride[3] = n4 * (
      stride[2] = Tag3::Size * (
      stride[1] = Tag2::Size * (
      stride[0] = Tag1::Size ))))));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n3 ,
                      const iType & n4 ,
                      const iType & n5 ,
                      const iType & n6 ,
                      const iType & n7 )
    {
      stride[7] = 0 ;
      stride[6] = n7 * (
      stride[5] = n6 * (
      stride[4] = n5 * (
      stride[3] = n4 * (
      stride[2] = n3 * (
      stride[1] = Tag2::Size * (
      stride[0] = Tag1::Size ))))));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n2 ,
                      const iType & n3 ,
                      const iType & n4 ,
                      const iType & n5 ,
                      const iType & n6 ,
                      const iType & n7 )
    {
      stride[7] = 0 ;
      stride[6] = n7 * (
      stride[5] = n6 * (
      stride[4] = n5 * (
      stride[3] = n4 * (
      stride[2] = n3 * (
      stride[1] = n2 * (
      stride[0] = Tag1::Size ))))));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n1 ,
                      const iType & n2 ,
                      const iType & n3 ,
                      const iType & n4 ,
                      const iType & n5 ,
                      const iType & n6 ,
                      const iType & n7 )
    {
      stride[7] = 0 ;
      stride[6] = n7 * (
      stride[5] = n6 * (
      stride[4] = n5 * (
      stride[3] = n4 * (
      stride[2] = n3 * (
      stride[1] = n2 * (
      stride[0] = n1 ))))));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType * const dims )
    {
      stride[7] = 0 ;
      stride[6] = dims[6] * (
      stride[5] = dims[5] * (
      stride[4] = dims[4] * (
      stride[3] = dims[3] * (
      stride[2] = dims[2] * (
      stride[1] = dims[1] * (
      stride[0] = dims[0] ))))));
   }
};

//----------------------------------------------------------------------
/** \brief  Rank 6 array types */

template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 >
struct Helper<Scalar,NaturalOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,void,void>
{
  typedef
    Array<Scalar,FortranOrder,Tag6,Tag5,Tag4,Tag3,Tag2,Tag1,void,void>
      reverse ;

  typedef
    Array<Scalar,NaturalOrder,Tag2,Tag3,Tag4,Tag5,Tag6,void,void,void>
      truncate ;

  template< class TagA >
  struct append {
    typedef
      Array<Scalar,NaturalOrder,TagA,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,void>
        natural ;

    typedef
      Array<Scalar,FortranOrder,Tag6,Tag5,Tag4,Tag3,Tag2,Tag1,TagA,void>
        fortran ;

    typedef  natural  type ;
    typedef  fortran  reverse ;
  };

  enum { Rank = 6 };

  static bool verify( int_t rank , const ArrayDimTag * tags[] )
    {
      return rank == Rank &&
             tags[0] == & Tag6::tag() &&
             tags[1] == & Tag5::tag() &&
             tags[2] == & Tag4::tag() &&
             tags[3] == & Tag3::tag() &&
             tags[4] == & Tag2::tag() &&
             tags[5] == & Tag1::tag();
    }

  static void assign_tags( const ArrayDimTag * tags[] )
    {
      tags[0] = & Tag6::tag();
      tags[1] = & Tag5::tag();
      tags[2] = & Tag4::tag();
      tags[3] = & Tag3::tag();
      tags[4] = & Tag2::tag();
      tags[5] = & Tag1::tag();
      tags[6] = NULL ;
      tags[7] = NULL ;
    }

  template< typename iType >
  static void assign( iType * stride )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = Tag1::Size * (
      stride[4] = Tag2::Size * (
      stride[3] = Tag3::Size * (
      stride[2] = Tag4::Size * (
      stride[1] = Tag5::Size * (
      stride[0] = Tag6::Size )))));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n1 )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = n1 * (
      stride[4] = Tag2::Size * (
      stride[3] = Tag3::Size * (
      stride[2] = Tag4::Size * (
      stride[1] = Tag5::Size * (
      stride[0] = Tag6::Size )))));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n1 ,
                      const iType & n2 )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = n1 * (
      stride[4] = n2 * (
      stride[3] = Tag3::Size * (
      stride[2] = Tag4::Size * (
      stride[1] = Tag5::Size * (
      stride[0] = Tag6::Size )))));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n1 ,
                      const iType & n2 ,
                      const iType & n3 )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = n1 * (
      stride[4] = n2 * (
      stride[3] = n3 * (
      stride[2] = Tag4::Size * (
      stride[1] = Tag5::Size * (
      stride[0] = Tag6::Size )))));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n1 ,
                      const iType & n2 ,
                      const iType & n3 ,
                      const iType & n4 )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = n1 * (
      stride[4] = n2 * (
      stride[3] = n3 * (
      stride[2] = n4 * (
      stride[1] = Tag5::Size * (
      stride[0] = Tag6::Size )))));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n1 ,
                      const iType & n2 ,
                      const iType & n3 ,
                      const iType & n4 ,
                      const iType & n5 )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = n1 * (
      stride[4] = n2 * (
      stride[3] = n3 * (
      stride[2] = n4 * (
      stride[1] = n5 * (
      stride[0] = Tag6::Size )))));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n1 ,
                      const iType & n2 ,
                      const iType & n3 ,
                      const iType & n4 ,
                      const iType & n5 ,
                      const iType & n6 )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = n1 * (
      stride[4] = n2 * (
      stride[3] = n3 * (
      stride[2] = n4 * (
      stride[1] = n5 * (
      stride[0] = n6 )))));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType * const dims )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = dims[0] * (
      stride[4] = dims[1] * (
      stride[3] = dims[2] * (
      stride[2] = dims[3] * (
      stride[1] = dims[4] * (
      stride[0] = dims[5] )))));
    }
};

template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 >
struct Helper<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,void,void>
{
  typedef
    Array<Scalar,NaturalOrder,Tag6,Tag5,Tag4,Tag3,Tag2,Tag1,void,void>
      reverse ;

  typedef
    Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,Tag5,void,void,void>
      truncate ;

  template< class TagA >
  struct append {
    typedef
      Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,TagA,void>
        fortran ;

    typedef
      Array<Scalar,NaturalOrder,TagA,Tag6,Tag5,Tag4,Tag3,Tag2,Tag1,void>
        natural ;

    typedef  fortran  type ;
    typedef  natural  reverse ;
  };

  enum { Rank = 6 };

  static bool verify( int_t rank , const ArrayDimTag * tags[] )
    {
      return rank == Rank &&
             tags[0] == & Tag1::tag() &&
             tags[1] == & Tag2::tag() &&
             tags[2] == & Tag3::tag() &&
             tags[3] == & Tag4::tag() &&
             tags[4] == & Tag5::tag() &&
             tags[5] == & Tag6::tag();
    }

  static void assign_tags( const ArrayDimTag * tags[] )
    {
      tags[0] = & Tag1::tag();
      tags[1] = & Tag2::tag();
      tags[2] = & Tag3::tag();
      tags[3] = & Tag4::tag();
      tags[4] = & Tag5::tag();
      tags[5] = & Tag6::tag();
      tags[6] = NULL ;
      tags[7] = NULL ;
    }

  template< typename iType >
  static void assign( iType * stride )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = Tag6::Size * (
      stride[4] = Tag5::Size * (
      stride[3] = Tag4::Size * (
      stride[2] = Tag3::Size * (
      stride[1] = Tag2::Size * (
      stride[0] = Tag1::Size )))));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n6 )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = n6 * (
      stride[4] = Tag5::Size * (
      stride[3] = Tag4::Size * (
      stride[2] = Tag3::Size * (
      stride[1] = Tag2::Size * (
      stride[0] = Tag1::Size )))));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n5 ,
                      const iType & n6 )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = n6 * (
      stride[4] = n5 * (
      stride[3] = Tag4::Size * (
      stride[2] = Tag3::Size * (
      stride[1] = Tag2::Size * (
      stride[0] = Tag1::Size )))));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n4 ,
                      const iType & n5 ,
                      const iType & n6 )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = n6 * (
      stride[4] = n5 * (
      stride[3] = n4 * (
      stride[2] = Tag3::Size * (
      stride[1] = Tag2::Size * (
      stride[0] = Tag1::Size )))));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n3 ,
                      const iType & n4 ,
                      const iType & n5 ,
                      const iType & n6 )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = n6 * (
      stride[4] = n5 * (
      stride[3] = n4 * (
      stride[2] = n3 * (
      stride[1] = Tag2::Size * (
      stride[0] = Tag1::Size )))));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n2 ,
                      const iType & n3 ,
                      const iType & n4 ,
                      const iType & n5 ,
                      const iType & n6 )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = n6 * (
      stride[4] = n5 * (
      stride[3] = n4 * (
      stride[2] = n3 * (
      stride[1] = n2 * (
      stride[0] = Tag1::Size )))));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n1 ,
                      const iType & n2 ,
                      const iType & n3 ,
                      const iType & n4 ,
                      const iType & n5 ,
                      const iType & n6 )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = n6 * (
      stride[4] = n5 * (
      stride[3] = n4 * (
      stride[2] = n3 * (
      stride[1] = n2 * (
      stride[0] = n1 )))));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType * const dims )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = dims[5] * (
      stride[4] = dims[4] * (
      stride[3] = dims[3] * (
      stride[2] = dims[2] * (
      stride[1] = dims[1] * (
      stride[0] = dims[0] )))));
   }
};

//----------------------------------------------------------------------
/** \brief  Rank 5 array types */

template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 , class Tag5 >
struct Helper<Scalar,NaturalOrder,Tag1,Tag2,Tag3,Tag4,Tag5,void,void,void>
{
  typedef
    Array<Scalar,FortranOrder,Tag5,Tag4,Tag3,Tag2,Tag1,void,void,void>
      reverse ;

  typedef
    Array<Scalar,NaturalOrder,Tag2,Tag3,Tag4,Tag5,void,void,void,void>
      truncate ;

  template< class TagA >
  struct append {
    typedef
      Array<Scalar,NaturalOrder,TagA,Tag1,Tag2,Tag3,Tag4,Tag5,void,void>
        natural ;

    typedef
      Array<Scalar,FortranOrder,Tag5,Tag4,Tag3,Tag2,Tag1,TagA,void,void>
        fortran ;

    typedef  natural  type ;
    typedef  fortran  reverse ;
  };

  enum { Rank = 5 };

  static bool verify( int_t rank , const ArrayDimTag * tags[] )
    {
      return rank == Rank &&
             tags[0] == & Tag5::tag() &&
             tags[1] == & Tag4::tag() &&
             tags[2] == & Tag3::tag() &&
             tags[3] == & Tag2::tag() &&
             tags[4] == & Tag1::tag();
    }

  static void assign_tags( const ArrayDimTag * tags[] )
    {
      tags[0] = & Tag5::tag();
      tags[1] = & Tag4::tag();
      tags[2] = & Tag3::tag();
      tags[3] = & Tag2::tag();
      tags[4] = & Tag1::tag();
      tags[5] = NULL ;
      tags[6] = NULL ;
      tags[7] = NULL ;
    }

  template< typename iType >
  static void assign( iType * stride )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = 0 ;
      stride[4] = Tag1::Size * (
      stride[3] = Tag2::Size * (
      stride[2] = Tag3::Size * (
      stride[1] = Tag4::Size * (
      stride[0] = Tag5::Size ))));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n1 )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = 0 ;
      stride[4] = n1 * (
      stride[3] = Tag2::Size * (
      stride[2] = Tag3::Size * (
      stride[1] = Tag4::Size * (
      stride[0] = Tag5::Size ))));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n1 ,
                      const iType & n2 )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = 0 ;
      stride[4] = n1 * (
      stride[3] = n2 * (
      stride[2] = Tag3::Size * (
      stride[1] = Tag4::Size * (
      stride[0] = Tag5::Size ))));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n1 ,
                      const iType & n2 ,
                      const iType & n3 )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = 0 ;
      stride[4] = n1 * (
      stride[3] = n2 * (
      stride[2] = n3 * (
      stride[1] = Tag4::Size * (
      stride[0] = Tag5::Size ))));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n1 ,
                      const iType & n2 ,
                      const iType & n3 ,
                      const iType & n4 )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = 0 ;
      stride[4] = n1 * (
      stride[3] = n2 * (
      stride[2] = n3 * (
      stride[1] = n4 * (
      stride[0] = Tag5::Size ))));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n1 ,
                      const iType & n2 ,
                      const iType & n3 ,
                      const iType & n4 ,
                      const iType & n5 )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = 0 ;
      stride[4] = n1 * (
      stride[3] = n2 * (
      stride[2] = n3 * (
      stride[1] = n4 * (
      stride[0] = n5 ))));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType * const dims )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = 0 ;
      stride[4] = dims[0] * (
      stride[3] = dims[1] * (
      stride[2] = dims[2] * (
      stride[1] = dims[3] * (
      stride[0] = dims[4] ))));
    }
};

template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 , class Tag5 >
struct Helper<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,Tag5,void,void,void>
{
  typedef
    Array<Scalar,NaturalOrder,Tag5,Tag4,Tag3,Tag2,Tag1,void,void,void>
      reverse ;

  typedef
    Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,void,void,void,void>
      truncate ;

  template< class TagA >
  struct append {
    typedef
      Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,Tag5,TagA,void,void>
        fortran ;

    typedef
      Array<Scalar,NaturalOrder,TagA,Tag5,Tag4,Tag3,Tag2,Tag1,void,void>
        natural ;

    typedef  fortran  type ;
    typedef  natural  reverse ;
  };

  enum { Rank = 5 };

  static bool verify( int_t rank , const ArrayDimTag * tags[] )
    {
      return rank == Rank &&
             tags[0] == & Tag1::tag() &&
             tags[1] == & Tag2::tag() &&
             tags[2] == & Tag3::tag() &&
             tags[3] == & Tag4::tag() &&
             tags[4] == & Tag5::tag();
    }

  static void assign_tags( const ArrayDimTag * tags[] )
    {
      tags[0] = & Tag1::tag();
      tags[1] = & Tag2::tag();
      tags[2] = & Tag3::tag();
      tags[3] = & Tag4::tag();
      tags[4] = & Tag5::tag();
      tags[5] = NULL ;
      tags[6] = NULL ;
      tags[7] = NULL ;
    }

  template< typename iType >
  static void assign( iType * stride )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = 0 ;
      stride[4] = Tag5::Size * (
      stride[3] = Tag4::Size * (
      stride[2] = Tag3::Size * (
      stride[1] = Tag2::Size * (
      stride[0] = Tag1::Size ))));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n5 )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = 0 ;
      stride[4] = n5 * (
      stride[3] = Tag4::Size * (
      stride[2] = Tag3::Size * (
      stride[1] = Tag2::Size * (
      stride[0] = Tag1::Size ))));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n4 ,
                      const iType & n5 )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = 0 ;
      stride[4] = n5 * (
      stride[3] = n4 * (
      stride[2] = Tag3::Size * (
      stride[1] = Tag2::Size * (
      stride[0] = Tag1::Size ))));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n3 ,
                      const iType & n4 ,
                      const iType & n5 )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = 0 ;
      stride[4] = n5 * (
      stride[3] = n4 * (
      stride[2] = n3 * (
      stride[1] = Tag2::Size * (
      stride[0] = Tag1::Size ))));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n2 ,
                      const iType & n3 ,
                      const iType & n4 ,
                      const iType & n5 )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = 0 ;
      stride[4] = n5 * (
      stride[3] = n4 * (
      stride[2] = n3 * (
      stride[1] = n2 * (
      stride[0] = Tag1::Size ))));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n1 ,
                      const iType & n2 ,
                      const iType & n3 ,
                      const iType & n4 ,
                      const iType & n5 )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = 0 ;
      stride[4] = n5 * (
      stride[3] = n4 * (
      stride[2] = n3 * (
      stride[1] = n2 * (
      stride[0] = n1 ))));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType * const dims )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = 0 ;
      stride[4] = dims[4] * (
      stride[3] = dims[3] * (
      stride[2] = dims[2] * (
      stride[1] = dims[1] * (
      stride[0] = dims[0] ))));
   }
};

//----------------------------------------------------------------------
/** \brief  Rank 4 array types */

template< typename Scalar , class Tag1 , class Tag2 , class Tag3 , class Tag4 >
struct Helper<Scalar,NaturalOrder,Tag1,Tag2,Tag3,Tag4,void,void,void,void>
{
  typedef
    Array<Scalar,FortranOrder,Tag4,Tag3,Tag2,Tag1,void,void,void,void>
      reverse ;

  typedef
    Array<Scalar,NaturalOrder,Tag2,Tag3,Tag4,void,void,void,void,void>
      truncate ;

  template< class TagA >
  struct append {
    typedef
      Array<Scalar,NaturalOrder,TagA,Tag1,Tag2,Tag3,Tag4,void,void,void>
        natural ;

    typedef
      Array<Scalar,FortranOrder,Tag4,Tag3,Tag2,Tag1,TagA,void,void,void>
        fortran ;

    typedef  natural  type ;
    typedef  fortran  reverse ;
  };

  enum { Rank = 4 };

  static bool verify( int_t rank , const ArrayDimTag * tags[] )
    {
      return rank == Rank &&
             tags[0] == & Tag4::tag() &&
             tags[1] == & Tag3::tag() &&
             tags[2] == & Tag2::tag() &&
             tags[3] == & Tag1::tag();
    }

  static void assign_tags( const ArrayDimTag * tags[] )
    {
      tags[0] = & Tag4::tag();
      tags[1] = & Tag3::tag();
      tags[2] = & Tag2::tag();
      tags[3] = & Tag1::tag();
      tags[4] = NULL ;
      tags[5] = NULL ;
      tags[6] = NULL ;
      tags[7] = NULL ;
    }

  template< typename iType >
  static void assign( iType * stride )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = 0 ;
      stride[4] = 0 ;
      stride[3] = Tag1::Size * (
      stride[2] = Tag2::Size * (
      stride[1] = Tag3::Size * (
      stride[0] = Tag4::Size )));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n1 )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = 0 ;
      stride[4] = 0 ;
      stride[3] = n1 * (
      stride[2] = Tag2::Size * (
      stride[1] = Tag3::Size * (
      stride[0] = Tag4::Size )));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n1 ,
                      const iType & n2 )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = 0 ;
      stride[4] = 0 ;
      stride[3] = n1 * (
      stride[2] = n2 * (
      stride[1] = Tag3::Size * (
      stride[0] = Tag4::Size )));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n1 ,
                      const iType & n2 ,
                      const iType & n3 )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = 0 ;
      stride[4] = 0 ;
      stride[3] = n1 * (
      stride[2] = n2 * (
      stride[1] = n3 * (
      stride[0] = Tag4::Size )));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n1 ,
                      const iType & n2 ,
                      const iType & n3 ,
                      const iType & n4 )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = 0 ;
      stride[4] = 0 ;
      stride[3] = n1 * (
      stride[2] = n2 * (
      stride[1] = n3 * (
      stride[0] = n4 )));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType * const dims )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = 0 ;
      stride[4] = 0 ;
      stride[3] = dims[0] * (
      stride[2] = dims[1] * (
      stride[1] = dims[2] * (
      stride[0] = dims[3] )));
    }
};

template< typename Scalar , class Tag1 , class Tag2 , class Tag3 , class Tag4 >
struct Helper<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,void,void,void,void>
{
  typedef
    Array<Scalar,NaturalOrder,Tag4,Tag3,Tag2,Tag1,void,void,void,void>
      reverse ;

  typedef
    Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,void,void,void,void,void>
      truncate ;

  template< class TagA >
  struct append {
    typedef
      Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,TagA,void,void,void>
        fortran ;

    typedef
      Array<Scalar,NaturalOrder,TagA,Tag4,Tag3,Tag2,Tag1,void,void,void>
        natural ;

    typedef  fortran  type ;
    typedef  natural  reverse ;
  };

  enum { Rank = 4 };

  static bool verify( int_t rank , const ArrayDimTag * tags[] )
    {
      return rank == Rank &&
             tags[0] == & Tag1::tag() &&
             tags[1] == & Tag2::tag() &&
             tags[2] == & Tag3::tag() &&
             tags[3] == & Tag4::tag();
    }

  static void assign_tags( const ArrayDimTag * tags[] )
    {
      tags[0] = & Tag1::tag();
      tags[1] = & Tag2::tag();
      tags[2] = & Tag3::tag();
      tags[3] = & Tag4::tag();
      tags[4] = NULL ;
      tags[5] = NULL ;
      tags[6] = NULL ;
      tags[7] = NULL ;
    }

  template< typename iType >
  static void assign( iType * stride )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = 0 ;
      stride[4] = 0 ;
      stride[3] = Tag4::Size * (
      stride[2] = Tag3::Size * (
      stride[1] = Tag2::Size * (
      stride[0] = Tag1::Size )));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n4 )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = 0 ;
      stride[4] = 0 ;
      stride[3] = n4 * (
      stride[2] = Tag3::Size * (
      stride[1] = Tag2::Size * (
      stride[0] = Tag1::Size )));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n3 ,
                      const iType & n4 )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = 0 ;
      stride[4] = 0 ;
      stride[3] = n4 * (
      stride[2] = n3 * (
      stride[1] = Tag2::Size * (
      stride[0] = Tag1::Size )));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n2 ,
                      const iType & n3 ,
                      const iType & n4 )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = 0 ;
      stride[4] = 0 ;
      stride[3] = n4 * (
      stride[2] = n3 * (
      stride[1] = n2 * (
      stride[0] = Tag1::Size )));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n1 ,
                      const iType & n2 ,
                      const iType & n3 ,
                      const iType & n4 )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = 0 ;
      stride[4] = 0 ;
      stride[3] = n4 * (
      stride[2] = n3 * (
      stride[1] = n2 * (
      stride[0] = n1 )));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType * const dims )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = 0 ;
      stride[4] = 0 ;
      stride[3] = dims[3] * (
      stride[2] = dims[2] * (
      stride[1] = dims[1] * (
      stride[0] = dims[0] )));
   }
};

//----------------------------------------------------------------------
/** \brief  Rank 3 array types */

template< typename Scalar , class Tag1 , class Tag2 , class Tag3 >
struct Helper<Scalar,NaturalOrder,Tag1,Tag2,Tag3,void,void,void,void,void>
{
  typedef
    Array<Scalar,FortranOrder,Tag3,Tag2,Tag1,void,void,void,void,void>
      reverse ;

  typedef
    Array<Scalar,NaturalOrder,Tag2,Tag3,void,void,void,void,void,void>
      truncate ;

  template< class TagA >
  struct append {
    typedef
      Array<Scalar,NaturalOrder,TagA,Tag1,Tag2,Tag3,void,void,void,void>
        natural ;

    typedef
      Array<Scalar,FortranOrder,Tag3,Tag2,Tag1,TagA,void,void,void,void>
        fortran ;

    typedef  natural  type ;
    typedef  fortran  reverse ;
  };

  enum { Rank = 3 };

  static bool verify( int_t rank , const ArrayDimTag * tags[] )
    {
      return rank == Rank &&
             tags[0] == & Tag3::tag() &&
             tags[1] == & Tag2::tag() &&
             tags[2] == & Tag1::tag();
    }

  static void assign_tags( const ArrayDimTag * tags[] )
    {
      tags[0] = & Tag3::tag();
      tags[1] = & Tag2::tag();
      tags[2] = & Tag1::tag();
      tags[3] = NULL ;
      tags[4] = NULL ;
      tags[5] = NULL ;
      tags[6] = NULL ;
      tags[7] = NULL ;
    }

  template< typename iType >
  static void assign( iType * stride )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = 0 ;
      stride[4] = 0 ;
      stride[3] = 0 ;
      stride[2] = Tag1::Size * (
      stride[1] = Tag2::Size * (
      stride[0] = Tag3::Size ));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n1 )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = 0 ;
      stride[4] = 0 ;
      stride[3] = 0 ;
      stride[2] = n1 * (
      stride[1] = Tag2::Size * (
      stride[0] = Tag3::Size ));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n1 ,
                      const iType & n2 )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = 0 ;
      stride[4] = 0 ;
      stride[3] = 0 ;
      stride[2] = n1 * (
      stride[1] = n2 * (
      stride[0] = Tag3::Size ));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n1 ,
                      const iType & n2 ,
                      const iType & n3 )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = 0 ;
      stride[4] = 0 ;
      stride[3] = 0 ;
      stride[2] = n1 * (
      stride[1] = n2 * (
      stride[0] = n3 ));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType * const dims )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = 0 ;
      stride[4] = 0 ;
      stride[3] = 0 ;
      stride[2] = dims[0] * (
      stride[1] = dims[1] * (
      stride[0] = dims[2] ));
    }
};

template< typename Scalar , class Tag1 , class Tag2 , class Tag3 >
struct Helper<Scalar,FortranOrder,Tag1,Tag2,Tag3,void,void,void,void,void>
{
  typedef
    Array<Scalar,NaturalOrder,Tag3,Tag2,Tag1,void,void,void,void,void>
      reverse ;

  typedef
    Array<Scalar,FortranOrder,Tag1,Tag2,void,void,void,void,void,void>
      truncate ;

  template< class TagA >
  struct append {
    typedef
      Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,TagA,void,void,void,void>
        fortran ;

    typedef
      Array<Scalar,NaturalOrder,TagA,Tag3,Tag2,Tag1,void,void,void,void>
        natural ;

    typedef  fortran  type ;
    typedef  natural  reverse ;
  };

  enum { Rank = 3 };

  static bool verify( int_t rank , const ArrayDimTag * tags[] )
    {
      return rank == Rank &&
             tags[0] == & Tag1::tag() &&
             tags[1] == & Tag2::tag() &&
             tags[2] == & Tag3::tag();
    }

  static void assign_tags( const ArrayDimTag * tags[] )
    {
      tags[0] = & Tag1::tag();
      tags[1] = & Tag2::tag();
      tags[2] = & Tag3::tag();
      tags[3] = NULL ;
      tags[4] = NULL ;
      tags[5] = NULL ;
      tags[6] = NULL ;
      tags[7] = NULL ;
    }

  template< typename iType >
  static void assign( iType * stride )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = 0 ;
      stride[4] = 0 ;
      stride[3] = 0 ;
      stride[2] = Tag3::Size * (
      stride[1] = Tag2::Size * (
      stride[0] = Tag1::Size ));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n3 )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = 0 ;
      stride[4] = 0 ;
      stride[3] = 0 ;
      stride[2] = n3 * (
      stride[1] = Tag2::Size * (
      stride[0] = Tag1::Size ));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n2 ,
                      const iType & n3 )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = 0 ;
      stride[4] = 0 ;
      stride[3] = 0 ;
      stride[2] = n3 * (
      stride[1] = n2 * (
      stride[0] = Tag1::Size ));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n1 ,
                      const iType & n2 ,
                      const iType & n3 )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = 0 ;
      stride[4] = 0 ;
      stride[3] = 0 ;
      stride[2] = n3 * (
      stride[1] = n2 * (
      stride[0] = n1 ));
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType * const dims )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = 0 ;
      stride[4] = 0 ;
      stride[3] = 0 ;
      stride[2] = dims[2] * (
      stride[1] = dims[1] * (
      stride[0] = dims[0] ));
   }
};

//----------------------------------------------------------------------
/** \brief  Rank 2 array types */

template< typename Scalar , class Tag1 , class Tag2 >
struct Helper<Scalar,NaturalOrder,Tag1,Tag2,void,void,void,void,void,void>
{
  typedef
    Array<Scalar,FortranOrder,Tag2,Tag1,void,void,void,void,void,void>
      reverse ;

  typedef
    Array<Scalar,NaturalOrder,Tag2,void,void,void,void,void,void,void>
      truncate ;

  template< class TagA >
  struct append {
    typedef
      Array<Scalar,NaturalOrder,TagA,Tag1,Tag2,void,void,void,void,void>
        natural ;

    typedef
      Array<Scalar,FortranOrder,Tag2,Tag1,TagA,void,void,void,void,void>
        fortran ;

    typedef  natural  type ;
    typedef  fortran  reverse ;
  };

  enum { Rank = 2 };

  static bool verify( int_t rank , const ArrayDimTag * tags[] )
    {
      return rank == Rank &&
             tags[0] == & Tag2::tag() &&
             tags[1] == & Tag1::tag();
    }

  static void assign_tags( const ArrayDimTag * tags[] )
    {
      tags[0] = & Tag2::tag();
      tags[1] = & Tag1::tag();
      tags[2] = NULL ;
      tags[3] = NULL ;
      tags[4] = NULL ;
      tags[5] = NULL ;
      tags[6] = NULL ;
      tags[7] = NULL ;
    }

  template< typename iType >
  static void assign( iType * stride )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = 0 ;
      stride[4] = 0 ;
      stride[3] = 0 ;
      stride[2] = 0 ;
      stride[1] = Tag1::Size * (
      stride[0] = Tag2::Size );
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n1 )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = 0 ;
      stride[4] = 0 ;
      stride[3] = 0 ;
      stride[2] = 0 ;
      stride[1] = n1 * (
      stride[0] = Tag2::Size );
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n1 ,
                      const iType & n2 )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = 0 ;
      stride[4] = 0 ;
      stride[3] = 0 ;
      stride[2] = 0 ;
      stride[1] = n1 * (
      stride[0] = n2 );
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType * const dims )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = 0 ;
      stride[4] = 0 ;
      stride[3] = 0 ;
      stride[2] = 0 ;
      stride[1] = dims[0] * (
      stride[0] = dims[1] );
    }
};

template< typename Scalar , class Tag1 , class Tag2 >
struct Helper<Scalar,FortranOrder,Tag1,Tag2,void,void,void,void,void,void>
{
  typedef
    Array<Scalar,NaturalOrder,Tag2,Tag1,void,void,void,void,void,void>
      reverse ;

  typedef
    Array<Scalar,FortranOrder,Tag1,void,void,void,void,void,void,void>
      truncate ;

  template< class TagA >
  struct append {
    typedef
      Array<Scalar,FortranOrder,Tag1,Tag2,TagA,void,void,void,void,void>
        fortran ;

    typedef
      Array<Scalar,NaturalOrder,TagA,Tag2,Tag1,void,void,void,void,void>
        natural ;

    typedef  fortran  type ;
    typedef  natural  reverse ;
  };

  enum { Rank = 2 };

  static bool verify( int_t rank , const ArrayDimTag * tags[] )
    {
      return rank == Rank &&
             tags[0] = & Tag1::tag() &&
             tags[1] = & Tag2::tag();
    }

  static void assign_tags( const ArrayDimTag * tags[] )
    {
      tags[0] = & Tag1::tag();
      tags[1] = & Tag2::tag();
      tags[2] = NULL ;
      tags[3] = NULL ;
      tags[4] = NULL ;
      tags[5] = NULL ;
      tags[6] = NULL ;
      tags[7] = NULL ;
    }

  template< typename iType >
  static void assign( iType * stride )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = 0 ;
      stride[4] = 0 ;
      stride[3] = 0 ;
      stride[2] = 0 ;
      stride[1] = Tag2::Size * (
      stride[0] = Tag1::Size );
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n2 )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = 0 ;
      stride[4] = 0 ;
      stride[3] = 0 ;
      stride[2] = 0 ;
      stride[1] = n2 * (
      stride[0] = Tag1::Size );
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n1 ,
                      const iType & n2 )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = 0 ;
      stride[4] = 0 ;
      stride[3] = 0 ;
      stride[2] = 0 ;
      stride[1] = n2 * (
      stride[0] = n1 );
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType * const dims )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = 0 ;
      stride[4] = 0 ;
      stride[3] = 0 ;
      stride[2] = 0 ;
      stride[1] = dims[1] * (
      stride[0] = dims[0] );
    }
};

//----------------------------------------------------------------------
/** \brief  Rank 1 array types */

template< typename Scalar , class Tag1 >
struct Helper<Scalar,NaturalOrder,Tag1,void,void,void,void,void,void,void>
{
  typedef
    Array<Scalar,FortranOrder,Tag1,void,void,void,void,void,void,void>
      reverse ;

  typedef
    Array<Scalar,RankZero,void,void,void,void,void,void,void,void>
      truncate ;

  template< class TagA >
  struct append {
    typedef
      Array<Scalar,NaturalOrder,TagA,Tag1,void,void,void,void,void,void>
        natural ;

    typedef
      Array<Scalar,FortranOrder,Tag1,TagA,void,void,void,void,void,void>
        fortran ;

    typedef  natural  type ;
    typedef  fortran  reverse ;
  };

  enum { Rank = 1 };

  static bool verify( int_t rank , const ArrayDimTag * tags[] )
    { return rank == Rank && tags[0] == & Tag1::tag(); }

  static void assign_tags( const ArrayDimTag * tags[] )
    {
      tags[0] = & Tag1::tag();
      tags[1] = NULL ;
      tags[2] = NULL ;
      tags[3] = NULL ;
      tags[4] = NULL ;
      tags[5] = NULL ;
      tags[6] = NULL ;
      tags[7] = NULL ;
    }

  template< typename iType >
  static void assign( iType * stride )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = 0 ;
      stride[4] = 0 ;
      stride[3] = 0 ;
      stride[2] = 0 ;
      stride[1] = 0 ;
      stride[0] = Tag1::Size ;
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType & n1 )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = 0 ;
      stride[4] = 0 ;
      stride[3] = 0 ;
      stride[2] = 0 ;
      stride[1] = 0 ;
      stride[0] = n1 ;
    }

  template< typename iType >
  static void assign( iType * stride ,
                      const iType * const dims )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = 0 ;
      stride[4] = 0 ;
      stride[3] = 0 ;
      stride[2] = 0 ;
      stride[1] = 0 ;
      stride[0] = dims[0] ;
    }
};

template< typename Scalar , class Tag1 >
struct Helper<Scalar,FortranOrder,Tag1,void,void,void,void,void,void,void>
{
  typedef
    Array<Scalar,NaturalOrder,Tag1,void,void,void,void,void,void,void>
      reverse ;

  typedef
    Array<Scalar,RankZero,void,void,void,void,void,void,void,void>
      truncate ;

  template< class TagA >
  struct append {
    typedef
      Array<Scalar,FortranOrder,Tag1,TagA,void,void,void,void,void,void>
        fortran ;

    typedef
      Array<Scalar,NaturalOrder,TagA,Tag1,void,void,void,void,void,void>
        natural ;

    typedef  fortran  type ;
    typedef  natural  reverse ;
  };

  enum { Rank = 1 };

  static bool verify( int_t rank , const ArrayDimTag * tags[] )
    { return rank == Rank && tags[0] == & Tag1::tag(); }

  static void assign_tags( const ArrayDimTag * tags[] )
    {
      tags[0] = & Tag1::tag();
      tags[1] = NULL ;
      tags[2] = NULL ;
      tags[3] = NULL ;
      tags[4] = NULL ;
      tags[5] = NULL ;
      tags[6] = NULL ;
      tags[7] = NULL ;
    }

  template< typename iType >
  static void assign( iType * stride )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = 0 ;
      stride[4] = 0 ;
      stride[3] = 0 ;
      stride[2] = 0 ;
      stride[1] = 0 ;
      stride[0] = Tag1::Size ;
    }

  template< typename iType >
  static void assign( iType * stride , const iType & n1 )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = 0 ;
      stride[4] = 0 ;
      stride[3] = 0 ;
      stride[2] = 0 ;
      stride[1] = 0 ;
      stride[0] = n1 ;
    }

  template< typename iType >
  static void assign( iType * stride , const iType * const dims )
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = 0 ;
      stride[4] = 0 ;
      stride[3] = 0 ;
      stride[2] = 0 ;
      stride[0] = dims[0] ;
    }
};

//----------------------------------------------------------------------
/** \brief  Rank Zero array types */

template< typename Scalar >
struct Helper<Scalar,RankZero,void,void,void,void,void,void,void,void>
{
  typedef
    Array<Scalar,RankZero,void,void,void,void,void,void,void,void>
      reverse ;

  template< class TagA >
  struct append {
    typedef
      Array<Scalar,NaturalOrder,TagA,void,void,void,void,void,void,void>
        natural ;

    typedef
      Array<Scalar,FortranOrder,TagA,void,void,void,void,void,void,void>
        fortran ;

    typedef natural type ;
    typedef fortran reverse ;
  };

  enum { Rank = 0 };

  template< typename iType >
  static void assign( iType * ) {}
};

//----------------------------------------------------------------------
/** \brief  Rank Dynamic array types */

template< typename Scalar >
struct Helper<Scalar,NaturalOrder,void,void,void,void,void,void,void,void>
{
  typedef
    Array<Scalar,FortranOrder,void,void,void,void,void,void,void,void>
      reverse ;
};

template< typename Scalar >
struct Helper<Scalar,FortranOrder,void,void,void,void,void,void,void,void>
{
  typedef
    Array<Scalar,NaturalOrder,void,void,void,void,void,void,void,void>
      reverse ;
};

//----------------------------------------------------------------------

} // namespace array_traits

template< class ArrayType , class TagA > struct ArrayAppend {};

template< typename Scalar , ArrayOrder Order ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class TagA >
struct ArrayAppend<
       Array<Scalar,Order,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,void> , TagA >
{
private:
  typedef typename
    array_traits::Helper<Scalar,Order,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,void>
      ::template append<TagA> helper ;
public:
  typedef typename helper::type    type ;
  typedef typename helper::natural natural ;
  typedef typename helper::fortran fortran ;
};

} // namespace shards

#endif /* DOXYGEN_COMPILE */

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace shards {

/** \addtogroup  shards_package_array
 *  \{
 */

//----------------------------------------------------------------------
/** \brief  The multi-dimensional Array interface with <b> runtime </b>
 *          user-defined dimension ordinates.
 *          Typically used when runtime-polymorphic arrays are passed to
 *          functions.
 *          
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
private:
  typedef typename array_traits::Offset<array_order> Offset ;
public:
  /** \name Array Attributes
   *  \{
   */

  /** \brief  Type of member data. */
  typedef Scalar  value_type ;

  /** \brief  Type for sizes. */
  typedef array_traits::int_t size_type ;

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
  size_type rank()   const { return m_rank ; }

  /** \brief  If the multidimension follows the natural ordering */
  bool natural()    const { return Natural ; }

  /** \brief  If the multidimension follows the reverse (Fortran) ordering */
  bool reverse()    const { return Reverse ; }

  /** \brief  If the member data storage is contiguous */
  bool contiguous() const { return Contiguous ; }

  //----------------------------------

  /** \brief  Access the dimension tag-singleton for a given ordinate. */
  tag_type tag( size_type ord ) const
    {
      array_traits::check_range( ord , m_rank );
      if ( Natural ) { ord = ( m_rank - 1 ) - ord ; }
      return m_tag[ord];
    }

  //----------------------------------

  /** \brief  Dimension of the given ordinate. */
  size_type dimension( size_type ord ) const
    {
      array_traits::check_range( ord , m_rank );
      if ( Natural ) { ord = ( m_rank - 1 ) - ord ; }
      return ord ? m_stride[ord] / m_stride[ord-1] : m_stride[ord] ;
    }

  /** \brief  Dimension of all ordinate. */
  template< typename iType >
  void dimensions( std::vector<iType> & n )
    {
      n.resize( m_rank );
      for ( size_type i = 0 ; i < m_rank ; ++i ) { n[i] = dimension(i); }
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
  template< typename iType >
  Array truncate( const iType & i ) const
    { return Array( *this , i ); }

  /** \brief Pointer to contiguous block of member data. */
  value_type * contiguous_data() const { return m_ptr ; }

  /** \brief Access member via full ordering of members. */
  template< typename iType >
  value_type & operator[]( const iType & i ) const
    {
      SHARDS_ARRAY_CHECK( array_traits::check_range(i,size()) );
      return m_ptr[ i ];
    }

  //----------------------------------
  /** \brief Access member via Rank 8 multi-index */
  template< typename iType >
  value_type & operator()( const iType & i1 , const iType & i2 ,
                           const iType & i3 , const iType & i4 ,
                           const iType & i5 , const iType & i6 ,
                           const iType & i7 , const iType & i8 ) const
    {
      SHARDS_ARRAY_CHECK( array_traits::check_rank( m_rank , 8 ) );
      return m_ptr[ Offset::op(m_stride,i1,i2,i3,i4,i5,i6,i7,i8) ];
    }

  template< typename iType >
  value_type & operator()( const iType & i1 , const iType & i2 ,
                           const iType & i3 , const iType & i4 ,
                           const iType & i5 , const iType & i6 ,
                           const iType & i7 ) const
    {
      SHARDS_ARRAY_CHECK( array_traits::check_rank( m_rank , 7 ) );
      return m_ptr[ Offset::op(m_stride,i1,i2,i3,i4,i5,i6,i7) ];
    }

  template< typename iType >
  value_type & operator()( const iType & i1 , const iType & i2 ,
                           const iType & i3 , const iType & i4 ,
                           const iType & i5 , const iType & i6 ) const
    {
      SHARDS_ARRAY_CHECK( array_traits::check_rank( m_rank , 6 ) );
      return m_ptr[ Offset::op(m_stride,i1,i2,i3,i4,i5,i6) ];
    }

  template< typename iType >
  value_type & operator()( const iType & i1 , const iType & i2 ,
                           const iType & i3 , const iType & i4 ,
                           const iType & i5 ) const
    {
      SHARDS_ARRAY_CHECK( array_traits::check_rank( m_rank , 5 ) );
      return m_ptr[ Offset::op(m_stride,i1,i2,i3,i4,i5) ];
    }

  template< typename iType >
  value_type & operator()( const iType & i1 , const iType & i2 ,
                           const iType & i3 , const iType & i4 ) const
    {
      SHARDS_ARRAY_CHECK( array_traits::check_rank( m_rank , 4 ) );
      return m_ptr[ Offset::op(m_stride,i1,i2,i3,i4) ];
    }

  template< typename iType >
  value_type & operator()( const iType & i1 , const iType & i2 ,
                           const iType & i3 ) const
    {
      SHARDS_ARRAY_CHECK( array_traits::check_rank( m_rank , 3 ) );
      return m_ptr[ Offset::op(m_stride,i1,i2,i3) ];
    }

  template< typename iType >
  value_type & operator()( const iType & i1 , const iType & i2 ) const
    {
      SHARDS_ARRAY_CHECK( array_traits::check_rank( m_rank , 2 ) );
      return m_ptr[ Offset::op(m_stride,i1,i2) ];
    }

  template< typename iType >
  value_type & operator()( const iType & i1 ) const
    {
      SHARDS_ARRAY_CHECK( array_traits::check_rank( m_rank , 1 ) );
      return m_ptr[ Offset::op(m_stride,i1) ];
    }

  /** \} */
  //----------------------------------
  /** \name Constructors and Assignment Operators
   * \{
   */

  typedef typename
    array_traits::Helper<Scalar,array_order,void,void,void,void,void,void,void,void>
      ::reverse  ReverseType ;

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
  // Class specific constructors:

  Array( value_type * ptr ,
         const size_type input_rank ,
         const size_type * const dims ,
         const tag_type  * const tags )
    : m_ptr( ptr ), m_rank( input_rank )
    {
      array_traits::init_dim( m_stride, dims, input_rank, Natural);
      array_traits::init_tags( m_tag,   tags, input_rank, Natural);
    }

  /** \} */

  template< class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
            class Tag5 , class Tag6 , class Tag7 , class Tag8 >
  Array & assign( value_type * ptr ,
                  size_type n1 , size_type n2 , size_type n3 , size_type n4 ,
                  size_type n5 , size_type n6 , size_type n7 , size_type n8 )
  {
    typedef
      array_traits::Helper<Scalar,array_order,
                           Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>
        helper ;
    m_ptr  = ptr ;
    m_rank = helper::Rank ;
    helper::assign_tags( m_tag );
    helper::assign( m_stride, n1, n2, n3, n4, n5, n6, n7, n8 );
    return *this ;
  }

  template< class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
            class Tag5 , class Tag6 , class Tag7 >
  Array & assign( value_type * ptr ,
                  size_type n1 , size_type n2 , size_type n3 , size_type n4 ,
                  size_type n5 , size_type n6 , size_type n7 )
  {
    typedef
      array_traits::Helper<Scalar,array_order,
                           Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,void>
        helper ;
    m_ptr = ptr ;
    m_rank = helper::Rank ;
    helper::assign_tags( m_tag );
    helper::assign( m_stride, n1, n2, n3, n4, n5, n6, n7 ); return *this ;
  }

  template< class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
            class Tag5 , class Tag6 >
  Array & assign( value_type * ptr ,
                  size_type n1 , size_type n2 , size_type n3 , size_type n4 ,
                  size_type n5 , size_type n6 )
  {
    typedef
      array_traits::Helper<Scalar,array_order,
                           Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,void,void>
        helper ;
    m_ptr = ptr ;
    m_rank = helper::Rank ;
    helper::assign_tags( m_tag );
    helper::assign( m_stride, n1, n2, n3, n4, n5, n6 );
    return *this ;
  }

  template< class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
            class Tag5 >
  Array & assign( value_type * ptr ,
                  size_type n1 , size_type n2 , size_type n3 , size_type n4 ,
                  size_type n5 )
  {
    typedef
      array_traits::Helper<Scalar,array_order,
                           Tag1,Tag2,Tag3,Tag4,Tag5,void,void,void>
        helper ;
    m_ptr = ptr ;
    m_rank = helper::Rank ;
    helper::assign_tags( m_tag );
    helper::assign( m_stride, n1, n2, n3, n4, n5 );
    return *this ;
  }

  template< class Tag1 , class Tag2 , class Tag3 , class Tag4 >
  Array & assign( value_type * ptr ,
                  size_type n1 , size_type n2 , size_type n3 , size_type n4 )
  {
    typedef
      array_traits::Helper<Scalar,array_order,
                           Tag1,Tag2,Tag3,Tag4,void,void,void,void>
        helper ;
    m_ptr = ptr ;
    m_rank = helper::Rank ;
    helper::assign_tags( m_tag );
    helper::assign( m_stride, n1, n2, n3, n4 );
    return *this ;
  }

  template< class Tag1 , class Tag2 , class Tag3 >
  Array & assign( value_type * ptr ,
                  size_type n1 , size_type n2 , size_type n3 )
  {
    typedef
      array_traits::Helper<Scalar,array_order,
                           Tag1,Tag2,Tag3,void,void,void,void,void>
        helper ;
    m_ptr = ptr ;
    m_rank = helper::Rank ;
    helper::assign_tags( m_tag );
    helper::assign( m_stride, n1, n2, n3 );
    return *this ;
  }

  template< class Tag1 , class Tag2 >
  Array & assign( value_type * ptr ,
                  size_type n1 , size_type n2 )
  {
    typedef
      array_traits::Helper<Scalar,array_order,
                           Tag1,Tag2,void,void,void,void,void,void>
        helper ;
    m_ptr = ptr ;
    m_rank = helper::Rank ;
    helper::assign_tags( m_tag );
    helper::assign( m_stride, n1, n2 );
    return *this ;
  }

  template< class Tag1 >
  Array & assign( value_type * ptr ,
                  size_type n1 )
  {
    typedef
      array_traits::Helper<Scalar,array_order,
                           Tag1,void,void,void,void,void,void,void>
        helper ;
    m_ptr = ptr ;
    m_rank = helper::Rank ;
    helper::assign_tags( m_tag );
    helper::assign( m_stride, n1 );
    return *this ;
  }

private:

  /** Truncation constructor */
  Array( const Array & rhs , const size_type i )
    : m_ptr( NULL ), m_rank( 0 )
    {
      if ( 1 < rhs.m_rank ) {
        Copy<8>( m_stride , rhs.m_stride );
        Copy<8>( m_tag , rhs.m_tag );
        m_rank = rhs.m_rank - 1 ;
        m_ptr  = rhs.m_ptr + ( m_rank ? m_stride[ m_rank - 1 ] * i : i );
        m_stride[ m_rank ] = 0 ;
        m_tag[ m_rank ] = 0 ;
      }
      else {
        Copy<8>( m_stride , (size_type) 0 );
        Copy<8>( m_tag , (tag_type) NULL );
      }
    }

  /** \brief Pointer to contiguous block of members */
  value_type * m_ptr ;

  /** \brief Rank of the array */
  size_type  m_rank ;

  /** \brief Array of strides, smallest to largest */
  size_type  m_stride[8];

  /** \brief Array of singleton tags, aligned with strides */
  tag_type   m_tag[8] ;

  template< typename , ArrayOrder ,
            class , class , class , class ,
            class , class , class , class >
  friend class shards::Array ;
};

//----------------------------------------------------------------------
/** \class  Array
 *  \brief  The <b> preferred </b> multi-dimensional Array interface
 *          with <b> compile-time </b> user-defined dimension ordinates.
 *  \nosubgrouping
 *
 *  \tparam Scalar  The "plain old data" type of the array's member data.
 *  \tparam array_order An <b> ArrayOrder </b> value that specifies whether to
 *                      use Natural (a.k.a. C-language) or Fortran ordering 
 *                      for the multi-dimensions and multi-indices.
 *  \tparam Tag#  The <b> Tag# </b> template parameters document the
 *                user-defiend purpose of each dimension ordinate.
 *                The <b> Rank </b> of the array (i.e. the number of dimensions)
 *                is the number of user-defined dimension tags, up to eight.
 *                A user-defined dimension <b> Tag# </b> must be derived from
 *                the <b> ArrayDimTag </b> template class.
 *
 *  \sa ArrayDimTag ArrayOrder
 */
template< typename Scalar , ArrayOrder array_order ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
class Array
{
private:

#ifndef DOXYGEN_COMPILE
  typedef
    array_traits::Helper<Scalar,array_order,
                         Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>
      helper ;
#endif /* DOXYGEN_COMPILE */

public:
  /** \name Array Attributes
   *  \{
   */

  /** \brief  Type of member data. */
  typedef Scalar  value_type ;

  /** \brief  Type for sizes. */
  typedef array_traits::int_t size_type ;

  /** \brief  Type of runtime dimension tags. */
  typedef const ArrayDimTag * tag_type ;

  //----------------------------------

  /** \brief  Rank of the array is the number of non-void dimension tags. */
  enum { Rank = helper::Rank };

  /** \brief  If the multidimension follows the natural ordering */
  enum { Natural = NaturalOrder == array_order };

  /** \brief  If the multidimension follows the reverse (Fortran) ordering */
  enum { Reverse = FortranOrder == array_order };

  /** \brief  If the member data storage is contiguous */
  enum { Contiguous = true };

  /** \brief  Rank of the array is the number of non-void dimension tags. */
  size_type rank()   const { return Rank ; }

  /** \brief  If the multidimension follows the natural ordering */
  bool natural()    const { return Natural ; }

  /** \brief  If the multidimension follows the reverse (Fortran) ordering */
  bool reverse()    const { return Reverse ; }

  /** \brief  If the member data storage is contiguous */
  bool contiguous() const { return Contiguous ; }

  //----------------------------------

#ifndef DOXYGEN_COMPILE
  /** \brief  Access the dimension tag-type for a given ordinate. */
  template < size_type ordinate >
  struct Tag {
    typedef typename array_traits::TagAt<Array,ordinate>::type type ;
  };
#endif

  /** \brief  Access the dimension tag-singleton for a given ordinate. */
  tag_type tag( const size_type ordinate ) const
    { return m_array.tag( ordinate ); }

  //----------------------------------
  /** \brief  Dimension of the given ordinate. */
  template < size_type ordinate > size_type dimension() const
    {
      typedef array_traits::StrideDim<array_order,Rank,ordinate> StrideDim ;
      return StrideDim::dimension(m_array.m_stride);
    }

  /** \brief  Dimension of the given ordinate. */
  size_type dimension( const size_type ordinate ) const
    {
      typedef array_traits::StrideDim<array_order,Rank> StrideDim ;
      return StrideDim::dimension(m_array.m_stride,ordinate);
    }

  /** \brief  Dimensions of all ordinates. */
  template< typename iType >
  void dimensions( std::vector<iType> & n )
    { m_array.template dimensions<iType>( n ); }

  /** \brief  Total number of member data items. */
  size_type size() const { return m_array.m_stride[ Rank - 1 ]; }

  /** \} */
  //----------------------------------
  /** \name Member data access operators
   *  \{
   */

  /** \brief  Subarray type that removes the slowest striding dimension
   *          (first natural or last fortran ordinate).
   */
  typedef typename helper::truncate TruncateType ;

  /** \brief  Generate a subarray view of the array with the
   *          slowest striding ordinate offset by <b> i </b>
   *          and removed.
   */
  template< typename iType >
  TruncateType truncate( const iType & i ) const
    { return TruncateType( m_array , i ); }

  //----------------------------------
  /** \brief Pointer to contiguous block of member data. */
  value_type * contiguous_data() const { return m_array.contiguous_data(); }

  /** \brief Access member via offset into contiguous block. */
  template< typename iType >
  value_type & operator[]( const iType & i ) const
    { return m_array[i]; }

  /** \brief Access member of a Rank 8 array */
  template< typename iType >
  value_type & operator()( const iType & i1 , const iType & i2 ,
                           const iType & i3 , const iType & i4 ,
                           const iType & i5 , const iType & i6 ,
                           const iType & i7 , const iType & i8 ) const
    {
      array_traits::CheckRank<8,Rank>::ok();
      return m_array(i1,i2,i3,i4,i5,i6,i7,i8);
    }

  /** \brief Access member of a Rank 7 array */
  template< typename iType >
  value_type & operator()( const iType & i1 , const iType & i2 ,
                           const iType & i3 , const iType & i4 ,
                           const iType & i5 , const iType & i6 ,
                           const iType & i7 ) const
    {
      array_traits::CheckRank<7,Rank>::ok();
      return m_array(i1,i2,i3,i4,i5,i6,i7);
    }

  /** \brief Access member of a Rank 6 array */
  template< typename iType >
  value_type & operator()( const iType & i1 , const iType & i2 ,
                           const iType & i3 , const iType & i4 ,
                           const iType & i5 , const iType & i6 ) const
    {
      array_traits::CheckRank<6,Rank>::ok();
      return m_array(i1,i2,i3,i4,i5,i6);
    }

  /** \brief Access member of a Rank 5 array */
  template< typename iType >
  value_type & operator()( const iType & i1 , const iType & i2 ,
                           const iType & i3 , const iType & i4 ,
                           const iType & i5 ) const
    {
      array_traits::CheckRank<5,Rank>::ok();
      return m_array(i1,i2,i3,i4,i5);
    }

  /** \brief Access member of a Rank 4 array */
  template< typename iType >
  value_type & operator()( const iType & i1 , const iType & i2 ,
                           const iType & i3 , const iType & i4 ) const
    {
      array_traits::CheckRank<4,Rank>::ok();
      return m_array(i1,i2,i3,i4);
    }

  /** \brief Access member of a Rank 3 array */
  template< typename iType >
  value_type & operator()( const iType & i1 , const iType & i2 ,
                           const iType & i3 ) const
    {
      array_traits::CheckRank<3,Rank>::ok();
      return m_array(i1,i2,i3);
    }

  /** \brief Access member of a Rank 2 array */
  template< typename iType >
  value_type & operator()( const iType & i1 , const iType & i2 ) const
    {
      array_traits::CheckRank<2,Rank>::ok();
      return m_array(i1,i2);
    }

  /** \brief Access member of a Rank 1 array */
  template< typename iType >
  value_type & operator()( const iType & i1 ) const
    {
      array_traits::CheckRank<1,Rank>::ok();
      return m_array(i1);
    }

  /** \} */
  //----------------------------------
  /** \name Constructors and Assignment Operators
   * \{
   */

  /** \brief  The compatible multidimensional array with
   *          reversed multi-index ordering and dimension tags.
   */
  typedef typename helper::reverse  ReverseType ;

  /** \brief Default constructor */
  Array() : m_array()
    { m_array.m_rank = Rank ; helper::assign_tags( m_array.m_tag ); }

  /** \brief Copy constructor */
  Array( const Array & rhs ) : m_array( rhs.m_array ) {}

  /** \brief Assignment operator */
  Array & operator = ( const Array & rhs )
    { m_array.operator=(rhs.m_array); return *this ; }

  /** \brief Copy constructor for compatible reverse type. */
  Array( const ReverseType & rhs ) : m_array( rhs.m_array ) {}

  /** \brief Assignment operator for compatible reverse type. */
  Array & operator = ( const ReverseType & rhs )
    { m_array.operator=(rhs.m_array); return *this ; }

  /** \brief Assign pointer and dimensions. */
  Array & assign( value_type * arg_ptr , const size_type * const dims )
    {
      m_array.m_ptr = arg_ptr ;
      array_traits::init_dim( m_array.m_stride , dims , Rank , Natural );
      return *this ;
    }

  /** \brief Construct with array of dimensions. */
  Array( value_type * arg_ptr , const size_type * const dims )
    : m_array()
    {
      m_array.m_rank = Rank ;
      helper::assign_tags( m_array.m_tag );
      assign( arg_ptr , dims );
    }

  /** \brief  Construct a Rank 8 array */
  Array & assign( value_type * arg_ptr ,
                  const size_type n1 , const size_type n2 ,
                  const size_type n3 , const size_type n4 ,
                  const size_type n5 , const size_type n6 ,
                  const size_type n7 , const size_type n8 )
    {
      array_traits::CheckRange<7,Rank>::ok();
      m_array.m_ptr  = arg_ptr ;
      helper::assign( m_array.m_stride, n1, n2, n3, n4, n5, n6, n7, n8 );
      return *this ;
    }

  /** \brief  Construct a Rank 8 array */
  Array( value_type * arg_ptr ,
         const size_type n1 , const size_type n2 ,
         const size_type n3 , const size_type n4 ,
         const size_type n5 , const size_type n6 ,
         const size_type n7 , const size_type n8 )
    : m_array()
    {
      m_array.m_rank = Rank ;
      helper::assign_tags( m_array.m_tag );
      assign( arg_ptr, n1, n2, n3, n4, n5, n6, n7, n8 );
    }

  /** \brief  Construct a Rank 7..8 array; use Tag#::Size for defaults.
   *          The input dimensions are the 7 slowest strides.
   */
  Array& assign( value_type * arg_ptr ,
                 const size_type n1 , const size_type n2 ,
                 const size_type n3 , const size_type n4 ,
                 const size_type n5 , const size_type n6 ,
                 const size_type n7 )
    {
      array_traits::CheckRange<6,Rank>::ok();
      m_array.m_ptr = arg_ptr ;
      helper::assign( m_array.m_stride, n1, n2, n3, n4, n5, n6, n7 );
      return *this ;
    }

  /** \brief  Construct a Rank 7..8 array; use Tag#::Size for defaults.
   *          The input dimensions are the 7 slowest strides.
   */
  Array( value_type * arg_ptr ,
         const size_type n1 , const size_type n2 ,
         const size_type n3 , const size_type n4 ,
         const size_type n5 , const size_type n6 ,
         const size_type n7 )
    : m_array()
    {
      m_array.m_rank = Rank ;
      helper::assign_tags( m_array.m_tag );
      assign( arg_ptr, n1, n2, n3, n4, n5, n6, n7 );
    }

  /** \brief  Construct a Rank 6..8 array; use Tag#::Size for defaults.
   *          The input dimensions are the 6 slowest strides.
   */
  Array & assign( value_type * arg_ptr ,
                  const size_type n1 , const size_type n2 ,
                  const size_type n3 , const size_type n4 ,
                  const size_type n5 , const size_type n6 )
    {
      array_traits::CheckRange<5,Rank>::ok();
      m_array.m_ptr = arg_ptr ;
      helper::assign( m_array.m_stride, n1, n2, n3, n4, n5, n6 );
      return *this ;
    }

  /** \brief  Construct a Rank 6..8 array; use Tag#::Size for defaults.
   *          The input dimensions are the 6 slowest strides.
   */
  Array( value_type * arg_ptr ,
         const size_type n1 , const size_type n2 ,
         const size_type n3 , const size_type n4 ,
         const size_type n5 , const size_type n6 )
    : m_array()
    {
      m_array.m_rank = Rank ;
      helper::assign_tags( m_array.m_tag );
      assign( arg_ptr, n1, n2, n3, n4, n5, n6 );
    }

  /** \brief  Construct a Rank 5..8 array; use Tag#::Size for defaults.
   *          The input dimensions are the 5 slowest strides.
   */
  Array & assign( value_type * arg_ptr ,
                  const size_type n1 , const size_type n2 ,
                  const size_type n3 , const size_type n4 ,
                  const size_type n5 )
    {
      array_traits::CheckRange<4,Rank>::ok();
      m_array.m_ptr  = arg_ptr ;
      helper::assign( m_array.m_stride, n1, n2, n3, n4, n5 );
      return *this ;
    }

  /** \brief  Construct a Rank 5..8 array; use Tag#::Size for defaults.
   *          The input dimensions are the 5 slowest strides.
   */
  Array( value_type * arg_ptr ,
         const size_type n1 , const size_type n2 ,
         const size_type n3 , const size_type n4 ,
         const size_type n5 )
    : m_array()
    {
      m_array.m_rank = Rank ;
      helper::assign_tags( m_array.m_tag );
      assign( arg_ptr, n1, n2, n3, n4, n5 );
    }

  /** \brief  Construct a Rank 4..8 array; use Tag#::Size for defaults.
   *          The input dimensions are the 4 slowest strides.
   */
  Array & assign( value_type * arg_ptr ,
                  const size_type n1 , const size_type n2 ,
                  const size_type n3 , const size_type n4 )
    {
      array_traits::CheckRange<3,Rank>::ok();
      m_array.m_ptr  = arg_ptr ;
      helper::assign( m_array.m_stride, n1, n2, n3, n4 );
      return *this ;
    }

  /** \brief  Construct a Rank 4..8 array; use Tag#::Size for defaults.
   *          The input dimensions are the 4 slowest strides.
   */
  Array( value_type * arg_ptr ,
         const size_type n1 , const size_type n2 ,
         const size_type n3 , const size_type n4 )
    : m_array()
    {
      m_array.m_rank = Rank ;
      helper::assign_tags( m_array.m_tag );
      assign( arg_ptr, n1, n2, n3, n4 );
    }

  /** \brief  Construct a Rank 3..8 array; use Tag#::Size for defaults.
   *          The input dimensions are the 3 slowest strides.
   */
  Array & assign( value_type * arg_ptr ,
                  const size_type n1 , const size_type n2 ,
                  const size_type n3 )
    {
      array_traits::CheckRange<2,Rank>::ok();
      m_array.m_ptr  = arg_ptr ;
      helper::assign( m_array.m_stride, n1, n2, n3 );
      return *this ;
    }

  /** \brief  Construct a Rank 3..8 array; use Tag#::Size for defaults.
   *          The input dimensions are the 3 slowest strides.
   */
  Array( value_type * arg_ptr ,
         const size_type n1 , const size_type n2 ,
         const size_type n3 )
    : m_array()
    {
      m_array.m_rank = Rank ;
      helper::assign_tags( m_array.m_tag );
      assign( arg_ptr , n1, n2, n3 );
    }

  /** \brief  Construct a Rank 2..8 array; use Tag#::Size for defaults.
   *          The input dimensions are the 2 slowest strides.
   */
  Array & assign( value_type * arg_ptr ,
                  const size_type n1 , const size_type n2 )
    {
      array_traits::CheckRange<1,Rank>::ok();
      m_array.m_ptr  = arg_ptr ;
      helper::assign( m_array.m_stride, n1, n2 );
      return *this ;
    }

  /** \brief  Construct a Rank 2..8 array; use Tag#::Size for defaults.
   *          The input dimensions are the 2 slowest strides.
   */
  Array( value_type * arg_ptr , const size_type n1 , const size_type n2 )
    : m_array()
    {
      m_array.m_rank = Rank ;
      helper::assign_tags( m_array.m_tag );
      assign( arg_ptr, n1, n2 );
    }

  /** \brief  Construct a Rank 1..8 array; use Tag#::Size for defaults.
   *          The input dimension is the slowest stride.
   */
  Array & assign( value_type * arg_ptr , const size_type n1 )
    {
      array_traits::CheckRange<0,Rank>::ok();
      m_array.m_ptr  = arg_ptr ;
      helper::assign( m_array.m_stride, n1 );
      return *this ;
    }

  /** \brief  Construct a Rank 1..8 array; use Tag#::Size for defaults.
   *          The input dimension is the slowest stride.
   */
  Array( value_type * arg_ptr , const size_type n1 )
    : m_array()
    {
      m_array.m_rank = Rank ;
      helper::assign_tags( m_array.m_tag );
      assign( arg_ptr, n1 );
    }

  /** \brief  Construct a Rank 1..8 array; use Tag#::Size for defaults. */
  Array & assign( value_type * arg_ptr )
    {
      m_array.m_ptr = arg_ptr ;
      helper::assign( m_array.m_stride );
      return *this ;
    }

  /** \brief  Construct a Rank 1..8 array; use Tag#::Size for defaults. */
  Array( value_type * arg_ptr )
    : m_array()
    {
      m_array.m_rank = Rank ;
      helper::assign_tags( m_array.m_tag );
      assign( arg_ptr );
    }

  /**  \brief  Construct compile-time array from run-time array */
  Array( const Array<Scalar,array_order> & rhs )
    : m_array( rhs )
    {
      if ( ! helper::verify( m_array.m_rank , m_array.m_tag ) ) {
        m_array.m_rank = Rank ;
        helper::assign_tags( m_array.m_tag );
        array_traits::throw_bad_conversion( m_array.m_rank , m_array.m_tag ,
                                            rhs.m_rank , rhs.m_tag );
      }
    }

  /** \brief  Return internal runtime implementation of the array */
  operator const Array<Scalar,array_order> & () const { return m_array ; }

  /** \brief  Return constructed reversed-ordering array */
  operator typename Array<Scalar,array_order>::ReverseType () const
    { return typename Array<Scalar,array_order>::ReverseType( m_array ); }

  /** \brief  Assign stride and pointer */
  void assign_stride( value_type      * arg_ptr ,
                      const size_type * arg_stride )
    {
      m_array.m_ptr = arg_ptr ;
      Copy<Rank>(   m_array.m_stride , arg_stride );
      Copy<8-Rank>( m_array.m_stride + Rank , size_type(0) );
    }

  /** \brief  Assign stride and pointer */
  void assign_stride( value_type      * arg_ptr ,
                      const size_type * arg_stride ,
                            size_type   arg_final_dim )
    {
      m_array.m_ptr = arg_ptr ;
      Copy<Rank-1>( m_array.m_stride , arg_stride );
      m_array.m_stride[Rank-1] = m_array.m_stride[Rank-2] * arg_final_dim ;
      Copy<8-Rank>( m_array.m_stride + Rank , size_type(0) );
    }

  /** \} */
protected:

  Array( const Array<Scalar,array_order> & rhs , size_type i )
    : m_array( rhs , i )
    {
      if ( ! helper::verify( m_array.m_rank , m_array.m_tag ) ) {
        m_array.m_rank = Rank ;
        helper::assign_tags( m_array.m_tag );
        array_traits::throw_bad_conversion( m_array.m_rank , m_array.m_tag ,
                                            rhs.m_rank - 1 , rhs.m_tag );
      }
    }

  Array<value_type,array_order> m_array ;

  template< typename , ArrayOrder ,
            class , class , class , class ,
            class , class , class , class >
  friend class shards::Array ;
};

//----------------------------------------------------------------------

/** \brief  Specialization for an array with Rank = 0.
 *  \nosubgrouping
 */
template< typename Scalar >
class Array<Scalar,RankZero,void,void,void,void,void,void,void,void>
{
public:
  /** \name Array Attributes
   *  \{
   */

  /** \brief  Type of member data. */
  typedef Scalar  value_type ;

  /** \brief  Type for sizes. */
  typedef array_traits::int_t size_type ;

  /** \brief  Type of runtime dimension tags. */
  typedef const ArrayDimTag * tag_type ;

  //----------------------------------

  /** \brief  Rank of the array is the number of non-void dimension tags. */
  enum { Rank = 0 };

  /** \brief  If the multidimension follows the natural ordering */
  enum { Natural = false };

  /** \brief  If the multidimension follows the reverse (Fortran) ordering */
  enum { Reverse = false };

  /** \brief  If the member data storage is contiguous */
  enum { Contiguous = true };

  /** \brief  Rank of the array is the number of non-void dimension tags. */
  size_type rank()   const { return Rank ; }

  /** \brief  If the multidimension follows the natural ordering */
  bool natural()    const { return Natural ; }

  /** \brief  If the multidimension follows the reverse (Fortran) ordering */
  bool reverse()    const { return Reverse ; }

  /** \brief  If the member data storage is contiguous */
  bool contiguous() const { return Contiguous ; }

  //----------------------------------

  /** \brief  Total number of member data items. */
  size_type size() const { return 1 ; }

  /** \} */
  //----------------------------------
  /** \name Member data access operators
   *  \{
   */

  /** \brief Pointer to contiguous block of member data. */
  value_type * contiguous_data() const { return m_ptr ; }

  /** \brief Access member via Rank 0 multi-index */
  value_type & operator()() const { return *m_ptr ; }

  /** \} */
  //----------------------------------
  /** \name Constructors and Assignment Operators
   * \{
   */

  Array() : m_ptr(NULL) {}

  Array( const Array & rhs ) : m_ptr( rhs.m_ptr ) {}

  Array & operator = ( const Array & rhs )
    { m_ptr = rhs.m_ptr ; return *this ; }

  //----------------------------------
  // Class specific constructors:

  Array( value_type * arg_ptr ) : m_ptr( arg_ptr ) {}

  /** \} */
protected:

#ifndef DOXYGEN_COMPILE
  value_type * m_ptr ;

  template< typename , ArrayOrder ,
            class , class , class , class ,
            class , class , class , class >
  friend class shards::Array ;

#endif /* DOXYGEN_COMPILE */
};


//----------------------------------------------------------------------
//----------------------------------------------------------------------

/** \} */

} // namespace shards

#undef SHARDS_ARRAY_CHECK

#endif /* Shards_Array_hpp */

