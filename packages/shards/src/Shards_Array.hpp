/*------------------------------------------------------------------------*/
/*                  shards : Shared Discretization Tools                  */
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
struct ArrayDimTag {

  /** \brief Name of the tag, typically the name of the derived class. */
  virtual const char * name() const = 0 ;

  /** \brief  Given a dimension and index produce a string for output.
   *
   *          Default to converting <b> index </b> to a string.
   */
  virtual std::string to_string( unsigned dimension ,
                                 unsigned index ) const ;

  /** \brief Given a dimension and input strige produce an index.
   *
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
 */
struct ArrayDimension : public ArrayDimTag {

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
  struct ADT : shards::ArrayDimTag { \
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

/** \brief  Return the total number of members from the array stride */
unsigned array_stride_size(
  const unsigned  rank ,
  const unsigned * const stride );

/** \brief  Generate natural dimension from array stride */
void array_stride_to_natural_dimensions(
  const unsigned   rank ,
  const unsigned * const stride ,
        unsigned * const dim );

/** \brief  Generate natural indices from array stride */
void array_stride_to_natural_indices(
  const unsigned   rank ,
  const unsigned * const stride ,
  const unsigned   offset ,
        unsigned * const indices );

/** \} */

} // namespace shards

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// Private implementation details for the array

#ifndef DOXYGEN_COMPILE

namespace shards {

void array_check_rank( const unsigned rank ,
                       const unsigned test_rank );

void array_check_ordinal( const unsigned rank ,
                          const unsigned test_ordinal );

void array_check_index( const unsigned dim ,
                        const unsigned index );

void array_check_offset(  const unsigned size ,
                          const unsigned test_offset );

void array_check_indices( const bool ,
                          const unsigned rank ,
                          const unsigned * const stride ,
                          const unsigned = 0 ,
                          const unsigned = 0 ,
                          const unsigned = 0 ,
                          const unsigned = 0 ,
                          const unsigned = 0 ,
                          const unsigned = 0 ,
                          const unsigned = 0 ,
                          const unsigned = 0 );

//----------------------------------------------------------------------

template< unsigned , unsigned > struct array_check_ordinal_is_less ;

template<> struct array_check_ordinal_is_less<0,8> {};
template<> struct array_check_ordinal_is_less<1,8> {};
template<> struct array_check_ordinal_is_less<2,8> {};
template<> struct array_check_ordinal_is_less<3,8> {};
template<> struct array_check_ordinal_is_less<4,8> {};
template<> struct array_check_ordinal_is_less<5,8> {};
template<> struct array_check_ordinal_is_less<6,8> {};
template<> struct array_check_ordinal_is_less<7,8> {};

template<> struct array_check_ordinal_is_less<0,7> {};
template<> struct array_check_ordinal_is_less<1,7> {};
template<> struct array_check_ordinal_is_less<2,7> {};
template<> struct array_check_ordinal_is_less<3,7> {};
template<> struct array_check_ordinal_is_less<4,7> {};
template<> struct array_check_ordinal_is_less<5,7> {};
template<> struct array_check_ordinal_is_less<6,7> {};

template<> struct array_check_ordinal_is_less<0,6> {};
template<> struct array_check_ordinal_is_less<1,6> {};
template<> struct array_check_ordinal_is_less<2,6> {};
template<> struct array_check_ordinal_is_less<3,6> {};
template<> struct array_check_ordinal_is_less<4,6> {};
template<> struct array_check_ordinal_is_less<5,6> {};

template<> struct array_check_ordinal_is_less<0,5> {};
template<> struct array_check_ordinal_is_less<1,5> {};
template<> struct array_check_ordinal_is_less<2,5> {};
template<> struct array_check_ordinal_is_less<3,5> {};
template<> struct array_check_ordinal_is_less<4,5> {};

template<> struct array_check_ordinal_is_less<0,4> {};
template<> struct array_check_ordinal_is_less<1,4> {};
template<> struct array_check_ordinal_is_less<2,4> {};
template<> struct array_check_ordinal_is_less<3,4> {};

template<> struct array_check_ordinal_is_less<0,3> {};
template<> struct array_check_ordinal_is_less<1,3> {};
template<> struct array_check_ordinal_is_less<2,3> {};

template<> struct array_check_ordinal_is_less<0,2> {};
template<> struct array_check_ordinal_is_less<1,2> {};

template<> struct array_check_ordinal_is_less<0,1> {};

//----------------------------------------------------------------------

template< class , unsigned > struct ArrayTagAt ;

template< typename Scalar , ArrayOrder order ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
struct ArrayTagAt<Array<Scalar,order,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>,0>
{ typedef Tag1 type ; };

template< typename Scalar , ArrayOrder order ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
struct ArrayTagAt<Array<Scalar,order,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>,1>
{ typedef Tag2 type ; };

template< typename Scalar , ArrayOrder order ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
struct ArrayTagAt<Array<Scalar,order,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>,2>
{ typedef Tag3 type ; };

template< typename Scalar , ArrayOrder order ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
struct ArrayTagAt<Array<Scalar,order,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>,3>
{ typedef Tag4 type ; };

template< typename Scalar , ArrayOrder order ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
struct ArrayTagAt<Array<Scalar,order,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>,4>
{ typedef Tag5 type ; };

template< typename Scalar , ArrayOrder order ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
struct ArrayTagAt<Array<Scalar,order,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>,5>
{ typedef Tag6 type ; };

template< typename Scalar , ArrayOrder order ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
struct ArrayTagAt<Array<Scalar,order,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>,6>
{ typedef Tag7 type ; };

template< typename Scalar , ArrayOrder order ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
struct ArrayTagAt<Array<Scalar,order,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>,7>
{ typedef Tag8 type ; };

//----------------------------------------------------------------------
//----------------------------------------------------------------------

template< ArrayOrder , unsigned Rank , unsigned Ordinal = 0 >
struct ArrayStrideDim ;

template< unsigned Rank , unsigned Ordinal >
struct ArrayStrideDim<RankZero,Rank,Ordinal> {

  static unsigned dimension( const unsigned * )
    { return 0 ; }

  static unsigned dimension( const unsigned * , const unsigned & )
    { return 0 ; }
};

template< unsigned Rank >
struct ArrayStrideDim<FortranOrder,Rank,0> {
  static unsigned dimension( const unsigned * stride )
    { return stride[0]; }

  static unsigned dimension( const unsigned * stride ,
                             const unsigned & ordinal )
    { return ordinal ? stride[ordinal] / stride[ordinal-1] : stride[0] ; }
};

template< unsigned Rank >
struct ArrayStrideDim<NaturalOrder,Rank,0> {
  static unsigned dimension( const unsigned * stride ) { return stride[0]; }

  static unsigned dimension( const unsigned * stride ,
                             const unsigned & ordinal )
    {
      const unsigned i = ( Rank - 1 ) - ordinal ;
      return i ? stride[i] / stride[i-1] : stride[0] ;
    }
};

template< unsigned Rank , unsigned Ordinal >
struct ArrayStrideDim<FortranOrder,Rank,Ordinal> {
  static unsigned dimension( const unsigned * stride )
    { return stride[Ordinal] / stride[Ordinal-1]; }
};

template< unsigned Rank , unsigned Ordinal >
struct ArrayStrideDim<NaturalOrder,Rank,Ordinal> {
  static unsigned dimension( const unsigned * stride )
    {
      enum { I = ( Rank - 1 ) - Ordinal };
      return stride[I] / stride[I-1];
    }
};

//----------------------------------------------------------------------

namespace {

template< class Tag > const ArrayDimTag * array_dim_tag();

template<>
inline
const ArrayDimTag * array_dim_tag<void>() { return NULL ; }

template< class Tag >
inline
const ArrayDimTag * array_dim_tag() { return & Tag::tag(); }

template< class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
const ArrayDimTag * const * array_dim_tags()
{
  static const ArrayDimTag * t[8] =
    {
      array_dim_tag< Tag1 >() ,
      array_dim_tag< Tag2 >() ,
      array_dim_tag< Tag3 >() ,
      array_dim_tag< Tag4 >() ,
      array_dim_tag< Tag5 >() ,
      array_dim_tag< Tag6 >() ,
      array_dim_tag< Tag7 >() ,
      array_dim_tag< Tag8 >()
    };

  return t ;
}

}

//----------------------------------------------------------------------

template< ArrayOrder array_order , unsigned rank >
unsigned array_offset(
  const unsigned * const ,
  const unsigned & i1 , const unsigned & i2 ,
  const unsigned & i3 , const unsigned & i4 ,
  const unsigned & i5 , const unsigned & i6 ,
  const unsigned & i7 , const unsigned & i8 );

template< ArrayOrder array_order , unsigned rank >
unsigned array_offset(
  const unsigned * const ,
  const unsigned & i1 , const unsigned & i2 ,
  const unsigned & i3 , const unsigned & i4 ,
  const unsigned & i5 , const unsigned & i6 ,
  const unsigned & i7 );

template< ArrayOrder array_order , unsigned rank >
unsigned array_offset(
  const unsigned * const ,
  const unsigned & i1 , const unsigned & i2 ,
  const unsigned & i3 , const unsigned & i4 ,
  const unsigned & i5 , const unsigned & i6 );

template< ArrayOrder array_order , unsigned rank >
unsigned array_offset(
  const unsigned * const ,
  const unsigned & i1 , const unsigned & i2 ,
  const unsigned & i3 , const unsigned & i4 ,
  const unsigned & i5 );

template< ArrayOrder array_order , unsigned rank >
unsigned array_offset(
  const unsigned * const ,
  const unsigned & i1 , const unsigned & i2 ,
  const unsigned & i3 , const unsigned & i4 );

template< ArrayOrder array_order , unsigned rank >
unsigned array_offset(
  const unsigned * const ,
  const unsigned & i1 , const unsigned & i2 ,
  const unsigned & i3 );

template< ArrayOrder array_order , unsigned rank >
unsigned array_offset(
  const unsigned * const ,
  const unsigned & i1 , const unsigned & i2 );

template< ArrayOrder array_order , unsigned rank >
unsigned array_offset(
  const unsigned * const ,
  const unsigned & i1 );

//----------------------------------------------------------------------

template<> inline
unsigned array_offset<FortranOrder,8>(
  const unsigned * const stride ,
  const unsigned & i1 , const unsigned & i2 ,
  const unsigned & i3 , const unsigned & i4 ,
  const unsigned & i5 , const unsigned & i6 ,
  const unsigned & i7 , const unsigned & i8 )
{
  SHARDS_ARRAY_CHECK(array_check_indices(false,8,stride,i1,i2,i3,i4,i5,i6,i7,i8));
  return i1             + i2 * stride[0] +
         i3 * stride[1] + i4 * stride[2] +
         i5 * stride[3] + i6 * stride[4] +
         i7 * stride[5] + i8 * stride[6] ;
}

template<> inline
unsigned array_offset<FortranOrder,7>(
  const unsigned * const stride ,
  const unsigned & i1 , const unsigned & i2 ,
  const unsigned & i3 , const unsigned & i4 ,
  const unsigned & i5 , const unsigned & i6 ,
  const unsigned & i7 )
{
  SHARDS_ARRAY_CHECK(array_check_indices(false,7,stride,i1,i2,i3,i4,i5,i6,i7));
  return i1             + i2 * stride[0] +
         i3 * stride[1] + i4 * stride[2] +
         i5 * stride[3] + i6 * stride[4] +
         i7 * stride[5] ;
}

template<> inline
unsigned array_offset<FortranOrder,6>(
  const unsigned * const stride ,
  const unsigned & i1 , const unsigned & i2 ,
  const unsigned & i3 , const unsigned & i4 ,
  const unsigned & i5 , const unsigned & i6 )
{
  SHARDS_ARRAY_CHECK(array_check_indices(false,6,stride,i1,i2,i3,i4,i5,i6));
  return i1             + i2 * stride[0] +
         i3 * stride[1] + i4 * stride[2] +
         i5 * stride[3] + i6 * stride[4] ;
}

template<> inline
unsigned array_offset<FortranOrder,5>(
  const unsigned * const stride ,
  const unsigned & i1 , const unsigned & i2 ,
  const unsigned & i3 , const unsigned & i4 ,
  const unsigned & i5 )
{
  SHARDS_ARRAY_CHECK(array_check_indices(false,5,stride,i1,i2,i3,i4,i5));
  return i1             + i2 * stride[0] +
         i3 * stride[1] + i4 * stride[2] +
         i5 * stride[3] ;
}

template<> inline
unsigned array_offset<FortranOrder,4>(
  const unsigned * const stride ,
  const unsigned & i1 , const unsigned & i2 ,
  const unsigned & i3 , const unsigned & i4 )
{
  SHARDS_ARRAY_CHECK(array_check_indices(false,4,stride,i1,i2,i3,i4));
  return i1             + i2 * stride[0] +
         i3 * stride[1] + i4 * stride[2] ;
}

template<> inline
unsigned array_offset<FortranOrder,3>(
  const unsigned * const stride ,
  const unsigned & i1 , const unsigned & i2 ,
  const unsigned & i3 )
{
  SHARDS_ARRAY_CHECK(array_check_indices(false,3,stride,i1,i2,i3));
  return i1 + i2 * stride[0] + i3 * stride[1] ;
}

template<> inline
unsigned array_offset<FortranOrder,2>(
  const unsigned * const stride ,
  const unsigned & i1 , const unsigned & i2 )
{
  SHARDS_ARRAY_CHECK(array_check_indices(false,2,stride,i1,i2));
  return i1 + i2 * stride[0] ;
}

template<> inline
unsigned array_offset<FortranOrder,1>(
  const unsigned * const SHARDS_ARRAY_CHECK( stride ) ,
  const unsigned & i1 )
{
  SHARDS_ARRAY_CHECK(array_check_indices(false,1,stride,i1));
  return i1 ;
}

//----------------------------------------------------------------------

template<> inline
unsigned array_offset<NaturalOrder,8>(
  const unsigned * const stride ,
  const unsigned & i1 , const unsigned & i2 ,
  const unsigned & i3 , const unsigned & i4 ,
  const unsigned & i5 , const unsigned & i6 ,
  const unsigned & i7 , const unsigned & i8 )
{
  SHARDS_ARRAY_CHECK(array_check_indices(true,8,stride,i1,i2,i3,i4,i5,i6,i7,i8));
  return i8             + i7 * stride[0] +
         i6 * stride[1] + i5 * stride[2] +
         i4 * stride[3] + i3 * stride[4] +
         i2 * stride[5] + i1 * stride[6] ;
}

template<> inline
unsigned array_offset<NaturalOrder,7>(
  const unsigned * const stride ,
  const unsigned & i1 , const unsigned & i2 ,
  const unsigned & i3 , const unsigned & i4 ,
  const unsigned & i5 , const unsigned & i6 ,
  const unsigned & i7 )
{
  SHARDS_ARRAY_CHECK(array_check_indices(true,7,stride,i1,i2,i3,i4,i5,i6,i7));
  return i7             + i6 * stride[0] +
         i5 * stride[1] + i4 * stride[2] +
         i3 * stride[3] + i2 * stride[4] +
         i1 * stride[5] ;
}

template<> inline
unsigned array_offset<NaturalOrder,6>(
  const unsigned * const stride ,
  const unsigned & i1 , const unsigned & i2 ,
  const unsigned & i3 , const unsigned & i4 ,
  const unsigned & i5 , const unsigned & i6 )
{
  SHARDS_ARRAY_CHECK(array_check_indices(true,6,stride,i1,i2,i3,i4,i5,i6));
  return i6             + i5 * stride[0] +
         i4 * stride[1] + i3 * stride[2] +
         i2 * stride[3] + i1 * stride[4] ;
}

template<> inline
unsigned array_offset<NaturalOrder,5>(
  const unsigned * const stride ,
  const unsigned & i1 , const unsigned & i2 ,
  const unsigned & i3 , const unsigned & i4 ,
  const unsigned & i5 )
{
  SHARDS_ARRAY_CHECK(array_check_indices(true,5,stride,i1,i2,i3,i4,i5));
  return i5             + i4 * stride[0] +
         i3 * stride[1] + i2 * stride[2] +
         i1 * stride[3] ;
}

template<> inline
unsigned array_offset<NaturalOrder,4>(
  const unsigned * const stride ,
  const unsigned & i1 , const unsigned & i2 ,
  const unsigned & i3 , const unsigned & i4 )
{
  SHARDS_ARRAY_CHECK(array_check_indices(true,4,stride,i1,i2,i3,i4));
  return i4             + i3 * stride[0] +
         i2 * stride[1] + i1 * stride[2] ;
}

template<> inline
unsigned array_offset<NaturalOrder,3>(
  const unsigned * const stride ,
  const unsigned & i1 , const unsigned & i2 ,
  const unsigned & i3 )
{
  SHARDS_ARRAY_CHECK(array_check_indices(true,3,stride,i1,i2,i3));
  return i3 + i2 * stride[0] + i1 * stride[1] ;
}

template<> inline
unsigned array_offset<NaturalOrder,2>(
  const unsigned * const stride ,
  const unsigned & i1 , const unsigned & i2 )
{
  SHARDS_ARRAY_CHECK(array_check_indices(true,2,stride,i1,i2));
  return i2 + i1 * stride[0] ;
}

template<> inline
unsigned array_offset<NaturalOrder,1>(
  const unsigned * const SHARDS_ARRAY_CHECK( stride ) ,
  const unsigned & i1 )
{
  SHARDS_ARRAY_CHECK(array_check_indices(true,1,stride,i1));
  return i1 ;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

template< typename Scalar , ArrayOrder Order ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
struct ArrayHelp ;

//----------------------------------------------------------------------
/** \brief  Rank 8 array types */

template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
struct ArrayHelp<Scalar,NaturalOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>
{
  typedef
    Array<Scalar,FortranOrder,Tag8,Tag7,Tag6,Tag5,Tag4,Tag3,Tag2,Tag1>
      reverse ;

  typedef
    Array<Scalar,NaturalOrder,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8,void>
      truncate ;

  enum { Rank = 8 };

  static void assign( unsigned * stride )
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

  static void assign( unsigned * stride ,
                      const unsigned & n1 )
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

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 )
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

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 ,
                      const unsigned & n3 )
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

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 ,
                      const unsigned & n3 ,
                      const unsigned & n4 )
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

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 ,
                      const unsigned & n3 ,
                      const unsigned & n4 ,
                      const unsigned & n5 )
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

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 ,
                      const unsigned & n3 ,
                      const unsigned & n4 ,
                      const unsigned & n5 ,
                      const unsigned & n6 )
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

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 ,
                      const unsigned & n3 ,
                      const unsigned & n4 ,
                      const unsigned & n5 ,
                      const unsigned & n6 ,
                      const unsigned & n7 )
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

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 ,
                      const unsigned & n3 ,
                      const unsigned & n4 ,
                      const unsigned & n5 ,
                      const unsigned & n6 ,
                      const unsigned & n7 ,
                      const unsigned & n8 )
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

  static void assign( unsigned * stride ,
                      const unsigned * const dims )
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
struct ArrayHelp<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>
{
  typedef
    Array<Scalar,NaturalOrder,Tag8,Tag7,Tag6,Tag5,Tag4,Tag3,Tag2,Tag1>
      reverse ;

  typedef
    Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,void>
      truncate ;

  enum { Rank = 8 };

  static void assign( unsigned * stride )
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

  static void assign( unsigned * stride ,
                      const unsigned & n8 )
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

  static void assign( unsigned * stride ,
                      const unsigned & n7 ,
                      const unsigned & n8 )
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

  static void assign( unsigned * stride ,
                      const unsigned & n6 ,
                      const unsigned & n7 ,
                      const unsigned & n8 )
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

  static void assign( unsigned * stride ,
                      const unsigned & n5 ,
                      const unsigned & n6 ,
                      const unsigned & n7 ,
                      const unsigned & n8 )
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

  static void assign( unsigned * stride ,
                      const unsigned & n4 ,
                      const unsigned & n5 ,
                      const unsigned & n6 ,
                      const unsigned & n7 ,
                      const unsigned & n8 )
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

  static void assign( unsigned * stride ,
                      const unsigned & n3 ,
                      const unsigned & n4 ,
                      const unsigned & n5 ,
                      const unsigned & n6 ,
                      const unsigned & n7 ,
                      const unsigned & n8 )
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

  static void assign( unsigned * stride ,
                      const unsigned & n2 ,
                      const unsigned & n3 ,
                      const unsigned & n4 ,
                      const unsigned & n5 ,
                      const unsigned & n6 ,
                      const unsigned & n7 ,
                      const unsigned & n8 )
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

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 ,
                      const unsigned & n3 ,
                      const unsigned & n4 ,
                      const unsigned & n5 ,
                      const unsigned & n6 ,
                      const unsigned & n7 ,
                      const unsigned & n8 )
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

  static void assign( unsigned * stride ,
                      const unsigned * const dims )
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
struct ArrayHelp<Scalar,NaturalOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,void>
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

  static void assign( unsigned * stride )
    {
        stride[6] = Tag1::Size * (
        stride[5] = Tag2::Size * (
        stride[4] = Tag3::Size * (
        stride[3] = Tag4::Size * (
        stride[2] = Tag5::Size * (
        stride[1] = Tag6::Size * (
        stride[0] = Tag7::Size ))))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 )
    {
        stride[6] = n1 * (
        stride[5] = Tag2::Size * (
        stride[4] = Tag3::Size * (
        stride[3] = Tag4::Size * (
        stride[2] = Tag5::Size * (
        stride[1] = Tag6::Size * (
        stride[0] = Tag7::Size ))))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 )
    {
        stride[6] = n1 * (
        stride[5] = n2 * (
        stride[4] = Tag3::Size * (
        stride[3] = Tag4::Size * (
        stride[2] = Tag5::Size * (
        stride[1] = Tag6::Size * (
        stride[0] = Tag7::Size ))))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 ,
                      const unsigned & n3 )
    {
        stride[6] = n1 * (
        stride[5] = n2 * (
        stride[4] = n3 * (
        stride[3] = Tag4::Size * (
        stride[2] = Tag5::Size * (
        stride[1] = Tag6::Size * (
        stride[0] = Tag7::Size ))))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 ,
                      const unsigned & n3 ,
                      const unsigned & n4 )
    {
        stride[6] = n1 * (
        stride[5] = n2 * (
        stride[4] = n3 * (
        stride[3] = n4 * (
        stride[2] = Tag5::Size * (
        stride[1] = Tag6::Size * (
        stride[0] = Tag7::Size ))))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 ,
                      const unsigned & n3 ,
                      const unsigned & n4 ,
                      const unsigned & n5 )
    {
        stride[6] = n1 * (
        stride[5] = n2 * (
        stride[4] = n3 * (
        stride[3] = n4 * (
        stride[2] = n5 * (
        stride[1] = Tag6::Size * (
        stride[0] = Tag7::Size ))))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 ,
                      const unsigned & n3 ,
                      const unsigned & n4 ,
                      const unsigned & n5 ,
                      const unsigned & n6 )
    {
        stride[6] = n1 * (
        stride[5] = n2 * (
        stride[4] = n3 * (
        stride[3] = n4 * (
        stride[2] = n5 * (
        stride[1] = n6 * (
        stride[0] = Tag7::Size ))))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 ,
                      const unsigned & n3 ,
                      const unsigned & n4 ,
                      const unsigned & n5 ,
                      const unsigned & n6 ,
                      const unsigned & n7 )
    {
        stride[6] = n1 * (
        stride[5] = n2 * (
        stride[4] = n3 * (
        stride[3] = n4 * (
        stride[2] = n5 * (
        stride[1] = n6 * (
        stride[0] = n7 ))))));
    }

  static void assign( unsigned * stride ,
                      const unsigned * const dims )
    {
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
struct ArrayHelp<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,void>
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

  static void assign( unsigned * stride )
    {
        stride[6] = Tag7::Size * (
        stride[5] = Tag6::Size * (
        stride[4] = Tag5::Size * (
        stride[3] = Tag4::Size * (
        stride[2] = Tag3::Size * (
        stride[1] = Tag2::Size * (
        stride[0] = Tag1::Size ))))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n7 )
    {
        stride[6] = n7 * (
        stride[5] = Tag6::Size * (
        stride[4] = Tag5::Size * (
        stride[3] = Tag4::Size * (
        stride[2] = Tag3::Size * (
        stride[1] = Tag2::Size * (
        stride[0] = Tag1::Size ))))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n6 ,
                      const unsigned & n7 )
    {
        stride[6] = n7 * (
        stride[5] = n6 * (
        stride[4] = Tag5::Size * (
        stride[3] = Tag4::Size * (
        stride[2] = Tag3::Size * (
        stride[1] = Tag2::Size * (
        stride[0] = Tag1::Size ))))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n5 ,
                      const unsigned & n6 ,
                      const unsigned & n7 )
    {
        stride[6] = n7 * (
        stride[5] = n6 * (
        stride[4] = n5 * (
        stride[3] = Tag4::Size * (
        stride[2] = Tag3::Size * (
        stride[1] = Tag2::Size * (
        stride[0] = Tag1::Size ))))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n4 ,
                      const unsigned & n5 ,
                      const unsigned & n6 ,
                      const unsigned & n7 )
    {
        stride[6] = n7 * (
        stride[5] = n6 * (
        stride[4] = n5 * (
        stride[3] = n4 * (
        stride[2] = Tag3::Size * (
        stride[1] = Tag2::Size * (
        stride[0] = Tag1::Size ))))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n3 ,
                      const unsigned & n4 ,
                      const unsigned & n5 ,
                      const unsigned & n6 ,
                      const unsigned & n7 )
    {
        stride[6] = n7 * (
        stride[5] = n6 * (
        stride[4] = n5 * (
        stride[3] = n4 * (
        stride[2] = n3 * (
        stride[1] = Tag2::Size * (
        stride[0] = Tag1::Size ))))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n2 ,
                      const unsigned & n3 ,
                      const unsigned & n4 ,
                      const unsigned & n5 ,
                      const unsigned & n6 ,
                      const unsigned & n7 )
    {
        stride[6] = n7 * (
        stride[5] = n6 * (
        stride[4] = n5 * (
        stride[3] = n4 * (
        stride[2] = n3 * (
        stride[1] = n2 * (
        stride[0] = Tag1::Size ))))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 ,
                      const unsigned & n3 ,
                      const unsigned & n4 ,
                      const unsigned & n5 ,
                      const unsigned & n6 ,
                      const unsigned & n7 )
    {
        stride[6] = n7 * (
        stride[5] = n6 * (
        stride[4] = n5 * (
        stride[3] = n4 * (
        stride[2] = n3 * (
        stride[1] = n2 * (
        stride[0] = n1 ))))));
    }

  static void assign( unsigned * stride ,
                      const unsigned * const dims )
    {
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
struct ArrayHelp<Scalar,NaturalOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,void,void>
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

  static void assign( unsigned * stride )
    {
        stride[5] = Tag1::Size * (
        stride[4] = Tag2::Size * (
        stride[3] = Tag3::Size * (
        stride[2] = Tag4::Size * (
        stride[1] = Tag5::Size * (
        stride[0] = Tag6::Size )))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 )
    {
        stride[5] = n1 * (
        stride[4] = Tag2::Size * (
        stride[3] = Tag3::Size * (
        stride[2] = Tag4::Size * (
        stride[1] = Tag5::Size * (
        stride[0] = Tag6::Size )))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 )
    {
        stride[5] = n1 * (
        stride[4] = n2 * (
        stride[3] = Tag3::Size * (
        stride[2] = Tag4::Size * (
        stride[1] = Tag5::Size * (
        stride[0] = Tag6::Size )))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 ,
                      const unsigned & n3 )
    {
        stride[5] = n1 * (
        stride[4] = n2 * (
        stride[3] = n3 * (
        stride[2] = Tag4::Size * (
        stride[1] = Tag5::Size * (
        stride[0] = Tag6::Size )))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 ,
                      const unsigned & n3 ,
                      const unsigned & n4 )
    {
        stride[5] = n1 * (
        stride[4] = n2 * (
        stride[3] = n3 * (
        stride[2] = n4 * (
        stride[1] = Tag5::Size * (
        stride[0] = Tag6::Size )))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 ,
                      const unsigned & n3 ,
                      const unsigned & n4 ,
                      const unsigned & n5 )
    {
        stride[5] = n1 * (
        stride[4] = n2 * (
        stride[3] = n3 * (
        stride[2] = n4 * (
        stride[1] = n5 * (
        stride[0] = Tag6::Size )))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 ,
                      const unsigned & n3 ,
                      const unsigned & n4 ,
                      const unsigned & n5 ,
                      const unsigned & n6 )
    {
        stride[5] = n1 * (
        stride[4] = n2 * (
        stride[3] = n3 * (
        stride[2] = n4 * (
        stride[1] = n5 * (
        stride[0] = n6 )))));
    }

  static void assign( unsigned * stride ,
                      const unsigned * const dims )
    {
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
struct ArrayHelp<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,void,void>
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

  static void assign( unsigned * stride )
    {
        stride[5] = Tag6::Size * (
        stride[4] = Tag5::Size * (
        stride[3] = Tag4::Size * (
        stride[2] = Tag3::Size * (
        stride[1] = Tag2::Size * (
        stride[0] = Tag1::Size )))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n6 )
    {
        stride[5] = n6 * (
        stride[4] = Tag5::Size * (
        stride[3] = Tag4::Size * (
        stride[2] = Tag3::Size * (
        stride[1] = Tag2::Size * (
        stride[0] = Tag1::Size )))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n5 ,
                      const unsigned & n6 )
    {
        stride[5] = n6 * (
        stride[4] = n5 * (
        stride[3] = Tag4::Size * (
        stride[2] = Tag3::Size * (
        stride[1] = Tag2::Size * (
        stride[0] = Tag1::Size )))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n4 ,
                      const unsigned & n5 ,
                      const unsigned & n6 )
    {
        stride[5] = n6 * (
        stride[4] = n5 * (
        stride[3] = n4 * (
        stride[2] = Tag3::Size * (
        stride[1] = Tag2::Size * (
        stride[0] = Tag1::Size )))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n3 ,
                      const unsigned & n4 ,
                      const unsigned & n5 ,
                      const unsigned & n6 )
    {
        stride[5] = n6 * (
        stride[4] = n5 * (
        stride[3] = n4 * (
        stride[2] = n3 * (
        stride[1] = Tag2::Size * (
        stride[0] = Tag1::Size )))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n2 ,
                      const unsigned & n3 ,
                      const unsigned & n4 ,
                      const unsigned & n5 ,
                      const unsigned & n6 )
    {
        stride[5] = n6 * (
        stride[4] = n5 * (
        stride[3] = n4 * (
        stride[2] = n3 * (
        stride[1] = n2 * (
        stride[0] = Tag1::Size )))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 ,
                      const unsigned & n3 ,
                      const unsigned & n4 ,
                      const unsigned & n5 ,
                      const unsigned & n6 )
    {
        stride[5] = n6 * (
        stride[4] = n5 * (
        stride[3] = n4 * (
        stride[2] = n3 * (
        stride[1] = n2 * (
        stride[0] = n1 )))));
    }

  static void assign( unsigned * stride ,
                      const unsigned * const dims )
    {
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
struct ArrayHelp<Scalar,NaturalOrder,Tag1,Tag2,Tag3,Tag4,Tag5,void,void,void>
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

  static void assign( unsigned * stride )
    {
        stride[4] = Tag1::Size * (
        stride[3] = Tag2::Size * (
        stride[2] = Tag3::Size * (
        stride[1] = Tag4::Size * (
        stride[0] = Tag5::Size ))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 )
    {
        stride[4] = n1 * (
        stride[3] = Tag2::Size * (
        stride[2] = Tag3::Size * (
        stride[1] = Tag4::Size * (
        stride[0] = Tag5::Size ))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 )
    {
        stride[4] = n1 * (
        stride[3] = n2 * (
        stride[2] = Tag3::Size * (
        stride[1] = Tag4::Size * (
        stride[0] = Tag5::Size ))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 ,
                      const unsigned & n3 )
    {
        stride[4] = n1 * (
        stride[3] = n2 * (
        stride[2] = n3 * (
        stride[1] = Tag4::Size * (
        stride[0] = Tag5::Size ))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 ,
                      const unsigned & n3 ,
                      const unsigned & n4 )
    {
        stride[4] = n1 * (
        stride[3] = n2 * (
        stride[2] = n3 * (
        stride[1] = n4 * (
        stride[0] = Tag5::Size ))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 ,
                      const unsigned & n3 ,
                      const unsigned & n4 ,
                      const unsigned & n5 )
    {
        stride[4] = n1 * (
        stride[3] = n2 * (
        stride[2] = n3 * (
        stride[1] = n4 * (
        stride[0] = n5 ))));
    }

  static void assign( unsigned * stride ,
                      const unsigned * const dims )
    {
        stride[4] = dims[0] * (
        stride[3] = dims[1] * (
        stride[2] = dims[2] * (
        stride[1] = dims[3] * (
        stride[0] = dims[4] ))));
    }
};

template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 , class Tag5 >
struct ArrayHelp<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,Tag5,void,void,void>
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

  static void assign( unsigned * stride )
    {
        stride[4] = Tag5::Size * (
        stride[3] = Tag4::Size * (
        stride[2] = Tag3::Size * (
        stride[1] = Tag2::Size * (
        stride[0] = Tag1::Size ))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n5 )
    {
        stride[4] = n5 * (
        stride[3] = Tag4::Size * (
        stride[2] = Tag3::Size * (
        stride[1] = Tag2::Size * (
        stride[0] = Tag1::Size ))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n4 ,
                      const unsigned & n5 )
    {
        stride[4] = n5 * (
        stride[3] = n4 * (
        stride[2] = Tag3::Size * (
        stride[1] = Tag2::Size * (
        stride[0] = Tag1::Size ))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n3 ,
                      const unsigned & n4 ,
                      const unsigned & n5 )
    {
        stride[4] = n5 * (
        stride[3] = n4 * (
        stride[2] = n3 * (
        stride[1] = Tag2::Size * (
        stride[0] = Tag1::Size ))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n2 ,
                      const unsigned & n3 ,
                      const unsigned & n4 ,
                      const unsigned & n5 )
    {
        stride[4] = n5 * (
        stride[3] = n4 * (
        stride[2] = n3 * (
        stride[1] = n2 * (
        stride[0] = Tag1::Size ))));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 ,
                      const unsigned & n3 ,
                      const unsigned & n4 ,
                      const unsigned & n5 )
    {
        stride[4] = n5 * (
        stride[3] = n4 * (
        stride[2] = n3 * (
        stride[1] = n2 * (
        stride[0] = n1 ))));
    }

  static void assign( unsigned * stride ,
                      const unsigned * const dims )
    {
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
struct ArrayHelp<Scalar,NaturalOrder,Tag1,Tag2,Tag3,Tag4,void,void,void,void>
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

  static void assign( unsigned * stride )
    {
        stride[3] = Tag1::Size * (
        stride[2] = Tag2::Size * (
        stride[1] = Tag3::Size * (
        stride[0] = Tag4::Size )));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 )
    {
        stride[3] = n1 * (
        stride[2] = Tag2::Size * (
        stride[1] = Tag3::Size * (
        stride[0] = Tag4::Size )));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 )
    {
        stride[3] = n1 * (
        stride[2] = n2 * (
        stride[1] = Tag3::Size * (
        stride[0] = Tag4::Size )));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 ,
                      const unsigned & n3 )
    {
        stride[3] = n1 * (
        stride[2] = n2 * (
        stride[1] = n3 * (
        stride[0] = Tag4::Size )));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 ,
                      const unsigned & n3 ,
                      const unsigned & n4 )
    {
        stride[3] = n1 * (
        stride[2] = n2 * (
        stride[1] = n3 * (
        stride[0] = n4 )));
    }

  static void assign( unsigned * stride ,
                      const unsigned * const dims )
    {
        stride[3] = dims[0] * (
        stride[2] = dims[1] * (
        stride[1] = dims[2] * (
        stride[0] = dims[3] )));
    }
};

template< typename Scalar , class Tag1 , class Tag2 , class Tag3 , class Tag4 >
struct ArrayHelp<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,void,void,void,void>
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

  static void assign( unsigned * stride )
    {
        stride[3] = Tag4::Size * (
        stride[2] = Tag3::Size * (
        stride[1] = Tag2::Size * (
        stride[0] = Tag1::Size )));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n4 )
    {
        stride[3] = n4 * (
        stride[2] = Tag3::Size * (
        stride[1] = Tag2::Size * (
        stride[0] = Tag1::Size )));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n3 ,
                      const unsigned & n4 )
    {
        stride[3] = n4 * (
        stride[2] = n3 * (
        stride[1] = Tag2::Size * (
        stride[0] = Tag1::Size )));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n2 ,
                      const unsigned & n3 ,
                      const unsigned & n4 )
    {
        stride[3] = n4 * (
        stride[2] = n3 * (
        stride[1] = n2 * (
        stride[0] = Tag1::Size )));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 ,
                      const unsigned & n3 ,
                      const unsigned & n4 )
    {
        stride[3] = n4 * (
        stride[2] = n3 * (
        stride[1] = n2 * (
        stride[0] = n1 )));
    }

  static void assign( unsigned * stride ,
                      const unsigned * const dims )
    {
        stride[3] = dims[3] * (
        stride[2] = dims[2] * (
        stride[1] = dims[1] * (
        stride[0] = dims[0] )));
   }
};

//----------------------------------------------------------------------
/** \brief  Rank 3 array types */

template< typename Scalar , class Tag1 , class Tag2 , class Tag3 >
struct ArrayHelp<Scalar,NaturalOrder,Tag1,Tag2,Tag3,void,void,void,void,void>
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

  static void assign( unsigned * stride )
    {
        stride[2] = Tag1::Size * (
        stride[1] = Tag2::Size * (
        stride[0] = Tag3::Size ));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 )
    {
        stride[2] = n1 * (
        stride[1] = Tag2::Size * (
        stride[0] = Tag3::Size ));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 )
    {
        stride[2] = n1 * (
        stride[1] = n2 * (
        stride[0] = Tag3::Size ));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 ,
                      const unsigned & n3 )
    {
        stride[2] = n1 * (
        stride[1] = n2 * (
        stride[0] = n3 ));
    }

  static void assign( unsigned * stride ,
                      const unsigned * const dims )
    {
        stride[2] = dims[0] * (
        stride[1] = dims[1] * (
        stride[0] = dims[2] ));
    }
};

template< typename Scalar , class Tag1 , class Tag2 , class Tag3 >
struct ArrayHelp<Scalar,FortranOrder,Tag1,Tag2,Tag3,void,void,void,void,void>
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

  static void assign( unsigned * stride )
    {
        stride[2] = Tag3::Size * (
        stride[1] = Tag2::Size * (
        stride[0] = Tag1::Size ));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n3 )
    {
        stride[2] = n3 * (
        stride[1] = Tag2::Size * (
        stride[0] = Tag1::Size ));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n2 ,
                      const unsigned & n3 )
    {
        stride[2] = n3 * (
        stride[1] = n2 * (
        stride[0] = Tag1::Size ));
    }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 ,
                      const unsigned & n3 )
    {
        stride[2] = n3 * (
        stride[1] = n2 * (
        stride[0] = n1 ));
    }

  static void assign( unsigned * stride ,
                      const unsigned * const dims )
    {
        stride[2] = dims[2] * (
        stride[1] = dims[1] * (
        stride[0] = dims[0] ));
   }
};

//----------------------------------------------------------------------
/** \brief  Rank 2 array types */

template< typename Scalar , class Tag1 , class Tag2 >
struct ArrayHelp<Scalar,NaturalOrder,Tag1,Tag2,void,void,void,void,void,void>
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

  static void assign( unsigned * stride )
    { stride[1] = Tag1::Size * ( stride[0] = Tag2::Size ); }

  static void assign( unsigned * stride ,
                      const unsigned & n1 )
    { stride[1] = n1 * ( stride[0] = Tag2::Size ); }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 )
    { stride[1] = n1 * ( stride[0] = n2 ); }

  static void assign( unsigned * stride ,
                      const unsigned * const dims )
    {
        stride[1] = dims[0] * (
        stride[0] = dims[1] );
    }
};

template< typename Scalar , class Tag1 , class Tag2 >
struct ArrayHelp<Scalar,FortranOrder,Tag1,Tag2,void,void,void,void,void,void>
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

  static void assign( unsigned * stride )
    { stride[1] = Tag2::Size * ( stride[0] = Tag1::Size ); }

  static void assign( unsigned * stride ,
                      const unsigned & n2 )
    { stride[1] = n2 * ( stride[0] = Tag1::Size ); }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ,
                      const unsigned & n2 )
    { stride[1] = n2 * ( stride[0] = n1 ); }

  static void assign( unsigned * stride ,
                      const unsigned * const dims )
    { stride[1] = dims[1] * ( stride[0] = dims[0] ); }
};

//----------------------------------------------------------------------
/** \brief  Rank 1 array types */

template< typename Scalar , class Tag1 >
struct ArrayHelp<Scalar,NaturalOrder,Tag1,void,void,void,void,void,void,void>
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

  static void assign( unsigned * stride ) { stride[0] = Tag1::Size ; }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ) { stride[0] = n1 ; }

  static void assign( unsigned * stride ,
                      const unsigned * const dims )
    { stride[0] = dims[0] ; }
};

template< typename Scalar , class Tag1 >
struct ArrayHelp<Scalar,FortranOrder,Tag1,void,void,void,void,void,void,void>
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

  static void assign( unsigned * stride ) { stride[0] = Tag1::Size ; }

  static void assign( unsigned * stride ,
                      const unsigned & n1 ) { stride[0] = n1 ; }

  static void assign( unsigned * stride ,
                      const unsigned * const dims )
    { stride[0] = dims[0] ; }
};

//----------------------------------------------------------------------
/** \brief  Rank Zero array types */

template< typename Scalar >
struct ArrayHelp<Scalar,RankZero,void,void,void,void,void,void,void,void>
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

  static void assign( unsigned * ) {}
};

//----------------------------------------------------------------------
/** \brief  Rank Dynamic array types */

template< typename Scalar >
struct ArrayHelp<Scalar,NaturalOrder,void,void,void,void,void,void,void,void>
{
  typedef
    Array<Scalar,FortranOrder,void,void,void,void,void,void,void,void>
      reverse ;
};

template< typename Scalar >
struct ArrayHelp<Scalar,FortranOrder,void,void,void,void,void,void,void,void>
{
  typedef
    Array<Scalar,NaturalOrder,void,void,void,void,void,void,void,void>
      reverse ;
};

//----------------------------------------------------------------------

template< class ArrayType , class TagA > struct ArrayAppend {};

template< typename Scalar , ArrayOrder Order ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class TagA >
struct ArrayAppend<
       Array<Scalar,Order,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,void> , TagA >
{
private:
  typedef typename
    ArrayHelp<Scalar,Order,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,void>
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
    ArrayHelp<Scalar,array_order,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>
      help_type ;
#endif /* DOXYGEN_COMPILE */

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
  enum { Rank = help_type::Rank };

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
  typedef typename help_type::truncate TruncateType ;

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
      SHARDS_ARRAY_CHECK( array_check_offset(size(),i) );
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
  typedef typename help_type::reverse  ReverseType ;

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
    : m_ptr( arg_ptr ) { help_type::assign( m_stride , dims ); }

  /** \brief  Construct a Rank 8 array */
  Array( value_type * arg_ptr ,
         const unsigned n1 , const unsigned n2 ,
         const unsigned n3 , const unsigned n4 ,
         const unsigned n5 , const unsigned n6 ,
         const unsigned n7 , const unsigned n8 )
    : m_ptr( arg_ptr )
    { help_type::assign( m_stride , n1 , n2 , n3 , n4 , n5 , n6 , n7 , n8 ); }

  /** \brief  Construct a Rank 7..8 array; use Tag#::Size for defaults.
   *          The input dimensions are the 7 slowest strides.
   */
  Array( value_type * arg_ptr ,
         const unsigned n1 , const unsigned n2 ,
         const unsigned n3 , const unsigned n4 ,
         const unsigned n5 , const unsigned n6 ,
         const unsigned n7 )
    : m_ptr( arg_ptr )
    { help_type::assign( m_stride , n1 , n2 , n3 , n4 , n5 , n6 , n7 ); }

  /** \brief  Construct a Rank 6..8 array; use Tag#::Size for defaults.
   *          The input dimensions are the 6 slowest strides.
   */
  Array( value_type * arg_ptr ,
         const unsigned n1 , const unsigned n2 ,
         const unsigned n3 , const unsigned n4 ,
         const unsigned n5 , const unsigned n6 )
    : m_ptr( arg_ptr )
    { help_type::assign( m_stride , n1 , n2 , n3 , n4 , n5 , n6 ); }

  /** \brief  Construct a Rank 5..8 array; use Tag#::Size for defaults.
   *          The input dimensions are the 5 slowest strides.
   */
  Array( value_type * arg_ptr ,
         const unsigned n1 , const unsigned n2 ,
         const unsigned n3 , const unsigned n4 ,
         const unsigned n5 )
    : m_ptr( arg_ptr )
    { help_type::assign( m_stride , n1 , n2 , n3 , n4 , n5 ); }

  /** \brief  Construct a Rank 4..8 array; use Tag#::Size for defaults.
   *          The input dimensions are the 4 slowest strides.
   */
  Array( value_type * arg_ptr ,
         const unsigned n1 , const unsigned n2 ,
         const unsigned n3 , const unsigned n4 )
    : m_ptr( arg_ptr )
    { help_type::assign( m_stride , n1 , n2 , n3 , n4 ); }

  /** \brief  Construct a Rank 3..8 array; use Tag#::Size for defaults.
   *          The input dimensions are the 3 slowest strides.
   */
  Array( value_type * arg_ptr ,
         const unsigned n1 , const unsigned n2 ,
         const unsigned n3 )
    : m_ptr( arg_ptr )
    { help_type::assign( m_stride , n1 , n2 , n3 ); }

  /** \brief  Construct a Rank 2..8 array; use Tag#::Size for defaults.
   *          The input dimensions are the 2 slowest strides.
   */
  Array( value_type * arg_ptr , const unsigned n1 , const unsigned n2 )
    : m_ptr( arg_ptr ) { help_type::assign( m_stride , n1 , n2 ); }

  /** \brief  Construct a Rank 1..8 array; use Tag#::Size for defaults.
   *          The input dimension is the slowest stride.
   */
  Array( value_type * arg_ptr , const unsigned n1 )
    : m_ptr( arg_ptr ) { help_type::assign( m_stride , n1 ); }

  /** \brief  Construct a Rank 1..8 array; use Tag#::Size for defaults. */
  Array( value_type * arg_ptr )
    : m_ptr( arg_ptr ) { help_type::assign( m_stride ); }

  /** \} */
protected:

  /** \brief Pointer to contiguous block of members */
  value_type * m_ptr ;

  /** \brief Array of strides, smallest to largest */
  size_type m_stride[ Rank ];

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
  typedef unsigned size_type ;

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
  unsigned rank()   const { return Rank ; }

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
/** \brief  The <b> not-preferred </b> multi-dimensional Array interface
 *          with <b> runtime </b> user-defined dimension ordinates.
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
      SHARDS_ARRAY_CHECK( array_check_offset(size(),i) );
      return m_ptr[ i ];
    }

  //----------------------------------
  /** \brief Access member via Rank 8 multi-index */
  value_type & operator()( const unsigned i1 , const unsigned i2 ,
                           const unsigned i3 , const unsigned i4 ,
                           const unsigned i5 , const unsigned i6 ,
                           const unsigned i7 , const unsigned i8 ) const
    {
      SHARDS_ARRAY_CHECK( array_check_rank( m_rank , 8 ) );
      return m_ptr[
        array_offset<array_order,8>(m_stride,i1,i2,i3,i4,i5,i6,i7,i8) ];
    }

  value_type & operator()( const unsigned i1 , const unsigned i2 ,
                           const unsigned i3 , const unsigned i4 ,
                           const unsigned i5 , const unsigned i6 ,
                           const unsigned i7 ) const
    {
      SHARDS_ARRAY_CHECK( array_check_rank( m_rank , 7 ) );
      return m_ptr[
        array_offset<array_order,7>(m_stride,i1,i2,i3,i4,i5,i6,i7) ];
    }

  value_type & operator()( const unsigned i1 , const unsigned i2 ,
                           const unsigned i3 , const unsigned i4 ,
                           const unsigned i5 , const unsigned i6 ) const
    {
      SHARDS_ARRAY_CHECK( array_check_rank( m_rank , 6 ) );
      return m_ptr[ array_offset<array_order,6>(m_stride,i1,i2,i3,i4,i5,i6) ];
    }

  value_type & operator()( const unsigned i1 , const unsigned i2 ,
                           const unsigned i3 , const unsigned i4 ,
                           const unsigned i5 ) const
    {
      SHARDS_ARRAY_CHECK( array_check_rank( m_rank , 5 ) );
      return m_ptr[ array_offset<array_order,5>(m_stride,i1,i2,i3,i4,i5) ];
    }

  value_type & operator()( const unsigned i1 , const unsigned i2 ,
                           const unsigned i3 , const unsigned i4 ) const
    {
      SHARDS_ARRAY_CHECK( array_check_rank( m_rank , 4 ) );
      return m_ptr[ array_offset<array_order,4>(m_stride,i1,i2,i3,i4) ];
    }

  value_type & operator()( const unsigned i1 , const unsigned i2 ,
                           const unsigned i3 ) const
    {
      SHARDS_ARRAY_CHECK( array_check_rank( m_rank , 3 ) );
      return m_ptr[ array_offset<array_order,3>(m_stride,i1,i2,i3) ];
    }

  value_type & operator()( const unsigned i1 , const unsigned i2 ) const
    {
      SHARDS_ARRAY_CHECK( array_check_rank( m_rank , 2 ) );
      return m_ptr[ array_offset<array_order,2>(m_stride,i1,i2) ];
    }

  value_type & operator()( const unsigned i1 ) const
    {
      SHARDS_ARRAY_CHECK( array_check_rank( m_rank , 1 ) );
      return m_ptr[ array_offset<array_order,1>(m_stride,i1) ];
    }

  /** \} */
  //----------------------------------
  /** \name Constructors and Assignment Operators
   * \{
   */

  typedef typename
    ArrayHelp<Scalar,array_order,void,void,void,void,void,void,void,void>
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


  /** \brief Pointer to contiguous block of members */
  value_type * m_ptr ;

  /** \brief Rank of the array */
  unsigned     m_rank ;

  /** \brief Array of strides, smallest to largest */
  size_type    m_stride[8];

  /** \brief Array of singleton tags, aligned with strides */
  tag_type     m_tag[8] ;

  template< typename , ArrayOrder ,
            class , class , class , class ,
            class , class , class , class >
  friend class shards::Array ;

};

//----------------------------------------------------------------------
//----------------------------------------------------------------------

/** \} */

} // namespace shards

#undef SHARDS_ARRAY_CHECK

#endif /* Shards_Array_hpp */

