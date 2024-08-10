// @HEADER
// *****************************************************************************
//                Shards : Shared Discretization Tools
//
// Copyright 2008-2011 NTESS and the Shards contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef Shards_TypeList_hpp
#define Shards_TypeList_hpp

namespace shards {

/** \ingroup shards_package
 *  \defgroup shards_package_typelist  Linked List of Types
 *  \brief     Compile-time linked-list of types and operations.
 *
 *  This type list capability was inspired by the type list
 *  described in Alexandrescu's "Modern C++ Design" book.
 *
 *  \author H. Carter Edwards  <hcedwar@sandia.gov>
 *
 *  \{
 */

//----------------------------------------------------------------------
/** \class  SameType
 *  \brief  Member <b> enum { value = ... }; </b>
 *          is true if <b> T1 </b> and <b> T2 </b> are the same type.
 */
template<typename T1, typename T2>
struct SameType { enum { value = false }; };

template<typename T>
struct SameType<T,T> { enum { value = true }; };

//----------------------------------------------------------------------

struct TypeListEnd {};

/** \class  TypeList
 *  \brief  A link within a linked list of types.
 *
 *  A linked list of types where <b> Tail </b> is required to either
 *  terminate the list with <b> TypeListEnd </b> or continue the list
 *  with another instantiation of <b> TypeList </b>.
 */
template< typename Value , class Tail = TypeListEnd > struct TypeList {};

//----------------------------------------------------------------------
/** \class  TypeListLength
 *  \brief  Member <b> enum { value = ... }; </b>
 *          is the length of the type list.
 */
template< class ListType > struct TypeListLength {};

template<>
struct TypeListLength< TypeListEnd >
{ enum { value = 0 }; };

template< typename Value , class Tail >
struct TypeListLength< TypeList< Value , Tail > >
{ enum { value = 1 + TypeListLength< Tail >::value }; };

//----------------------------------------------------------------------
/** \class  TypeListAt
 *  \brief  Member <b> typedef ... type ; </b> is the type of the member
 *          of <b> ListType </b> at location <b> ordinal </b>
 *          if <b> ordinal </b> is less than the type list length.
 */
template< class ListType, unsigned ordinal > struct TypeListAt {};

template< unsigned ordinal >
struct TypeListAt< TypeListEnd , ordinal >
{ typedef TypeListEnd type ; };

template< typename Value , class Tail >
struct TypeListAt< TypeList< Value , Tail > , 0 >
{ typedef Value type ; };

template< typename Value , class Tail , unsigned ordinal >
struct TypeListAt< TypeList< Value , Tail > , ordinal >
{ typedef typename TypeListAt< Tail , ordinal - 1 >::type type ; };

//----------------------------------------------------------------------
/** \class  TypeListIndex
 *  \brief  Member <b> enum { value = ... }; </b>
 *          is the location within <b> ListType </b>
 *          of occurance <b> I </b> of type <b> TestValue </b>.
 *          If this occurance does not exist then <b> value = -1 </b>.
 */
template< class ListType , typename TestValue , unsigned ordinal = 0 >
struct TypeListIndex {};

template< typename TestValue , unsigned ordinal >
struct TypeListIndex< TypeListEnd , TestValue , ordinal >
{
  enum { value = -1 };
};

template< typename Value , class Tail , typename TestValue , unsigned ordinal >
struct TypeListIndex< TypeList< Value , Tail > , TestValue , ordinal >
{
private:
  enum { match = SameType< Value , TestValue >::value };
  enum { J = match && 0 < ordinal ? ordinal - 1 : ordinal };
  enum { N = TypeListIndex< Tail , TestValue , J >::value };
public:
  enum { value = match && 0 == ordinal ? 0 : ( -1 == N ? -1 : N + 1 ) };
};

//----------------------------------------------------------------------
/** \class  TypeListCount
 *  \brief  Member <b> enum { value = ... }; </b>
 *          is the number of occurances of <b> TestValue </b>
 *          within <b> ListType </b>.
 */
template< class ListType , typename TestValue >
struct TypeListCount {};

template< typename TestValue >
struct TypeListCount< TypeListEnd , TestValue >
{ enum { value = 0 }; };

template< typename Value , class Tail , typename TestValue >
struct TypeListCount< TypeList< Value , Tail > , TestValue >
{
  enum { value = TypeListCount< Tail , TestValue >::value +
                 ( SameType< Value , TestValue >::value ? 1 : 0 ) };
};

//----------------------------------------------------------------------
/** \class  TypeListMember
 *  \brief  Member <b> enum { value = ... }; </b> is true
 *          if <b> TestValue </b> is a member of <b> ListType </b>.
 */
template< class ListType , typename TestValue > struct TypeListMember {};

template< typename TestValue >
struct TypeListMember< TypeListEnd , TestValue >
{ enum { value = false }; };

template< typename Value , class Tail , typename TestValue >
struct TypeListMember< TypeList< Value , Tail > , TestValue >
{
  enum { value = SameType< Value , TestValue >::value ||
                 TypeListMember< Tail , TestValue >::value };
};

//----------------------------------------------------------------------
/** \class  TypeListUnique
 *  \brief  Member <b> enum { value = ... }; </b> is true
 *          if each member of <b> ListType </b> appears exactly once.
 */
template< class ListType > struct TypeListUnique {};

template<>
struct TypeListUnique< TypeListEnd >
{ enum { value = true }; };

template< typename Value , class Tail >
struct TypeListUnique< TypeList< Value , Tail > >
{
  enum { value = ! TypeListMember< Tail , Value >::value &&
                 TypeListUnique< Tail >::value };
};

//----------------------------------------------------------------------
/** \class  TypeListDisjoint
 *  \brief  Member <b> enum { value = ... }; </b> is true
 *          if all members of <b> ListA </b>
 *          are not a member <b> ListB </b>.
 */
template< class ListA , class ListB > struct TypeListDisjoint {};

template< class ListB >
struct TypeListDisjoint< TypeListEnd , ListB >
{ enum { value = true }; };

template< typename Value , class Tail , class ListB >
struct TypeListDisjoint< TypeList< Value , Tail > , ListB >
{
  enum { value = ! TypeListMember< ListB , Value >::value &&
                 TypeListDisjoint< Tail , ListB >::value };
};

//----------------------------------------------------------------------
/** \class  TypeListFirst
 *  \brief  Member <b> typedef ... type ; </b> is the first member
 *          of <b> ListType </b>.
 */
template< class ListType > struct TypeListFirst {};

template<>
struct TypeListFirst< TypeListEnd >
{ typedef TypeListEnd type ; };

template< typename Value , class Tail >
struct TypeListFirst< TypeList< Value , Tail > >
{ typedef Value type ; };

//----------------------------------------------------------------------
/** \class  TypeListLast
 *  \brief  Member <b> typedef ... type ; </b> is the last member
 *          of <b> ListType </b>.
 */
template< class ListType > struct TypeListLast {};

template<>
struct TypeListLast< TypeListEnd >
{ typedef TypeListEnd type ; };

template< typename Value >
struct TypeListLast< TypeList< Value , TypeListEnd > >
{ typedef Value type ; };

template< typename Value , class Tail >
struct TypeListLast< TypeList< Value , Tail > >
{ typedef typename TypeListLast< Tail >::type type ; };

//----------------------------------------------------------------------
/** \class  TypeListAppend
 *  \brief  Member <b> typedef ... type ; </b> is defined
 *          by appending <b> T </b> to the end of <b> ListA </b>.
 */
template< class ListA , typename T > struct TypeListAppend {};

template<>
struct TypeListAppend< TypeListEnd , TypeListEnd >
{ typedef TypeListEnd type ; };

template< typename T >
struct TypeListAppend< TypeListEnd , T >
{ typedef TypeList< T > type ; };

template< typename Value , class Tail , typename T >
struct TypeListAppend< TypeList< Value , Tail > , T >
{
  typedef TypeList< Value , typename TypeListAppend< Tail , T >::type > type ;
};

//----------------------------------------------------------------------
/** \class  TypeListJoin
 *  \brief  Member <b> typedef ... type ; </b> is defined
 *          by joining <b> ListB </b> to the end of <b> ListA </b>.
 */
template< class ListA , class ListB > struct TypeListJoin {};

template<>
struct TypeListJoin< TypeListEnd , TypeListEnd >
{ typedef TypeListEnd type ; };

template< typename Value , class Tail >
struct TypeListJoin< TypeListEnd , TypeList< Value , Tail > >
{ typedef TypeList< Value , Tail > type ; };

template< typename ValueA , class TailA , typename ValueB , class TailB >
struct TypeListJoin< TypeList< ValueA , TailA > ,
                     TypeList< ValueB , TailB > >
{
private:
  typedef typename
    TypeListJoin< TailA , TypeList< ValueB , TailB > >::type Tail ;
public:
  typedef TypeList< ValueA , Tail > type ;
};

//----------------------------------------------------------------------
/** \class  TypeListEraseAt
 *  \brief  Member <b> typedef ... type ; </b> is defined
 *          by erasing member at <b> ordinal </b> from <b> ListType </b>.
 */
template< class ListType, unsigned ordinal > struct TypeListEraseAt {};

template< typename Value , class Tail >
struct TypeListEraseAt< TypeList< Value , Tail > , 0 >
{ typedef Tail type ; };

template< typename Value , class Tail , unsigned ordinal >
struct TypeListEraseAt< TypeList< Value , Tail > , ordinal >
{
  typedef TypeList< Value ,
                    typename TypeListEraseAt<Tail,ordinal-1>::type > type ;
};

//----------------------------------------------------------------------
/** \class  TypeListClean
 *  \brief  Member <b> typedef ... type ; </b> is defined
 *          by truncating <b> ListType </b> at the first
 *          occurance of <b> TypeListEnd </b>.
 *          Used by <b> MakeTypeList </b> to generate a clean type list.
 */
template< class ListType > struct TypeListClean {};

template<>
struct TypeListClean< TypeListEnd >
{ typedef TypeListEnd type ; };

template< class Tail >
struct TypeListClean< TypeList< TypeListEnd , Tail > >
{ typedef TypeListEnd type ; };

template< typename Value , class Tail >
struct TypeListClean< TypeList< Value , Tail > >
{
  typedef TypeList< Value , typename TypeListClean< Tail >::type > type ;
};

//----------------------------------------------------------------------
/** \class  MakeTypeList
 *  \brief  Member <b> typedef ... type ; </b>
 *          is a type list constructed from the template arguments.
 */
template< typename T00 = TypeListEnd ,
          typename T01 = TypeListEnd ,
          typename T02 = TypeListEnd ,
          typename T03 = TypeListEnd ,
          typename T04 = TypeListEnd ,
          typename T05 = TypeListEnd ,
          typename T06 = TypeListEnd ,
          typename T07 = TypeListEnd ,
          typename T08 = TypeListEnd ,
          typename T09 = TypeListEnd ,
          typename T10 = TypeListEnd ,
          typename T11 = TypeListEnd ,
          typename T12 = TypeListEnd ,
          typename T13 = TypeListEnd ,
          typename T14 = TypeListEnd ,
          typename T15 = TypeListEnd ,
          typename T16 = TypeListEnd ,
          typename T17 = TypeListEnd ,
          typename T18 = TypeListEnd ,
          typename T19 = TypeListEnd ,
          typename T20 = TypeListEnd ,
          typename T21 = TypeListEnd ,
          typename T22 = TypeListEnd ,
          typename T23 = TypeListEnd ,
          typename T24 = TypeListEnd ,
          typename T25 = TypeListEnd ,
          typename T26 = TypeListEnd ,
          typename T27 = TypeListEnd ,
          typename T28 = TypeListEnd ,
          typename T29 = TypeListEnd ,
          typename T30 = TypeListEnd ,
          typename T31 = TypeListEnd ,
          typename T32 = TypeListEnd ,
          typename T33 = TypeListEnd ,
          typename T34 = TypeListEnd ,
          typename T35 = TypeListEnd ,
          typename T36 = TypeListEnd ,
          typename T37 = TypeListEnd ,
          typename T38 = TypeListEnd ,
          typename T39 = TypeListEnd ,
          typename T40 = TypeListEnd ,
          typename T41 = TypeListEnd ,
          typename T42 = TypeListEnd ,
          typename T43 = TypeListEnd ,
          typename T44 = TypeListEnd ,
          typename T45 = TypeListEnd ,
          typename T46 = TypeListEnd ,
          typename T47 = TypeListEnd ,
          typename T48 = TypeListEnd ,
          typename T49 = TypeListEnd ,
          typename T50 = TypeListEnd ,
          typename T51 = TypeListEnd ,
          typename T52 = TypeListEnd ,
          typename T53 = TypeListEnd ,
          typename T54 = TypeListEnd ,
          typename T55 = TypeListEnd ,
          typename T56 = TypeListEnd ,
          typename T57 = TypeListEnd ,
          typename T58 = TypeListEnd ,
          typename T59 = TypeListEnd ,
          typename T60 = TypeListEnd ,
          typename T61 = TypeListEnd ,
          typename T62 = TypeListEnd ,
          typename T63 = TypeListEnd >
struct MakeTypeList
{
#ifndef DOXYGEN_COMPILE
private:
  typedef  TypeList< T00 ,
           TypeList< T01 ,
           TypeList< T02 ,
           TypeList< T03 ,
           TypeList< T04 ,
           TypeList< T05 ,
           TypeList< T06 ,
           TypeList< T07 ,
           TypeList< T08 ,
           TypeList< T09 ,
           TypeList< T10 ,
           TypeList< T11 ,
           TypeList< T12 ,
           TypeList< T13 ,
           TypeList< T14 ,
           TypeList< T15 ,
           TypeList< T16 ,
           TypeList< T17 ,
           TypeList< T18 ,
           TypeList< T19 ,
           TypeList< T20 ,
           TypeList< T21 ,
           TypeList< T22 ,
           TypeList< T23 ,
           TypeList< T24 ,
           TypeList< T25 ,
           TypeList< T26 ,
           TypeList< T27 ,
           TypeList< T28 ,
           TypeList< T29 ,
           TypeList< T30 ,
           TypeList< T31 ,
           TypeList< T32 ,
           TypeList< T33 ,
           TypeList< T34 ,
           TypeList< T35 ,
           TypeList< T36 ,
           TypeList< T37 ,
           TypeList< T38 ,
           TypeList< T39 ,
           TypeList< T40 ,
           TypeList< T41 ,
           TypeList< T42 ,
           TypeList< T43 ,
           TypeList< T44 ,
           TypeList< T45 ,
           TypeList< T46 ,
           TypeList< T47 ,
           TypeList< T48 ,
           TypeList< T49 ,
           TypeList< T50 ,
           TypeList< T51 ,
           TypeList< T52 ,
           TypeList< T53 ,
           TypeList< T54 ,
           TypeList< T55 ,
           TypeList< T56 ,
           TypeList< T57 ,
           TypeList< T58 ,
           TypeList< T59 ,
           TypeList< T60 ,
           TypeList< T61 ,
           TypeList< T62 ,
           TypeList< T63 ,
           TypeListEnd > > > > > > > > > > > > > > > >
                       > > > > > > > > > > > > > > > >
                       > > > > > > > > > > > > > > > >
                       > > > > > > > > > > > > > > > > dirty_type ;
#endif /* DOXYGEN_COMPILE */
public:

  /** \brief The constructed type list. */
  typedef typename TypeListClean< dirty_type >::type type ;

  /** \brief Length of the constructed type list. */
  enum { length = TypeListLength<type>::value };

  /** \brief If every member of the constructed type list is unique. */
  enum { unique = TypeListUnique<type>::value };
};

/** \} */
} // namespace shards

#endif // Shards_TypeList_hpp

