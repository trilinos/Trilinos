/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#ifndef stk_util_util_TypeList_h
#define stk_util_util_TypeList_h

#include <stk_util/util/TypeUtil.hpp>
#include <stk_util/util/SameType.hpp>

namespace stk {

/** \ingroup  util_module
 *  \defgroup  typelist_module  TypeList: linked list of types
 *  \brief     Linked-list of compile-time types and
 *             supporting compile-time linked list operations.
 *  \author H. Carter Edwards  <hcedwar@sandia.gov>
 *
 *  'TypeList' templates significantly enhanced from
 *  Alexandrescu's "Modern C++ Design" book.
 */

//----------------------------------------------------------------------

struct TypeListEnd {};

/** \class  TypeList
 *  \brief  A link within a linked list of types.
 *  \ingroup typelist_module
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
 *  \ingroup typelist_module
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
 *  \ingroup typelist_module
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
 *  \ingroup typelist_module
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
  enum { match = stk::SameType< Value , TestValue >::value };
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
 *  \ingroup typelist_module
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
                 ( stk::SameType< Value , TestValue >::value ? 1 : 0 ) };
};

//----------------------------------------------------------------------
/** \class  TypeListMember
 *  \brief  Member <b> enum { value = ... }; </b> is true
 *          if <b> TestValue </b> is a member of <b> ListType </b>.
 *  \ingroup typelist_module
 */
template< class ListType , typename TestValue > struct TypeListMember {};

template< typename TestValue >
struct TypeListMember< TypeListEnd , TestValue >
{ enum { value = false }; };

template< typename Value , class Tail , typename TestValue >
struct TypeListMember< TypeList< Value , Tail > , TestValue >
{
  enum { value = stk::SameType< Value , TestValue >::value ||
                 TypeListMember< Tail , TestValue >::value };
};

//----------------------------------------------------------------------
/** \class  TypeListUnique
 *  \brief  Member <b> enum { value = ... }; </b> is true
 *          if each member of <b> ListType </b> appears exactly once.
 *  \ingroup typelist_module
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
 *  \ingroup typelist_module
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
 *  \ingroup typelist_module
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
 *  \ingroup typelist_module
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
 *  \ingroup typelist_module
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
 *  \ingroup typelist_module
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
 *  \ingroup typelist_module
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
 *  \ingroup typelist_module
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

} // namespace stk

namespace sierra {

//----------------------------------------------------------------------

/** @class  sierra::TypeList
 *  @brief  Linked list of types.
 */
template<typename ValueType, typename ListType>
struct TypeList {
  /** Value for the current entry */
  typedef ValueType TypeListValue ;

  /** Remainder of the list */
  typedef ListType  TypeListTail ;
};

/** The end of a TypeList is the first encounter with an entry of TypeListEnd.
 *  Thus all entries after a TypeListEnd entry are ignored.  This allows
 *  The MakeTypeList<...> template to define TypeLists of the desired length.
 */
struct TypeListEnd {};

//----------------------------------------------------------------------

/** Length of a TypeList */
template< class ListType>
struct TypeListLength /* { enum { value = <> }; } */ ;

/** Location of a ValueType in the TypeList */
template< class ListType, typename ValueType, unsigned Ordinal = 0>
struct TypeListIndex /* { enum { value = <> }; } */ ;

/** Count of appearances of ValueType in the TypeList */
template< class ListType, typename ValueType>
struct TypeListCount /* { enum { value = <> }; } */ ;

/** Last ValueType in the TypeList */
template< class ListType >
struct TypeListLast /* { typedef <> type ; } */ ;

/** ValueType and sub-ListType at a location in a TypeList */
template< class ListType, unsigned I>
struct TypeListAt /* { typedef <> type ; typedef <> list_type ; } */ ;

/** TypeList member of a ValueType in a TypeList */
template< class ListType, typename ValueType>
struct TypeListMember /* { typedef <> list_type ; } */ ;

/** Erase type from TypeList at I */
template< class ListType, unsigned I >
struct TypeListEraseAt /* { typedef <> list_type ; } */ ;

// /** Erase type from TypeList at I */
// template< class ListType, int I >
// struct TypeListErase /* { typedef <> list_type ; } */ ;

/** Check for uniqueness of a TypeList */
template< class ListType >
struct TypeListUnique /* { enum { value = <> }; } */ ;

/** Check if SuperList contains SubList */
template< class SuperList , class SubList >
struct TypeListContains /* { enum { value = <> }; } */ ;

/** Check if ListA is disjoint from ListB */
template< class ListA , class ListB >
struct TypeListDisjoint /* { enum { value = <> }; } */ ;

/** Truncate a TypeList at the first appearance of TypeListEnd */
template<class ListType>
struct TypeListClean /* { typedef <> list_type ; } */ ;

/** Make a TypeList from a sequence of type entries.
 *  Implemented to support list of up to thirtytwo (32) types.
 */
template< typename T0 = TypeListEnd ,
	  typename T1 = TypeListEnd ,
	  typename T2 = TypeListEnd ,
	  typename T3 = TypeListEnd ,
	  typename T4 = TypeListEnd ,
	  typename T5 = TypeListEnd ,
	  typename T6 = TypeListEnd ,
	  typename T7 = TypeListEnd ,
	  typename T8 = TypeListEnd ,
	  typename T9 = TypeListEnd ,
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
	  typename T31 = TypeListEnd >
struct MakeTypeList
/* { typedef <> type ; enum { length = <> , unique = <> }; } */ ;

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// Clean: End -> End
// Clean: <End,  Tail> -> End
// Clean: <Value,Tail> -> < Value , Clean<Tail> >

template<>
struct TypeListClean<TypeListEnd> {
  typedef TypeListEnd list_type ;
};

template<typename Tail>
struct TypeListClean< TypeList< TypeListEnd , Tail > > {
  typedef TypeListEnd list_type ;
};

template<class ListType>
struct TypeListClean {
private:
  typedef typename ListType::TypeListValue                 ValueType ;
  typedef typename ListType::TypeListTail                  InputTailType ;
  typedef typename TypeListClean<InputTailType>::list_type TailType ;
public:
  typedef TypeList< ValueType , TailType > list_type ;
};

//----------
// Length: End          -> 0
// Length: <End,  Tail> -> 0
// Length: <Value,Tail> -> 1 + Length<Tail>

template<> struct TypeListLength<TypeListEnd>
{ enum { value = 0 }; };

template<class ListType>
struct TypeListLength< TypeList<TypeListEnd,ListType> >
{ enum { value = 0 }; };

template<class ListType>
struct TypeListLength
{
private: typedef typename ListType::TypeListTail TailType ;
public:  enum { value = 1 + TypeListLength<TailType>::value };
};

//----------
// Index:  < End ,                  ValueType >  ->  -1
// Index:  < List ,                 End >        ->  -1
// Index:  < < End ,       Tail > , ValueType >  ->  -1
// Index:  < < ValueType , Tail > , ValueType >  ->   0
// Index:  < < OtherType , Tail > , ValueType >  ->
//         ( I = Index<List::Tail,ValueType> , I == -1 ? -1 : I + 1 )

template<typename ValueType, unsigned Ordinal>
struct TypeListIndex< TypeListEnd , ValueType , Ordinal> {
  enum { value = -1 };
  typedef TypeListEnd tail_type ;
};

template<class ListType, unsigned Ordinal>
struct TypeListIndex< ListType , TypeListEnd , Ordinal > {
  enum { value = -1 };
  typedef TypeListEnd tail_type ;
};

template<class Tail, typename ValueType, unsigned Ordinal >
struct TypeListIndex< TypeList<TypeListEnd,Tail> , ValueType, Ordinal >
{
  enum { value = -1 };
  typedef TypeListEnd tail_type ;
};

template<typename ValueType, class Tail>
struct TypeListIndex< TypeList<ValueType,Tail> , ValueType , 0 >
{
  enum { value = 0 };
  typedef Tail tail_type ;
};

// Condition: ValueType == ListType::TypeListValue && Ordinal matches

template<class ListType, typename ValueType, unsigned Ordinal>
struct TypeListIndex
{
private:
  enum { same = stk::SameType<ValueType,typename ListType::TypeListValue>::value };
  enum { ord = Ordinal == 0 ? 0 : ( same ? Ordinal - 1 : Ordinal ) };
  typedef typename ListType::TypeListTail TailType ;
  typedef TypeListIndex< TailType , ValueType , ord > type_list_index ;
  enum { temp = type_list_index::value };
public:
  enum { value = temp == -1 ? -1 : 1 + temp };
  typedef typename type_list_index::tail_type tail_type ;
};

//----------
// Count :  < End , ValueType >  ->  0
// Count :  < List , End >       ->  0
// Count :  < < End ,       Tail > , ValueType >  ->  0
// Count :  < < ValueType , Tail > , ValueType >  ->  Count<Tail,ValueType> + 1
// Count :  < < OtherType , Tail > , ValueType >  ->  Count<Tail,ValueType>

template<typename ValueType>
struct TypeListCount< TypeListEnd , ValueType > { enum { value = 0 }; };

template<class ListType>
struct TypeListCount< ListType , TypeListEnd > { enum { value = 0 }; };

template<class Tail, typename ValueType>
struct TypeListCount< TypeList<TypeListEnd,Tail>,ValueType>
{ enum { value = 0 }; };

template<typename ValueType, class Tail>
struct TypeListCount< TypeList<ValueType,Tail> , ValueType>
{ enum { value = 1 + TypeListCount< Tail , ValueType >::value }; };

template<class ListType, typename ValueType>
struct TypeListCount
{
private: typedef typename ListType::TypeListTail TailType ;
public:  enum { value = TypeListCount< TailType , ValueType >::value };
};

//----------
// At :  < End ,                  0 >  ->  { End , End }
// At :  < End ,                  I >  ->  { End , End }
// At :  < < End ,       Tail > , I >  ->  { End , End }
// At :  < < ValueType , Tail > , 0 >  ->  { ValueType , < ValueType , Tail > }
// At :  < < ValueType , Tail > , I >  ->  At< Tail , I - 1 >

template<>
struct TypeListAt< TypeListEnd, 0>
{
  typedef TypeListEnd type ;
  typedef TypeListEnd list_type ;
};

template<unsigned I>
struct TypeListAt< TypeListEnd, I>
{
  typedef TypeListEnd type ;
  typedef TypeListEnd list_type ;
};

template< class ListType >
struct TypeListAt< ListType , 0 >
{
private:
  typedef typename ListType::TypeListTail Tail ;
public:
  typedef typename ListType::TypeListValue type ;
  typedef TypeList< type , Tail >          list_type ;
};

template<class Tail, unsigned I>
struct TypeListAt< TypeList<TypeListEnd,Tail>, I>
{
  typedef TypeListEnd type ;
  typedef TypeListEnd list_type ;
};

template<class ListType, unsigned I>
struct TypeListAt
{
private:
  typedef typename ListType::TypeListTail Tail ;
  typedef TypeListAt<Tail,I-1>            AtType ;
public:
  typedef typename AtType::type      type ;
  typedef typename AtType::list_type list_type ;
};

//----------
// Last : End -> End
// Last : < ValueType , End >             ->  ValueType
// Last : < ValueType , < End , Tail > >  ->  ValueType
// Last : < ValueType , Tail >            ->  Last< Tail >

template<>
struct TypeListLast< TypeListEnd >
{ typedef TypeListEnd type ; };

template<class ValueType>
struct TypeListLast< TypeList<ValueType,TypeListEnd> >
{ typedef ValueType type ; };

template<class ValueType,class Tail>
struct TypeListLast< TypeList<ValueType,TypeList<TypeListEnd,Tail> > >
{ typedef ValueType type ; };

template<class ValueType, class Tail>
struct TypeListLast< TypeList<ValueType,Tail> >
{ typedef typename TypeListLast<Tail>::type type ; };

//----------
// Member :
// Member :
// Member :
//

template< typename ValueType >
struct TypeListMember< TypeListEnd , ValueType >
{ typedef TypeListEnd list_type ; };

template< class Tail , typename ValueType >
struct TypeListMember< TypeList<TypeListEnd,Tail> , ValueType >
{ typedef TypeListEnd list_type ; };

template< typename ValueType , class ListType>
struct TypeListMember< TypeList<ValueType,ListType> , ValueType >
{ typedef TypeList<ValueType,ListType> list_type ; };

template< class ListType, typename ValueType>
struct TypeListMember
{
  private: typedef typename ListType::TypeListTail Tail ;
  public:  typedef typename TypeListMember<Tail,ValueType>::list_type list_type;
};

// //----------
// // Erase :
// // Erase :
// // Erase :
// //

// template< class ValueType>
// struct TypeListErase<TypeListEnd, ValueType>
// { typedef TypeListEnd list_type ; };

// template< class ValueType, class Tail>
// struct TypeListErase< TypeList<ValueType, Tail> , ValueType >
// { typedef Tail list_type ; };

// template< class ListType, class ValueType>
// struct TypeListErase< TypeList<Head, Tail>, ValueType>
// {
// private: typedef typename ListType::TypeListTail Tail ;
// public:  typedef TypeList<Head, typename TypeListErase<Tail, ValueType>::list_type> list_type;
// };

//----------
// EraseAt :  < End ,                  0 >  ->  { End , End }
// EraseAt :  < End ,                  I >  ->  { End , End }
// EraseAt :  < < End ,       Tail > , I >  ->  { End , End }
// EraseAt :  < < ListType ,           0 >  ->  { TypeList < ValueType , Tail > }
// EraseAt :  < < ListType , Tail > ,  I >  ->  { EraseAt< Tail , I - 1 > }

template<>
struct TypeListEraseAt< TypeListEnd, 0>
{
  typedef TypeListEnd				list_type ;
};

template<unsigned I>
struct TypeListEraseAt< TypeListEnd, I>
{
  typedef TypeListEnd				list_type ;
};

template<class Tail, unsigned I>
struct TypeListEraseAt< TypeList<TypeListEnd,Tail>, I>
{
  typedef TypeListEnd				list_type ;
};

template< class ListType >
struct TypeListEraseAt< ListType , 0 >
{
private:
  typedef typename ListType::TypeListTail	Tail ;
public:
  typedef Tail					list_type ;
};

template<class ListType, unsigned I>
struct TypeListEraseAt
{
private:
  typedef typename ListType::TypeListTail	Tail ;
  typedef TypeListEraseAt<Tail, I - 1>		EraseAtType ;
public:
  typedef TypeList<typename ListType::TypeListValue, typename EraseAtType::list_type>	list_type ;
};

//----------
// Unique :  End  ->  true
// Unique :  < End , Tail >  -> true
// Unique :  < ValueType , Tail >  ->
//           Index<Tail,ValueType> == -1 && Unique<Tail>

template<>
struct TypeListUnique<TypeListEnd> { enum { value = true }; };

template<class Tail>
struct TypeListUnique< TypeList<TypeListEnd,Tail> >
{ enum { value = true }; };

template< class ListType >
struct TypeListUnique
{
private:
  typedef typename ListType::TypeListValue ValueType ;
  typedef typename ListType::TypeListTail  TailType ;
public:
  // This ValueType does not appear in the remainder of the TypeList and
  // the remainder of the TypeList is also unique.
  enum { value = ( TypeListIndex<TailType,ValueType>::value == -1 ) &&
		   TypeListUnique<TailType>::value };
};

//----------
// Contains : < SuperList , End >  -> true
// Contains : < SuperList , < End , Tail > >  -> true
// Contains : < SuperList , SubList >  ->
//            Index<   SuperList,SubList::Value> != -1 &&
//            Contains<SuperList,SubList::Tail>

template<class SuperList>
struct TypeListContains<SuperList,TypeListEnd>
{ enum { value = true }; };

template<class SuperList,typename Tail>
struct TypeListContains<SuperList,TypeList<TypeListEnd,Tail> >
{ enum { value = true }; };

template<class SuperList, class SubList >
struct TypeListContains
{
private:
  typedef typename SubList::TypeListValue ValueType ;
  typedef typename SubList::TypeListTail  TailType ;
public:
  // The SuperList contains this ValueType and the remainder of the SubList
  enum { value = ( TypeListIndex<SuperList,ValueType>::value != -1 ) &&
		   TypeListContains<SuperList,TailType>::value };
};

//----------
// Disjoint : < ListA , End >  ->  true
// Disjoint : < ListA , < End , Tail > >  ->  true
// Disjoint : < ListA , ListB >  ->
//            Index<   ListA,ListB::Value> == -1 &&
//            Disjoint<ListA,ListB::Tail>

template<class SuperList>
struct TypeListDisjoint<SuperList,TypeListEnd>
{ enum { value = true }; };

template<class SuperList,typename Tail>
struct TypeListDisjoint<SuperList,TypeList<TypeListEnd,Tail> >
{ enum { value = true }; };

template<class ListA, class ListB>
struct TypeListDisjoint
{
private:
  typedef typename ListB::TypeListValue ValueType ;
  typedef typename ListB::TypeListTail  TailType ;
public:
  // ListA does not contain this ValueType and does not contain the remainder
  enum { value = ( TypeListIndex<ListA,ValueType>::value == -1 ) &&
		   TypeListDisjoint<ListA,TailType>::value };
};

//----------------------------------------------------------------------

template< typename T1 ,
	  typename T2 ,
	  typename T3 ,
	  typename T4 ,
	  typename T5 ,
	  typename T6 ,
	  typename T7 ,
	  typename T8 ,
	  typename T9 ,
	  typename T10 ,
	  typename T11 ,
	  typename T12 ,
	  typename T13 ,
	  typename T14 ,
	  typename T15 ,
	  typename T16 ,
	  typename T17 ,
	  typename T18 ,
	  typename T19 ,
	  typename T20 ,
	  typename T21 ,
	  typename T22 ,
	  typename T23 ,
	  typename T24 ,
	  typename T25 ,
	  typename T26 ,
	  typename T27 ,
	  typename T28 ,
	  typename T29 ,
	  typename T30 ,
	  typename T31 >
struct MakeTypeList<TypeListEnd,
		    T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,
		    T11,T12,T13,T14,T15,T16,T17,T18,T19,T20,
		    T21,T22,T23,T24,T25,T26,T27,T28,T29,T30,T31> {
  typedef TypeListEnd type ;
  enum { length = 0 };
  enum { unique = true };
};

template< typename T0 ,
	  typename T1 ,
	  typename T2 ,
	  typename T3 ,
	  typename T4 ,
	  typename T5 ,
	  typename T6 ,
	  typename T7 ,
	  typename T8 ,
	  typename T9 ,
	  typename T10 ,
	  typename T11 ,
	  typename T12 ,
	  typename T13 ,
	  typename T14 ,
	  typename T15 ,
	  typename T16 ,
	  typename T17 ,
	  typename T18 ,
	  typename T19 ,
	  typename T20 ,
	  typename T21 ,
	  typename T22 ,
	  typename T23 ,
	  typename T24 ,
	  typename T25 ,
	  typename T26 ,
	  typename T27 ,
	  typename T28 ,
	  typename T29 ,
	  typename T30 ,
	  typename T31 >
struct MakeTypeList {
  typedef typename TypeListClean<
	  TypeList< T0 ,
	  TypeList< T1,
	  TypeList< T2 ,
	  TypeList< T3 ,
	  TypeList< T4 ,
	  TypeList< T5 ,
	  TypeList< T6 ,
	  TypeList< T7 ,
	  TypeList< T8 ,
	  TypeList< T9 ,
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
	  TypeListEnd > > > > > > > > > > > > > > > >
		      > > > > > > > > > > > > > > > >
  >::list_type type ;
  enum { length = TypeListLength<type>::value };
  enum { unique = TypeListUnique<type>::value };
};

} // namespace sierra

#endif // stk_util_util_TypeList_h

