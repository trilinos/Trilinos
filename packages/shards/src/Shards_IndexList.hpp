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

#ifndef Shards_IndexList_hpp
#define Shards_IndexList_hpp

namespace shards {

/** \ingroup shards_package
 *  \defgroup  shards_package_index_list  Compile-time List of Indices
 *  \brief  Compile-time list of indices and access templates.
 *
 *  Compile-time list of indices and access operations.
 *
 *  \author H. Carter Edwards  <hcedwar@sandia.gov>
 *
 *  \{
 */

/** \brief Compile-time list of indices.  */
template< int  I0 = -1 , int  I1 = -1 , int  I2 = -1 , int  I3 = -1 ,
          int  I4 = -1 , int  I5 = -1 , int  I6 = -1 , int  I7 = -1 ,
          int  I8 = -1 , int  I9 = -1 , int I10 = -1 , int I11 = -1 ,
          int I12 = -1 , int I13 = -1 , int I14 = -1 , int I15 = -1 ,
          int I16 = -1 , int I17 = -1 , int I18 = -1 , int I19 = -1 ,
          int I20 = -1 , int I21 = -1 , int I22 = -1 , int I23 = -1 ,
          int I24 = -1 , int I25 = -1 , int I26 = -1 , int I27 = -1 ,
          int I28 = -1 , int I29 = -1 , int I30 = -1 , int I31 = -1 >
struct IndexList {};

/** \brief  Access member of compile-time list of indices. <br>
 *          Defines <b> enum { value = Jth member }; </b> 
 */
template< class List , unsigned J > struct IndexListAt {};

/** \brief  Find member of compile-time list of indices. <br>
 *          Defines <b> enum { value = index of member equal to J }; </b> 
 */
template< class List , unsigned J > struct IndexListFind ;

/** \brief  Length of list.
 *          Defines <b> enum { value }; </b>
 */
template< class List > struct IndexListLength {};

/** \brief  Inverse of list containing [0..N].
 *          Defines <b> typedef IndexList<...> type ; </b>
 */
template< class List > struct IndexListInverse {};

#ifndef DOXYGEN_COMPILE

//----------------------------------------------------------------------

#define SHARDS_INDEX_LIST_AT_SPECIALIZATION( J , K )	\
  template< int  I0 , int  I1 ,	int  I2 , int  I3 ,	\
            int  I4 , int  I5 ,	int  I6 , int  I7 ,	\
            int  I8 , int  I9 ,	int I10 , int I11 ,	\
            int I12 , int I13 ,	int I14 , int I15 ,	\
            int I16 , int I17 ,	int I18 , int I19 ,	\
            int I20 , int I21 ,	int I22 , int I23 ,	\
            int I24 , int I25 ,	int I26 , int I27 ,	\
            int I28 , int I29 ,	int I30 , int I31 >	\
struct IndexListAt<	\
  IndexList< I0 ,  I1 ,  I2 ,  I3 ,  I4 ,  I5 ,  I6 ,  I7 ,	\
             I8 ,  I9 , I10 , I11 , I12 , I13 , I14 , I15 ,	\
            I16 , I17 , I18 , I19 , I20 , I21 , I22 , I23 ,	\
            I24 , I25 , I26 , I27 , I28 , I29 , I30 , I31 > , J >	\
{ enum { value = K }; };

SHARDS_INDEX_LIST_AT_SPECIALIZATION(  0 ,  I0 )
SHARDS_INDEX_LIST_AT_SPECIALIZATION(  1 ,  I1 )
SHARDS_INDEX_LIST_AT_SPECIALIZATION(  2 ,  I2 )
SHARDS_INDEX_LIST_AT_SPECIALIZATION(  3 ,  I3 )
SHARDS_INDEX_LIST_AT_SPECIALIZATION(  4 ,  I4 )
SHARDS_INDEX_LIST_AT_SPECIALIZATION(  5 ,  I5 )
SHARDS_INDEX_LIST_AT_SPECIALIZATION(  6 ,  I6 )
SHARDS_INDEX_LIST_AT_SPECIALIZATION(  7 ,  I7 )
SHARDS_INDEX_LIST_AT_SPECIALIZATION(  8 ,  I8 )
SHARDS_INDEX_LIST_AT_SPECIALIZATION(  9 ,  I9 )
SHARDS_INDEX_LIST_AT_SPECIALIZATION( 10 , I10 )
SHARDS_INDEX_LIST_AT_SPECIALIZATION( 11 , I11 )
SHARDS_INDEX_LIST_AT_SPECIALIZATION( 12 , I12 )
SHARDS_INDEX_LIST_AT_SPECIALIZATION( 13 , I13 )
SHARDS_INDEX_LIST_AT_SPECIALIZATION( 14 , I14 )
SHARDS_INDEX_LIST_AT_SPECIALIZATION( 15 , I15 )
SHARDS_INDEX_LIST_AT_SPECIALIZATION( 16 , I16 )
SHARDS_INDEX_LIST_AT_SPECIALIZATION( 17 , I17 )
SHARDS_INDEX_LIST_AT_SPECIALIZATION( 18 , I18 )
SHARDS_INDEX_LIST_AT_SPECIALIZATION( 19 , I19 )
SHARDS_INDEX_LIST_AT_SPECIALIZATION( 20 , I20 )
SHARDS_INDEX_LIST_AT_SPECIALIZATION( 21 , I21 )
SHARDS_INDEX_LIST_AT_SPECIALIZATION( 22 , I22 )
SHARDS_INDEX_LIST_AT_SPECIALIZATION( 23 , I23 )
SHARDS_INDEX_LIST_AT_SPECIALIZATION( 24 , I24 )
SHARDS_INDEX_LIST_AT_SPECIALIZATION( 25 , I25 )
SHARDS_INDEX_LIST_AT_SPECIALIZATION( 26 , I26 )
SHARDS_INDEX_LIST_AT_SPECIALIZATION( 27 , I27 )
SHARDS_INDEX_LIST_AT_SPECIALIZATION( 28 , I28 )
SHARDS_INDEX_LIST_AT_SPECIALIZATION( 29 , I29 )
SHARDS_INDEX_LIST_AT_SPECIALIZATION( 30 , I30 )
SHARDS_INDEX_LIST_AT_SPECIALIZATION( 31 , I31 )

#undef SHARDS_INDEX_LIST_AT_SPECIALIZATION

//----------------------------------------------------------------------

template< class List , unsigned J >
struct IndexListFind { enum { value = -1 }; };

#define SHARDS_INDEX_LIST_FIND_SPECIALIZATION( J , K )	\
  template< int  I0 , int  I1 ,	int  I2 , int  I3 ,	\
            int  I4 , int  I5 ,	int  I6 , int  I7 ,	\
            int  I8 , int  I9 ,	int I10 , int I11 ,	\
            int I12 , int I13 ,	int I14 , int I15 ,	\
            int I16 , int I17 ,	int I18 , int I19 ,	\
            int I20 , int I21 ,	int I22 , int I23 ,	\
            int I24 , int I25 ,	int I26 , int I27 ,	\
            int I28 , int I29 ,	int I30 , int I31 >	\
struct IndexListFind<	\
  IndexList< I0 ,  I1 ,  I2 ,  I3 ,  I4 ,  I5 ,  I6 ,  I7 ,	\
             I8 ,  I9 , I10 , I11 , I12 , I13 , I14 , I15 ,	\
            I16 , I17 , I18 , I19 , I20 , I21 , I22 , I23 ,	\
            I24 , I25 , I26 , I27 , I28 , I29 , I30 , I31 > , K >	\
{ enum { value = J }; };

SHARDS_INDEX_LIST_FIND_SPECIALIZATION(  0 ,  I0 )
SHARDS_INDEX_LIST_FIND_SPECIALIZATION(  1 ,  I1 )
SHARDS_INDEX_LIST_FIND_SPECIALIZATION(  2 ,  I2 )
SHARDS_INDEX_LIST_FIND_SPECIALIZATION(  3 ,  I3 )
SHARDS_INDEX_LIST_FIND_SPECIALIZATION(  4 ,  I4 )
SHARDS_INDEX_LIST_FIND_SPECIALIZATION(  5 ,  I5 )
SHARDS_INDEX_LIST_FIND_SPECIALIZATION(  6 ,  I6 )
SHARDS_INDEX_LIST_FIND_SPECIALIZATION(  7 ,  I7 )
SHARDS_INDEX_LIST_FIND_SPECIALIZATION(  8 ,  I8 )
SHARDS_INDEX_LIST_FIND_SPECIALIZATION(  9 ,  I9 )
SHARDS_INDEX_LIST_FIND_SPECIALIZATION( 10 , I10 )
SHARDS_INDEX_LIST_FIND_SPECIALIZATION( 11 , I11 )
SHARDS_INDEX_LIST_FIND_SPECIALIZATION( 12 , I12 )
SHARDS_INDEX_LIST_FIND_SPECIALIZATION( 13 , I13 )
SHARDS_INDEX_LIST_FIND_SPECIALIZATION( 14 , I14 )
SHARDS_INDEX_LIST_FIND_SPECIALIZATION( 15 , I15 )
SHARDS_INDEX_LIST_FIND_SPECIALIZATION( 16 , I16 )
SHARDS_INDEX_LIST_FIND_SPECIALIZATION( 17 , I17 )
SHARDS_INDEX_LIST_FIND_SPECIALIZATION( 18 , I18 )
SHARDS_INDEX_LIST_FIND_SPECIALIZATION( 19 , I19 )
SHARDS_INDEX_LIST_FIND_SPECIALIZATION( 20 , I20 )
SHARDS_INDEX_LIST_FIND_SPECIALIZATION( 21 , I21 )
SHARDS_INDEX_LIST_FIND_SPECIALIZATION( 22 , I22 )
SHARDS_INDEX_LIST_FIND_SPECIALIZATION( 23 , I23 )
SHARDS_INDEX_LIST_FIND_SPECIALIZATION( 24 , I24 )
SHARDS_INDEX_LIST_FIND_SPECIALIZATION( 25 , I25 )
SHARDS_INDEX_LIST_FIND_SPECIALIZATION( 26 , I26 )
SHARDS_INDEX_LIST_FIND_SPECIALIZATION( 27 , I27 )
SHARDS_INDEX_LIST_FIND_SPECIALIZATION( 28 , I28 )
SHARDS_INDEX_LIST_FIND_SPECIALIZATION( 29 , I29 )
SHARDS_INDEX_LIST_FIND_SPECIALIZATION( 30 , I30 )
SHARDS_INDEX_LIST_FIND_SPECIALIZATION( 31 , I31 )

#undef SHARDS_INDEX_LIST_FIND_SPECIALIZATION

//----------------------------------------------------------------------

template<>
struct IndexListLength< IndexList<> > { enum { value = 0 }; };

template< int  I0 , int  I1 , int  I2 , int  I3 ,
          int  I4 , int  I5 , int  I6 , int  I7 ,
          int  I8 , int  I9 , int I10 , int I11 ,
          int I12 , int I13 , int I14 , int I15 ,
          int I16 , int I17 , int I18 , int I19 ,
          int I20 , int I21 , int I22 , int I23 ,
          int I24 , int I25 , int I26 , int I27 ,
          int I28 , int I29 , int I30 , int I31 >
struct IndexListLength<
  IndexList<  I0 ,  I1 ,  I2 ,  I3 ,  I4 ,  I5 ,  I6 ,  I7 ,
              I8 ,  I9 , I10 , I11 , I12 , I13 , I14 , I15 ,
             I16 , I17 , I18 , I19 , I20 , I21 , I22 , I23 ,
             I24 , I25 , I26 , I27 , I28 , I29 , I30 , I31 > >
{
private:
  typedef IndexList<        I1 ,  I2 ,  I3 ,  I4 ,  I5 ,  I6 ,  I7 ,
                      I8 ,  I9 , I10 , I11 , I12 , I13 , I14 , I15 ,
                     I16 , I17 , I18 , I19 , I20 , I21 , I22 , I23 ,
                     I24 , I25 , I26 , I27 , I28 , I29 , I30 , I31 , -1 >
    shift_type ;

public:
  enum { value = 1 + IndexListLength< shift_type >::value };
};

//----------------------------------------------------------------------

template< int  I0 , int  I1 , int  I2 , int  I3 ,
          int  I4 , int  I5 , int  I6 , int  I7 ,
          int  I8 , int  I9 , int I10 , int I11 ,
          int I12 , int I13 , int I14 , int I15 ,
          int I16 , int I17 , int I18 , int I19 ,
          int I20 , int I21 , int I22 , int I23 ,
          int I24 , int I25 , int I26 , int I27 ,
          int I28 , int I29 , int I30 , int I31 >
struct IndexListInverse<
  IndexList<  I0 ,  I1 ,  I2 ,  I3 ,  I4 ,  I5 ,  I6 ,  I7 ,
              I8 ,  I9 , I10 , I11 , I12 , I13 , I14 , I15 ,
             I16 , I17 , I18 , I19 , I20 , I21 , I22 , I23 ,
             I24 , I25 , I26 , I27 , I28 , I29 , I30 , I31 > >
{
private:
  typedef IndexList<  I0 ,  I1 ,  I2 ,  I3 ,  I4 ,  I5 ,  I6 ,  I7 ,
                      I8 ,  I9 , I10 , I11 , I12 , I13 , I14 , I15 ,
                     I16 , I17 , I18 , I19 , I20 , I21 , I22 , I23 ,
                     I24 , I25 , I26 , I27 , I28 , I29 , I30 , I31 > list ;

  enum { length = IndexListLength< list >::value };

public:

#define TYPE_LIST_SPECIALIZATION_INVERSE( J )	\
  ( J < length ? IndexListFind< list , J >::value : -1 )

  typedef IndexList< TYPE_LIST_SPECIALIZATION_INVERSE(  0 ) ,
                     TYPE_LIST_SPECIALIZATION_INVERSE(  1 ) ,
                     TYPE_LIST_SPECIALIZATION_INVERSE(  2 ) ,
                     TYPE_LIST_SPECIALIZATION_INVERSE(  3 ) ,
                     TYPE_LIST_SPECIALIZATION_INVERSE(  4 ) ,
                     TYPE_LIST_SPECIALIZATION_INVERSE(  5 ) ,
                     TYPE_LIST_SPECIALIZATION_INVERSE(  6 ) ,
                     TYPE_LIST_SPECIALIZATION_INVERSE(  7 ) ,
                     TYPE_LIST_SPECIALIZATION_INVERSE(  8 ) ,
                     TYPE_LIST_SPECIALIZATION_INVERSE(  9 ) ,
                     TYPE_LIST_SPECIALIZATION_INVERSE( 10 ) ,
                     TYPE_LIST_SPECIALIZATION_INVERSE( 11 ) ,
                     TYPE_LIST_SPECIALIZATION_INVERSE( 12 ) ,
                     TYPE_LIST_SPECIALIZATION_INVERSE( 13 ) ,
                     TYPE_LIST_SPECIALIZATION_INVERSE( 14 ) ,
                     TYPE_LIST_SPECIALIZATION_INVERSE( 15 ) ,
                     TYPE_LIST_SPECIALIZATION_INVERSE( 16 ) ,
                     TYPE_LIST_SPECIALIZATION_INVERSE( 17 ) ,
                     TYPE_LIST_SPECIALIZATION_INVERSE( 18 ) ,
                     TYPE_LIST_SPECIALIZATION_INVERSE( 19 ) ,
                     TYPE_LIST_SPECIALIZATION_INVERSE( 20 ) ,
                     TYPE_LIST_SPECIALIZATION_INVERSE( 21 ) ,
                     TYPE_LIST_SPECIALIZATION_INVERSE( 22 ) ,
                     TYPE_LIST_SPECIALIZATION_INVERSE( 23 ) ,
                     TYPE_LIST_SPECIALIZATION_INVERSE( 24 ) ,
                     TYPE_LIST_SPECIALIZATION_INVERSE( 25 ) ,
                     TYPE_LIST_SPECIALIZATION_INVERSE( 26 ) ,
                     TYPE_LIST_SPECIALIZATION_INVERSE( 27 ) ,
                     TYPE_LIST_SPECIALIZATION_INVERSE( 28 ) ,
                     TYPE_LIST_SPECIALIZATION_INVERSE( 29 ) ,
                     TYPE_LIST_SPECIALIZATION_INVERSE( 30 ) ,
                     TYPE_LIST_SPECIALIZATION_INVERSE( 31 ) > type ;

#undef TYPE_LIST_SPECIALIZATION_INVERSE

};

#endif /* DOXYGEN_COMPILE */

/** \} */

} // namespace shards


#endif // Shards_IndexList_hpp

