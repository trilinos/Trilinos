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
template< unsigned  I0 = 0 , unsigned  I1 = 0 ,
          unsigned  I2 = 0 , unsigned  I3 = 0 ,
          unsigned  I4 = 0 , unsigned  I5 = 0 ,
          unsigned  I6 = 0 , unsigned  I7 = 0 ,
          unsigned  I8 = 0 , unsigned  I9 = 0 ,
          unsigned I10 = 0 , unsigned I11 = 0 ,
          unsigned I12 = 0 , unsigned I13 = 0 ,
          unsigned I14 = 0 , unsigned I15 = 0 ,
          unsigned I16 = 0 , unsigned I17 = 0 ,
          unsigned I18 = 0 , unsigned I19 = 0 ,
          unsigned I20 = 0 , unsigned I21 = 0 ,
          unsigned I22 = 0 , unsigned I23 = 0 ,
          unsigned I24 = 0 , unsigned I25 = 0 ,
          unsigned I26 = 0 , unsigned I27 = 0 ,
          unsigned I28 = 0 , unsigned I29 = 0 ,
          unsigned I30 = 0 , unsigned I31 = 0 >
struct IndexList {};

/** \brief Access member of compile-time list of indices. <br>
 *         Defines <b> enum { value = Jth index in the List }; </b> 
 */
template< class List , unsigned J > struct IndexListAt {};

#ifndef DOXYGEN_COMPILE

#define SHARDS_INDEX_LIST_AT_SPECIALIZATION( J , K )	\
  template< unsigned  I0 , unsigned  I1 ,	\
            unsigned  I2 , unsigned  I3 ,	\
            unsigned  I4 , unsigned  I5 ,	\
            unsigned  I6 , unsigned  I7 ,	\
            unsigned  I8 , unsigned  I9 ,	\
            unsigned I10 , unsigned I11 ,	\
            unsigned I12 , unsigned I13 ,	\
            unsigned I14 , unsigned I15 ,	\
            unsigned I16 , unsigned I17 ,	\
            unsigned I18 , unsigned I19 ,	\
            unsigned I20 , unsigned I21 ,	\
            unsigned I22 , unsigned I23 ,	\
            unsigned I24 , unsigned I25 ,	\
            unsigned I26 , unsigned I27 ,	\
            unsigned I28 , unsigned I29 ,	\
            unsigned I30 , unsigned I31 >	\
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

#endif /* DOXYGEN_COMPILE */

/** \} */

} // namespace shards


#endif // Shards_IndexList_hpp

