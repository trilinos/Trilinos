/*
//@HEADER
// ************************************************************************
//
//                             KokkosArray
//         Manycore Performance-Portable Multidimensional Arrays
//
//              Copyright (2012) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions?  Contact  H. Carter Edwards (hcedwar@sandia.gov)
// 
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOSARRAY_CRSARRAY_HPP
#define KOKKOSARRAY_CRSARRAY_HPP

#include <string>
#include <vector>

#include <KokkosArray_View.hpp>

namespace KokkosArray {

//----------------------------------------------------------------------------
/** \brief  Compressed row storage array.
 *
 *  A row has a range of entries:
 *
 *    row_map[i0] <= entry < row_map[i0+1]
 *    0 <= i1 < row_map[i0+1] - row_map[i0]
 *
 *  entries( entry ,            i2 , i3 , ... );
 *  entries( row_map[i0] + i1 , i2 , i3 , ... );
 */

template< class DataType ,
          class Arg1Type ,
          class Arg2Type = void ,
          typename SizeType = typename ViewTraits< DataType* , Arg1Type , Arg2Type , void >::size_type >
class CrsArray {
private:

  typedef ViewTraits< DataType* , Arg1Type , Arg2Type , void > traits ;

public:

  typedef DataType                                            data_type ;
  typedef typename traits::array_layout                       array_layout ;
  typedef typename traits::device_type                        device_type ;
  typedef SizeType                                            size_type ;

  typedef CrsArray< DataType , Arg1Type , Arg2Type , SizeType > crsarray_type ;
  typedef CrsArray< DataType , array_layout , Host , SizeType > HostMirror ;
  typedef View< const size_type* , array_layout, device_type >  row_map_type ;
  typedef View<       DataType*  , array_layout, device_type >  entries_type ;

  entries_type entries ;
  row_map_type row_map ;

  /** \brief  Construct a NULL view */
  CrsArray() : entries(), row_map() {}

  /** \brief  Construct a view of the array */
  CrsArray( const CrsArray & rhs )
    : entries( rhs.entries )
    , row_map( rhs.row_map )
    {}

  /** \brief  Assign to a view of the rhs array.
   *          If the old view is the last view
   *          then allocated memory is deallocated.
   */
  CrsArray & operator = ( const CrsArray & rhs )
    {
      entries   = rhs.entries ;
      row_map   = rhs.row_map ;
      return *this ;
    }

  /**  \brief  Destroy this view of the array.
   *           If the last view then allocated memory is deallocated.
   */
  ~CrsArray() {}
};

//----------------------------------------------------------------------------

template< class CrsArrayType , class InputSizeType >
typename CrsArrayType::crsarray_type
create_crsarray( const std::string & label ,
                 const std::vector< InputSizeType > & input );

template< class CrsArrayType , class InputSizeType >
typename CrsArrayType::crsarray_type
create_crsarray( const std::string & label ,
                 const std::vector< std::vector< InputSizeType > > & input );

//----------------------------------------------------------------------------

template< class DataType ,
          class Arg1Type ,
          class Arg2Type ,
          typename SizeType >
typename CrsArray< DataType , Arg1Type , Arg2Type , SizeType >::HostMirror
create_mirror_view( const CrsArray<DataType,Arg1Type,Arg2Type,SizeType > & input );

template< class DataType ,
          class Arg1Type ,
          class Arg2Type ,
          typename SizeType >
typename CrsArray< DataType , Arg1Type , Arg2Type , SizeType >::HostMirror
create_mirror( const CrsArray<DataType,Arg1Type,Arg2Type,SizeType > & input );

} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#include <impl/KokkosArray_CrsArray_factory.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOSARRAY_CRSARRAY_HPP */

