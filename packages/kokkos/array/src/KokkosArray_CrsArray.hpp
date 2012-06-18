/*
//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_CRSARRAY_HPP
#define KOKKOS_CRSARRAY_HPP

#include <string>

#include <KokkosArray_PrefixSum.hpp>
#include <KokkosArray_Array.hpp>

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

template< class ArrayType , class DeviceType ,
          typename SizeType = typename DeviceType::size_type >
class CrsArray {
public:
  typedef DeviceType                            device_type ;
  typedef PrefixSum< SizeType ,  device_type >  row_map_type ;
  typedef Array<     ArrayType , device_type >  entries_type ;

  row_map_type row_map ;
  entries_type entries ;

  typedef CrsArray< ArrayType ,
                    typename HostMapped< device_type >::type ,
                    SizeType > HostMirror ;

  /** \brief  Construct a NULL view */
  CrsArray() : row_map() , entries() {}

  /** \brief  Construct a view of the array */
  CrsArray( const CrsArray & rhs )
    : row_map( rhs.row_map ), entries( rhs.entries ) {}

  /** \brief  Assign to a view of the rhs array.
   *          If the old view is the last view
   *          then allocated memory is deallocated.
   */
  CrsArray & operator = ( const CrsArray & rhs )
    { row_map = rhs.row_map ; entries = rhs.entries ; return *this ; }

  /**  \brief  Destroy this view of the array.
   *           If the last view then allocated memory is deallocated.
   */
  ~CrsArray() {}
};

//----------------------------------------------------------------------------

template< class CrsArrayType , class InputType >
inline
CrsArray< typename CrsArrayType::entries_type::array_type ,
          typename CrsArrayType::device_type ,
          typename CrsArrayType::row_map_type::size_type >
create_crsarray( const std::string & label ,
                 const InputType & input )
{
  typedef CrsArray< typename CrsArrayType::entries_type::array_type ,
                    typename CrsArrayType::device_type ,
                    typename CrsArrayType::row_map_type::size_type >
    output_type ;

  return Impl::Factory< output_type , InputType >
             ::create( label , input );
}

template< class CrsArrayType , class InputType >
inline
CrsArray< typename CrsArrayType::entries_type::array_type ,
          typename CrsArrayType::device_type ,
          typename CrsArrayType::row_map_type::size_type >
create_crsarray( const InputType & input )
{
  typedef CrsArray< typename CrsArrayType::entries_type::array_type ,
                    typename CrsArrayType::device_type ,
                    typename CrsArrayType::row_map_type::size_type >
    output_type ;

  return Impl::Factory< output_type , InputType >
             ::create( std::string() , input );
}

//----------------------------------------------------------------------------

template< class ArrayType ,
          class DeviceDst ,
          class DeviceSrc ,
          typename SizeType >
inline
void deep_copy(       CrsArray<ArrayType,DeviceDst,SizeType> & dst ,
                const CrsArray<ArrayType,DeviceSrc,SizeType> & src )
{
  deep_copy( dst.entries , src.entries );
  deep_copy( dst.row_map , src.row_map );
}

} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#include <impl/KokkosArray_CrsArray_factory.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_CRSARRAY_HPP */

