/** \HEADER
 *************************************************************************
 *
 *                            Kokkos
 *                 Copyright 2010 Sandia Corporation
 *
 *  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 *  the U.S. Government retains certain rights in this software.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are
 *  met:
 *
 *  1. Redistributions of source code must retain the above copyright
 *  notice, this list of conditions and the following disclaimer.
 *
 *  2. Redistributions in binary form must reproduce the above copyright
 *  notice, this list of conditions and the following disclaimer in the
 *  documentation and/or other materials provided with the distribution.
 *
 *  3. Neither the name of the Corporation nor the names of the
 *  contributors may be used to endorse or promote products derived from
 *  this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 *  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 *  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 *  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 *  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 *  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *************************************************************************
 */

#ifndef KOKKOS_CUDAMAPPEDARRAY_HPP
#define KOKKOS_CUDAMAPPEDARRAY_HPP

#if defined(__CUDA_ARCH__)
#include <cuda.h>
#endif

#include <list>

namespace Kokkos {

class CudaDevice ;
class CudaMap ;
class CudaMappedMemory ;

/*--------------------------------------------------------------------------*/
/** \brief  Handle for an array allocated and mapped onto a cuda device.
 *
 *  The array is a simple rank-2 container of simple scalar values.
 *  The first rank corresponds to an index into a block of scalars.
 *  The second rank corresponds to the parallel work index.
 *
 *  No assumptions should be made as to the mapping, contiguity, or strides
 *  of the storage of these arrays.  The mapping will vary according to the
 *  underlying device.  The usage model is for algorithms to be parameterized
 *  with respect to the type of the mapped array and thus achieve portability
 *  across compute devices.
 */
class CudaMappedArray {
public:
  typedef CudaMap     map_type ;
  typedef CudaDevice  device_type ;
  typedef int         size_type ;

  /** \brief  Access a value, only valid on the device */
  template< typename ValueType >
  __device__
  ValueType & value( size_type offset , size_type iwork ) const
    {
#if defined(__CUDA_ARCH__)
      return ((ValueType*)( m_ptr_on_device ))[ iwork + offset * m_work_count ];
#endif
    }

  /*------------------------------------------------------------------*/
  CudaMappedArray()
    : m_map_on_host( NULL )
    , m_next_on_host( this )
    , m_ptr_on_device( NULL )
    , m_work_count( 0 )
    , m_chunk_count( 0 )
    {}

  CudaMappedArray( const CudaMappedArray & rhs )
    : m_map_on_host(   NULL )
    , m_next_on_host(  rhs.m_next_on_host )
    , m_ptr_on_device( rhs.m_ptr_on_device )
    , m_work_count(    rhs.m_work_count )
    , m_chunk_count(   rhs.m_chunk_count )
    { const_cast<CudaMappedArray&>(rhs).m_next_on_host = this ; }

  CudaMappedArray & operator = ( CudaMappedArray & rhs )
    {
      prev_on_host()->m_next_on_host = m_next_on_host ;

      m_next_on_host   = rhs.m_next_on_host ; 
      m_ptr_on_device  = rhs.m_ptr_on_device ;
      m_work_count     = rhs.m_work_count ;
      m_chunk_count    = rhs.m_chunk_count ;
      const_cast<CudaMappedArray&>(rhs).m_next_on_host = this ;
      return *this ;
    }

  ~CudaMappedArray()
    {
      prev_on_host()->m_next_on_host = m_next_on_host ;

      m_map_on_host   = NULL ;
      m_next_on_host  = this ;
      m_ptr_on_device = NULL ;
      m_work_count    = 0 ;
      m_chunk_count   = 0 ;
    }

  /*------------------------------------------------------------------*/
  /** \brief  Query the map for which the array is allocated.
   *
   *  If the array view has not been allocated
   *  on a map then NULL is returned.
   */
  map_type * map() const ;

private:

  friend class CudaMap ;

  CudaMappedArray * prev_on_host()
    {
      CudaMappedArray * a = m_next_on_host ;
      while ( this != a->m_next_on_host ) { a = a->m_next_on_host ; }
      return a ;
    }

  /** \brief  Query the view of this array that is owned by a host map.
   *
   *  At most one array view is owned by a host map.
   *  That view is private to the host map.
   */ 
  CudaMappedArray * owned_view();
  
  /** \brief  Assign this array view to be owned by the host map
   *          an all associated array views to reference the
   *          allocated memory.
   */ 
  void assign( map_type * map , void * pointer ,
               size_type work_count ,
               size_type chunk_count );
    
  // The 'm_map' member is only non-null for the host map's owned view
  // A linked list is used to assign and clear allocation information
  // for all views.

  map_type        * m_map_on_host ;   ///< Cuda map for this array
  CudaMappedArray * m_next_on_host ;  ///< Ring of views
  void            * m_ptr_on_device ; ///< Pointer to array memory
  size_type         m_work_count ;    ///< Parallel work count
  size_type         m_chunk_count ;   ///< Count of items per chunk
};

/*--------------------------------------------------------------------------*/

class CudaMap {
public:
  typedef CudaDevice                   device_type ;
  typedef CudaMappedArray              mapped_array_type ;
  typedef mapped_array_type::size_type size_type ;

  CudaMap( device_type & device ,
           size_type parallel_work_length );

  ~CudaMap();

  void allocate( mapped_array_type & array ,
                 size_type values_per_chunk ,
                 size_type sizeof_value );
    
  void deallocate( mapped_array_type & array );

private:

  static void require_array_is_not_owned( mapped_array_type & );

  std::list< mapped_array_type >::iterator
    require_array_is_member( mapped_array_type & array );

  device_type                   & m_device ;
  std::list< mapped_array_type >  m_allocated_arrays ;
  size_type                       m_parallel_work_length ;

};

/*--------------------------------------------------------------------------*/

class CudaDevice {
public:
  typedef CudaMap map_type ;
};

/*--------------------------------------------------------------------------*/

} /* namespace Kokkos */

#endif /* KOKKOS_CUDAMAPPEDARRAY_HPP */

