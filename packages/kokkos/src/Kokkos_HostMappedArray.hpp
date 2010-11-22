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

#ifndef KOKKOS_HOSTMAPPEDARRAY_HPP
#define KOKKOS_HOSTMAPPEDARRAY_HPP

#include <cstddef>
#include <list>

namespace Kokkos {

class HostDevice ;
class HostMap ;
class HostMappedArray ;

/*--------------------------------------------------------------------------*/
/** \brief  Handle for an array allocated and mapped onto a host device.
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
class HostMappedArray {
public:
  typedef HostMap     map_type ;
  typedef HostDevice  device_type ;
  typedef size_t      size_type ;

  /*------------------------------------------------------------------*/
  /** \brief  Access member value of the mapped array.  */
  template< typename ValueType >
  inline
  ValueType & value( size_type offset , size_type iwork ) const
    { return ((ValueType*)( m_pointer ))[ offset + iwork * m_number ]; }

  /*------------------------------------------------------------------*/
  /** \brief  Construct array view to NULL */
  HostMappedArray()
    : m_map( NULL )
    , m_next( this )
    , m_pointer( NULL )
    , m_number( 0 )
    {}

  /** \brief  Construct another view to the 'rhs' array view */
  HostMappedArray( const HostMappedArray & rhs )
    : m_map( NULL )
    , m_next(    rhs.m_next )
    , m_pointer( rhs.m_pointer )
    , m_number(  rhs.m_number )
    { const_cast<HostMappedArray&>(rhs).m_next = this ; }

  /** \brief  Assign this array view to the 'rhs' array view */
  HostMappedArray & operator = ( const HostMappedArray & rhs )
    {
      // Remove this view from the current family of views
      prev()->m_next = m_next ;

      // Add this family to the 'rhs' family of views
      m_next     = rhs.m_next ;
      m_pointer  = rhs.m_pointer ;
      m_number   = rhs.m_number ;
      const_cast<HostMappedArray&>(rhs).m_next = this ;
      return *this ;
    }

  /** \brief  Destroy this view to the array,
   *          the array is not deallocated.
   */
  ~HostMappedArray()
    {
      // Remove this view from the current family of views
      prev()->m_next = m_next ;

      // Clear all member data for safety
      m_map     = NULL ;
      m_next    = this ;
      m_pointer = NULL ;
      m_number  = 0 ;
    }

  /*------------------------------------------------------------------*/
  /** \brief  Query the host map for which the array is allocated.
   *
   *  If the array view has not been allocated
   *  on a host map then NULL is returned.
   */
  map_type * map() const ;

private:
  /** \brief  A host mapped array is owned, allocated, and deallocated
   *          by a host map.
   */
  friend class HostMap ;

  HostMappedArray * prev() const
    {
      HostMappedArray * a = m_next ;
      while ( this != a->m_next ) { a = a->m_next ; }
      return a ;
    }

  /** \brief  Query the view of this array that is owned by a host map.
   *
   *  At most one array view is owned by a host map.
   *  That view is private to the host map.
   */
  HostMappedArray * owned_view();

  /** \brief  Assign this array view to be owned by the host map
   *          an all associated array views to reference the
   *          allocated memory.
   */
  void assign( map_type * map , void * pointer , size_type number );

  // The 'm_map' member is only non-null for the host map's owned view
  // A linked list is used to assign and clear allocation information
  // for all views.

  map_type        * m_map ;      ///< Host map that owns the memory
  HostMappedArray * m_next ;     ///< Linked list to manage multiple views
  void            * m_pointer ;  ///< Memory allocated by host map
  size_type         m_number ;   ///< Number of values per work dimension
};

/*--------------------------------------------------------------------------*/

class HostDevice {
public:
  typedef HostMap map_type ;
  HostDevice() {}
private:
  HostDevice( const HostDevice & );
  HostDevice operator = ( const HostDevice & );
};

class HostMap {
public:
  typedef HostDevice                    device_type ;
  typedef HostMappedArray               mapped_array_type ;
  typedef mapped_array_type::size_type  size_type ;

  size_type parallel_work_length() const
    { return m_parallel_work_length ; }

  /** \brief  Map arrays of parallel_work_length onto device.  */
  HostMap( device_type & device ,
           size_type parallel_work_length );

  /** \brief  Destroy the host map and all of its allocated arrays. */
  ~HostMap();

  /** \brief  Allocate a mapped array owned by this host map object.
   *
   *  The array will have 'parallel_work_length * values_per_chunk'
   *  simple scalar values where each value is 'sizeof_value' bytes.
   *
   *  Precondition: The array cannot already be allocated.
   */
  void allocate( mapped_array_type & array ,
                 size_type values_per_chunk ,
                 size_type sizeof_value );

  /** \brief  Deallocate a mapped array owned by this host map object.
   *
   *  A warning will be issued if there exists a view to the array
   *  other than the input view.
   */
  void deallocate( mapped_array_type & array );

private:

  /** \brief  Require that the array is not owned by any host map */
  static void require_array_is_not_owned( mapped_array_type & );

  /** \brief  Require that the array is owned by this host map
   *          and return the owned member from 'm_allocated_arrays'.
   */
  std::list< mapped_array_type >::iterator
    require_array_is_member( mapped_array_type & );

  HostMap();
  HostMap( const HostMap & );
  HostMap & operator = ( const HostMap & );

  device_type                  & m_device ;
  std::list<mapped_array_type>   m_allocated_arrays ;
  size_type                      m_parallel_work_length ;
};

/*--------------------------------------------------------------------------*/

} // namespace Kokkos

#endif /* KOKKOS_HOSTMAPPEDARRAY_HPP */


