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

#ifndef KOKKOS_BASEMAPPEDARRAY_HPP
#define KOKKOS_BASEMAPPEDARRAY_HPP

#include <cstddef>
#include <list>

#define KOKKOS_MDARRAY_CHECK( EXPR ) /* */

namespace Kokkos {

/** \brief  Interface for device map */
class BaseMapInterface {
public:
  typedef size_t size_type ;
  virtual size_type parallel_work_count() const = 0 ;
  virtual ~BaseMapInterface() {}
};

/** \brief  Base data for contiguous multidimensional array */
class BaseMappedArray {
public:
  typedef size_t size_type ;

  /** \brief  Construct a NULL view */
  BaseMappedArray();

  /** \brief  Construct another view of the array */
  BaseMappedArray( const BaseMappedArray & rhs );

  /** \brief  Assign a view of the array */
  BaseMappedArray & operator = ( const BaseMappedArray & rhs );

  /**  Destroy this view of the array */
  ~BaseMappedArray();

  /** \brief  Require the array view is not already allocated */
  void require_not_allocated() const ;

  /** \brief  Clear this array and all of its views of the
   *          allocated memory, return the pointer to the
   *          allocated memory.  Remove the owned view from the list.
   *
   *  Required for the array view to be a member of the input mapped arrays.
   */
  void * clear( std::list< BaseMappedArray > & mapped_arrays );

  /** \brief  Assign this array and all of its views to the
   *          allocated memory and insert an owned view into
   *          the list.
   */
  void assign( std::list< BaseMappedArray > & mapped_arrays ,
               BaseMapInterface * const map ,
               void * const ptr_on_device ,
               const size_type rank ,
               const size_type * const dimension );

  /*------------------------------------------------------------------*/

private:

  template< typename , class > friend class MDArray ;

  /*------------------------------------------------------------------*/

  void require_in_bounds( ptrdiff_t i0 , ptrdiff_t i1 ,
                          ptrdiff_t i2 , ptrdiff_t i3 ,
                          ptrdiff_t i4 , ptrdiff_t i5 ,
                          ptrdiff_t i6 , ptrdiff_t iP );

  void require_in_bounds( ptrdiff_t i0 , ptrdiff_t i1 ,
                          ptrdiff_t i2 , ptrdiff_t i3 ,
                          ptrdiff_t i4 , ptrdiff_t i5 ,
                          ptrdiff_t iP );

  void require_in_bounds( ptrdiff_t i0 , ptrdiff_t i1 ,
                          ptrdiff_t i2 , ptrdiff_t i3 ,
                          ptrdiff_t i4 , ptrdiff_t iP );

  void require_in_bounds( ptrdiff_t i0 , ptrdiff_t i1 ,
                          ptrdiff_t i2 , ptrdiff_t i3 ,
                          ptrdiff_t iP );

  void require_in_bounds( ptrdiff_t i0 , ptrdiff_t i1 ,
                          ptrdiff_t i2 , ptrdiff_t iP );

  void require_in_bounds( ptrdiff_t i0 , ptrdiff_t i1 ,
                          ptrdiff_t iP );

  void require_in_bounds( ptrdiff_t i0 , ptrdiff_t iP );

  void require_in_bounds( ptrdiff_t iP );

  /*------------------------------------------------------------------*/
  /** \brief  Query the map interface from the owned view */
  BaseMapInterface * base_map() const ;

  /** \brief  Query the map interface from the owned view */
  template< class MapType >
  MapType * map() const { return dynamic_cast<MapType*>( base_map() ); }

  /** \brief  Require this is not the owned view */
  void require_not_owned_view() const ;

  /** \brief  Find previous view: this == prev_on_host()->m_next_on_host */
  BaseMappedArray * prev_on_host() const ;

  /*------------------------------------------------------------------*/

  BaseMapInterface * m_map_on_host ;
  BaseMappedArray  * m_next_on_host ;
  void             * m_ptr_on_device ;
  size_type          m_rank ;
  size_type          m_dimension[8] ;

  /*------------------------------------------------------------------*/
};

/*------------------------------------------------------------------------*/

} // namespace Kokkos

#endif /* KOKKOS_MAPPEDARRAY_HPP */

