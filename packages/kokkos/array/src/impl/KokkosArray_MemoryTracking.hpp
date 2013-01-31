/*
//@HEADER
// ************************************************************************
// 
//   KokkosArray: Manycore Performance-Portable Multidimensional Arrays
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
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOSARRAY_MEMORY_TRACKING_HPP
#define KOKKOSARRAY_MEMORY_TRACKING_HPP

#include <cstddef>
#include <utility>
#include <vector>
#include <string>
#include <typeinfo>
#include <ostream>

namespace KokkosArray {
namespace Impl {

class MemoryTracking {
public:

  struct Info {
    std::string            label ;
    const void           * begin ;
    const std::type_info * type ;
    size_t                 size ;
    size_t                 length ;
    size_t                 count ;

    Info()
    : label(), begin(0), type(0)
    , size(0), length(0), count(0)
    {}

    Info( const Info & rhs )
    : label(rhs.label), begin(rhs.begin), type(rhs.type)
    , size(rhs.size), length(rhs.length), count(rhs.count)
    {}

    Info & operator = ( const Info & rhs )
    {
      label = rhs.label ; begin = rhs.begin ; type = rhs.type ;
      size  = rhs.size ;  length = rhs.length ; count = rhs.count ;
      return *this ;
    }

    ~Info()
    { begin = 0 ; type = 0 ; size  = 0 ; length = 0 ; count = 0 ; }

    void print( std::ostream & ) const ;
  };

  /** \brief  Track a pointer. */
  void track( const void           * ptr ,
              const std::type_info * type ,
              const size_t           size ,
              const size_t           length ,
              const std::string      label );

  /** \brief  Track a pointer. */
  template< typename Type >
  void track( const Type      * ptr ,
              const size_t      length ,
              const std::string label )
  { track( ptr , & typeid(Type) , sizeof(Type) , length , label ); }

  /** \brief  Increment the tracking count.  */
  void increment( const void * ptr );

  /** \brief  Decrement the tracking count.
   *          If zero then the entry is deleted and the
   *          allocated pointer is returned.
   */
  void * decrement( const void * ptr );

  /** \brief  Query a tracked pointer */
  Info query( const void * ptr ) const ;

  /** \brief  Print tracked pointer information */
  void print( std::ostream & , const std::string & lead ) const ;

  /** \brief  Query if empty of tracked pointers.
   *
   *  Intent: A memory manager destructor queries if non-empty
   *  which would indicate memory leaks.
   */
  bool empty() const ;

  explicit MemoryTracking( const std::string & space );

  ~MemoryTracking();

private:
  MemoryTracking();
  MemoryTracking( const MemoryTracking & );
  MemoryTracking & operator = ( const MemoryTracking & );

  std::string                m_space ;
  std::vector<Info>          m_tracking ;
  std::vector<const void *>  m_tracking_end ;
};


} /* namespace Impl */
} /* namespace KokkosArray */

#endif

