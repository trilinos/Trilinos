/*************************************************************************
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

#include <impl/Kokkos_ViewTracker.hpp>

/** TODO : pthread safety when HAVE_PTHREAD is defined */

namespace Kokkos {
namespace Impl {

bool ViewTracker::insert( const ViewTracker & rhs )
{
  const bool ok = 0 == next ;

  if ( ok ) { // 'this' is not in a ring.

    if ( 0 != rhs.next ) { // 'rhs' is in a ring.

      // Insert in the ring.
      next = rhs.next->next ; 
      rhs.next->next = this ;
    }
    else if ( & rhs == this ) {
      next = this ; // Start a ring with this->next = this
    }
  }

  return ok ;
}

bool ViewTracker::remove_and_query_is_last()
{
  bool is_last = false ;

  if ( 0 != next ) { // This is in a ring.

    is_last = this == next ; // Is the last member in the ring.

    if ( ! is_last ) {
      // Find previous member in the ring and remove 'this'
      // A corrupted ring will iterate bad memory locations.
      ViewTracker * tmp = next ;
      for ( ; this != tmp->next ; tmp = tmp->next );
      tmp->next = next ;
    }

    next = 0 ;
  }

  return is_last ;
}

unsigned ViewTracker::test_support_view_count() const
{
  unsigned count = 0 ;

  if ( 0 != next ) {
    count = 1 ;

    for ( ViewTracker * tmp = next ;
          this != tmp ; tmp = tmp->next ) { ++count ; }
  }
  
  return count ;
}

} // namespace Impl
} // namespace Kokkos

