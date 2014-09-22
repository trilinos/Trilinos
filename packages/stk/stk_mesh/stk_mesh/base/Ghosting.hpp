/*------------------------------------------------------------------------*/
/*                 Copyright (c) 2013, Sandia Corporation.
/*                 Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
/*                 the U.S. Governement retains certain rights in this software.
/*                 
/*                 Redistribution and use in source and binary forms, with or without
/*                 modification, are permitted provided that the following conditions are
/*                 met:
/*                 
/*                     * Redistributions of source code must retain the above copyright
/*                       notice, this list of conditions and the following disclaimer.
/*                 
/*                     * Redistributions in binary form must reproduce the above
/*                       copyright notice, this list of conditions and the following
/*                       disclaimer in the documentation and/or other materials provided
/*                       with the distribution.
/*                 
/*                     * Neither the name of Sandia Corporation nor the names of its
/*                       contributors may be used to endorse or promote products derived
/*                       from this software without specific prior written permission.
/*                 
/*                 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
/*                 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
/*                 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
/*                 A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
/*                 OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
/*                 SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
/*                 LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
/*                 DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
/*                 THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
/*                 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
/*                 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
/*                 
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#ifndef stk_mesh_Ghosting_hpp
#define stk_mesh_Ghosting_hpp

#include <stddef.h>                     // for size_t
#include <iosfwd>                       // for ostream
#include <stk_mesh/base/Types.hpp>      // for EntityProc
#include <string>                       // for string
#include <vector>                       // for vector
namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { struct EntityKey; } }

namespace stk {
namespace mesh {

/** \brief  Data for ghosting mesh entities.
 *
 *  This class is a member of the BulkData 'aggregate'.
 *  As such the BulkData class is the only one allowed to modify it,
 *  the "aggregate owner modifies" rule.  Thus all public methods
 *  are const and the BulkData (owner) class is a friend.
 */
class Ghosting {
public:

  /** \brief  Text name for printing purposes only */
  const std::string & name() const { return m_name ; }

  /** \brief  Ordinal to identify the ghosting subset */
  unsigned ordinal() const { return m_ordinal ; }

  /** \brief  Bulk data synchronization count when this
   *          ghosting object was last modified.
   */
  size_t synchronized_count() const { return m_sync_count ; }

  /** \brief  Locally owned entities ghosted on other processors.
   *
   *          This generated communication list for sending updates
   *          is sorted by entity key and processor rank.
   */
  void send_list( std::vector< EntityProc > & ) const ;

  /** \brief  Entities ghosted on this processor from the owner.
   *
   *          This generated communication list for receiving updates
   *          is sorted by entity key.
   */
  void receive_list( std::vector<EntityKey> & ) const ;

  /** \brief  Print the details of this object for debugging
   */
  std::ostream& operator<<(std::ostream& out) const;

private:
  /** \brief  A Ghosting object is owned by a BulkData object,
   *          and as such can only be modified by its owner.
   */

  BulkData & bulk_data() const { return m_mesh ; }
  friend class BulkData ;

  BulkData                & m_mesh ; ///< Owner
  const std::string         m_name ; ///< Name for printing purposes
  size_t                    m_sync_count ; ///< Bulk data sync count
  unsigned                  m_ordinal ;

  Ghosting( BulkData & M , const std::string & n , unsigned ord , size_t count )
    : m_mesh( M ) , m_name( n ), m_sync_count( count ), m_ordinal( ord ) {}

  ~Ghosting() {}

  // None of the following are implemented:
  Ghosting();
  Ghosting( const Ghosting & );
  Ghosting & operator = ( const Ghosting & );
};

std::ostream& operator<<(std::ostream& out, const Ghosting& rhs);

}
}

#endif

