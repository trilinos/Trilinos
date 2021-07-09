// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

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
  void receive_list( std::vector<EntityKey> & keys) const ;
  void receive_list( std::vector<Entity> & entities) const ;

  /** \brief  Print the details of this object for debugging
   */
  std::ostream& operator<<(std::ostream& out) const;

  const BulkData & mesh() const { return m_mesh ; }
        BulkData & mesh()       { return m_mesh ; }
private:
  /** \brief  A Ghosting object is owned by a BulkData object,
   *          and as such can only be modified by its owner.
   */

  friend class BulkData ;

  BulkData                & m_mesh ; ///< Owner
  const std::string         m_name ; ///< Name for printing purposes
  unsigned                  m_ordinal ;

  Ghosting( BulkData & M , const std::string & n , unsigned ord )
    : m_mesh( M ) , m_name( n ), m_ordinal( ord ) {}

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

