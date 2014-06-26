/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#ifndef stk_mesh_Ghosting_hpp
#define stk_mesh_Ghosting_hpp

#include <vector>
#include <string>
#include <iosfwd>
#include <stk_mesh/base/Types.hpp>

namespace stk_classic {
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
  void receive_list( std::vector< Entity * > & ) const ;

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

