
#ifndef stk_mesh_Ghosting_hpp
#define stk_mesh_Ghosting_hpp

#include <vector>
#include <string>
#include <stk_mesh/base/Types.hpp>

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
  BulkData & mesh() const { return m_mesh ; }

  /** \brief  Text name for printing purposes only */
  const std::string & name() const { return m_name ; }

  /** \brief  Bulk data synchronization count when this
   *          ghosting object was last modified.
   */
  size_t synchronized_count() const { return m_sync_count ; }

  /** \brief  Locally owned entities ghosted on other processors.
   *
   *          This communication list for sending updates
   *          is sorted by entity key and processor rank.
   *
   * \todo REFACTOR rename to 'send_list'
   */
  const std::vector< EntityProc > & send() const { return m_send ; }

  /** \brief  Entities ghosted on this processor from the owner.
   *
   *          This communication list for receiving updates
   *          is sorted by entity key.
   *
   * \todo REFACTOR rename to 'receive_list'
   */
  const std::vector< Entity * > & receive() const { return m_recv ; }

private:
  /** \brief  A Ghosting object is owned by a BulkData object,
   *          and as such can only be modified by its owner.
   */
  friend class BulkData ;

  BulkData                & m_mesh ; ///< Owner
  const std::string         m_name ; ///< Name for printing purposes
  std::vector< EntityProc > m_send ; ///< Locally owned entities to send
  std::vector< Entity * >   m_recv ; ///< Ghosted entities
  size_t                    m_sync_count ; ///< Bulk data sync count

  Ghosting( BulkData & M , const std::string & n , size_t count )
    : m_mesh( M ) , m_name( n ), m_send(), m_recv(), m_sync_count( count ) {}

  ~Ghosting() {}

  // None of the following are implemented:
  Ghosting();
  Ghosting( const Ghosting & );
  Ghosting & operator = ( const Ghosting & );
};

}
}

#endif

