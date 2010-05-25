/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_mesh_EntityImpl_hpp
#define stk_mesh_EntityImpl_hpp

//#include <utility>

#include <stk_mesh/base/Types.hpp>

//#include <stk_util/util/NamedPair.hpp>
#include <stk_util/util/PairIter.hpp>
#include <stk_mesh/base/Types.hpp>
//#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Relation.hpp>

namespace stk {
namespace mesh {
namespace impl {

/** \addtogroup stk_mesh_module
 * \{
 */

class EntityImpl {
public:

  ~EntityImpl();
  EntityImpl( const EntityKey & arg_key );

  // Exposed in external interface:
  EntityRank entity_rank() const { return stk::mesh::entity_rank( m_key ); }
  EntityId identifier() const { return stk::mesh::entity_id( m_key ); }
  const EntityKey & key() const { return m_key ; }
  PairIterRelation relations() const { return PairIterRelation( m_relation ); }
  PairIterRelation relations( unsigned type ) const ;
  PairIterEntityComm comm() const { return PairIterEntityComm( m_comm ); }
  PairIterEntityComm sharing() const ;
  PairIterEntityComm comm( const Ghosting & sub ) const ;
  const Bucket & bucket() const { return *m_bucket ; }
  unsigned bucket_ordinal() const { return m_bucket_ord ; }
  unsigned owner_rank() const { return m_owner_rank ; }
  size_t synchronized_count() const { return m_sync_count ; }

  // Exposed in internal interface:

  /** Change log to reflect change from before 'modification_begin'
   *  to the current status.
   */
  enum ModificationLog { LogNoChange = 0 ,
                         LogCreated  = 1 ,
                         LogModified = 2 };


  static void declare_relation( Entity & e_from, Entity & e_to, const unsigned local_id, unsigned sync_count);
  static void destroy_relation( Entity & e_from, Entity & e_to);

  // Communication info access:
  bool insert( const EntityCommInfo & );
  bool erase(  const EntityCommInfo & ); ///< Erase this entry
  bool erase(  const Ghosting & );       ///< Erase this ghosting info.
  void comm_clear_ghosting(); ///< Clear ghosting
  void comm_clear(); ///< Clear everything
  // Miscellaneous accessors:
  Bucket * get_bucket() const { return m_bucket; }

  void set_bucket_and_ordinal( Bucket * bucket, unsigned ordinal ) {
    m_bucket = bucket;
    m_bucket_ord = ordinal;
    log_modified();
  }

  void set_owner_rank( unsigned owner_rank ) {
    m_owner_rank = owner_rank;
    log_modified();
  }

  void set_sync_count( size_t sync_count ) {
    m_sync_count = sync_count;
    log_modified();
  }

  // Change log access:
  ModificationLog log_query() const { return m_mod_log ; }

  void log_clear() {
    m_mod_log = LogNoChange;
  }
  void log_created() {
    m_mod_log = LogCreated;
  }
  void log_modified() {
    if ( LogCreated != m_mod_log) {
      m_mod_log = LogModified;
    }
  }

  bool marked_for_destruction() const;

 private:

  const EntityKey         m_key ;        ///< Globally unique key
  RelationVector          m_relation ;   ///< This entity's relationships
  EntityCommInfoVector    m_comm ;       ///< This entity's communications
  Bucket                * m_bucket ;     ///< Bucket for the entity's field data
  unsigned                m_bucket_ord ; ///< Ordinal within the bucket
  unsigned                m_owner_rank ; ///< Owner processors' rank
  size_t                  m_sync_count ; ///< Last membership change
  ModificationLog         m_mod_log ;

private:
  EntityImpl();
  EntityImpl( const EntityImpl & ); ///< Copy constructor not allowed
  EntityImpl & operator = ( const EntityImpl & ); ///< Assignment operator not allowed

};

//----------------------------------------------------------------------

/** \} */

} // namespace impl
} // namespace mesh
} // namespace stk

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif /* stk_mesh_EntityImpl_hpp */

