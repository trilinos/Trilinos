// NOTE: This is not a stand alone header file. It is meant to help clean up BulkData.hpp
// This is currently just housing the mutli-lined inlined methods.

#ifndef BULKDATA_INLINED_METHODS_HPP
#define BULKDATA_INLINED_METHODS_HPP


namespace stk {
namespace mesh {


///////////////////////////////////////////////////////////////////////////////
// The following get_connectivity API is designed for *internal* STK_Mesh usage
// to support algorithms that must always be able to get connectivity even
// when it is disabled in the connectivity map. These functions are designed
// to support connectivity-retrieval callsites that need to work regardless of
// the connectivity map. The scratch vectors will be used for allocation when it
// is needed; otherwise, they are ignored. The overloads are provided for the
// common cases where not all connectivity data is needed.
///////////////////////////////////////////////////////////////////////////////

size_t get_connectivity( const BulkData & mesh,
                         Entity entity,
                         EntityRank to_rank,
                         EntityVector & entity_scratch_storage );

size_t get_connectivity( const BulkData & mesh,
                         Entity entity, EntityRank to_rank,
                         EntityVector & entity_scratch_storage,
                         std::vector<ConnectivityOrdinal> & ordinal_scratch_storage );

size_t get_connectivity( const BulkData & mesh,
                         Entity entity,
                         EntityRank to_rank,
                         EntityVector & entity_scratch_storage,
                         std::vector<Permutation> & permutation_scratch_storage );

size_t get_connectivity( const BulkData & mesh,
                         Entity entity,
                         EntityRank to_rank,
                         EntityVector & entity_scratch_storage,
                         std::vector<ConnectivityOrdinal> & ordinal_scratch_storage,
                         std::vector<Permutation> & permutation_scratch_storage );


 /** \brief  Comparator functor for entities compares the entities' keys */
inline
EntityLess::EntityLess(const BulkData& mesh) : m_mesh(&mesh) {}

/** \brief  Comparison operator */
inline
bool EntityLess::operator()(const Entity lhs, const Entity rhs) const
{
  const EntityKey lhs_key = m_mesh->in_index_range(lhs) ? m_mesh->entity_key(lhs) : EntityKey();
  const EntityKey rhs_key = m_mesh->in_index_range(rhs) ? m_mesh->entity_key(rhs) : EntityKey();
  return lhs_key < rhs_key;
}

inline
bool EntityLess::operator()(const Entity lhs, const EntityKey & rhs) const
{
  const EntityKey lhs_key = m_mesh->in_index_range(lhs) ? m_mesh->entity_key(lhs) : EntityKey();
  return lhs_key < rhs;
}

inline
bool EntityLess::operator()( const EntityProc & lhs, const EntityProc & rhs) const
{
  const EntityKey lhs_key = m_mesh->in_index_range(lhs.first) ? m_mesh->entity_key(lhs.first) : EntityKey() ;
  const EntityKey rhs_key = m_mesh->in_index_range(rhs.first) ? m_mesh->entity_key(rhs.first) : EntityKey() ;
  return lhs_key != rhs_key ? lhs_key < rhs_key : lhs.second < rhs.second ;
}

inline
bool EntityLess::operator()( const EntityProc & lhs, const Entity rhs) const
{
  const EntityKey lhs_key = m_mesh->in_index_range(lhs.first) ? m_mesh->entity_key(lhs.first) : EntityKey();
  const EntityKey rhs_key = m_mesh->in_index_range(rhs)       ? m_mesh->entity_key(rhs)       : EntityKey();
  return lhs_key < rhs_key;
}

inline
bool EntityLess::operator()( const EntityProc & lhs, const EntityKey & rhs) const
{
  const EntityKey lhs_key = m_mesh->in_index_range(lhs.first) ? m_mesh->entity_key(lhs.first) : EntityKey();
  return lhs_key < rhs ;
}

inline
EntityLess& EntityLess::operator=(const EntityLess& rhs)
{
  m_mesh = rhs.m_mesh;
  return *this;
}

inline
BulkData & BulkData::get( const Bucket & bucket) {
  return bucket.bulk_data();
}

inline
BulkData & BulkData::get( const Ghosting & ghost) {
  return ghost.bulk_data();
}

inline
unsigned BulkData::num_connectivity(Entity entity, EntityRank rank) const
{
  ThrowAssert(bucket_ptr(entity));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->num_connectivity(mesh_idx.bucket_ordinal, rank);
}

inline
unsigned BulkData::find_ordinal(Entity entity, EntityRank rank, ConnectivityOrdinal ordinal) const
{
  ThrowAssert(bucket_ptr(entity));
  const MeshIndex &mesh_idx = mesh_index(entity);
  unsigned num_rels = mesh_idx.bucket->num_connectivity(mesh_idx.bucket_ordinal, rank);
  ConnectivityOrdinal const *ords = mesh_idx.bucket->begin_ordinals(mesh_idx.bucket_ordinal, rank);

  unsigned i = 0;
  for (; i < num_rels; ++i)
  {
    if (ords[i] == ordinal)
      break;
  }
  return i;
}

inline
Entity const* BulkData::begin(Entity entity, EntityRank rank) const
{
  ThrowAssert(is_valid_connectivity(entity, rank));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->begin(mesh_idx.bucket_ordinal, rank);
}

inline
Entity const* BulkData::begin_nodes(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::NODE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->begin_nodes(mesh_idx.bucket_ordinal);
}

inline
Entity const* BulkData::begin_edges(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::EDGE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->begin_edges(mesh_idx.bucket_ordinal);
}

inline
Entity const* BulkData::begin_faces(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::FACE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->begin_faces(mesh_idx.bucket_ordinal);
}

inline
Entity const* BulkData::begin_elements(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::ELEMENT_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->begin_elements(mesh_idx.bucket_ordinal);
}

inline
ConnectivityOrdinal const* BulkData::begin_ordinals(Entity entity, EntityRank rank) const
{
  ThrowAssert(is_valid_connectivity(entity, rank));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->begin_ordinals(mesh_idx.bucket_ordinal, rank);
}

inline
ConnectivityOrdinal const* BulkData::begin_node_ordinals(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::NODE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->begin_node_ordinals(mesh_idx.bucket_ordinal);
}

inline
ConnectivityOrdinal const* BulkData::begin_edge_ordinals(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::EDGE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->begin_edge_ordinals(mesh_idx.bucket_ordinal);
}

inline
ConnectivityOrdinal const* BulkData::begin_face_ordinals(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::FACE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->begin_face_ordinals(mesh_idx.bucket_ordinal);
}

inline
ConnectivityOrdinal const* BulkData::begin_element_ordinals(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::ELEMENT_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->begin_element_ordinals(mesh_idx.bucket_ordinal);
}

inline
Permutation const* BulkData::begin_permutations(Entity entity, EntityRank rank) const
{
  ThrowAssert(is_valid_connectivity(entity, rank));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->begin_permutations(mesh_idx.bucket_ordinal, rank);
}

inline
Permutation const* BulkData::begin_node_permutations(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::NODE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->begin_node_permutations(mesh_idx.bucket_ordinal);
}

inline
Permutation const* BulkData::begin_edge_permutations(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::EDGE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->begin_edge_permutations(mesh_idx.bucket_ordinal);
}

inline
Permutation const* BulkData::begin_face_permutations(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::FACE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->begin_face_permutations(mesh_idx.bucket_ordinal);
}

inline
Permutation const* BulkData::begin_element_permutations(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::ELEMENT_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->begin_element_permutations(mesh_idx.bucket_ordinal);
}

inline
unsigned BulkData::num_nodes(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::NODE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->num_nodes(mesh_idx.bucket_ordinal);
}

inline
unsigned BulkData::num_edges(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::EDGE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->num_edges(mesh_idx.bucket_ordinal);
}

inline
unsigned BulkData::num_faces(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::FACE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->num_faces(mesh_idx.bucket_ordinal);
}

inline
unsigned BulkData::num_elements(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::ELEMENT_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->num_elements(mesh_idx.bucket_ordinal);
}

inline
Entity const* BulkData::end(Entity entity, EntityRank rank) const
{
  ThrowAssert(is_valid_connectivity(entity, rank));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->end(mesh_idx.bucket_ordinal, rank);
}

inline
Entity const* BulkData::end_nodes(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::NODE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->end_nodes(mesh_idx.bucket_ordinal);
}

inline
Entity const* BulkData::end_edges(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::EDGE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->end_edges(mesh_idx.bucket_ordinal);
}

inline
Entity const* BulkData::end_faces(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::FACE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->end_faces(mesh_idx.bucket_ordinal);
}

inline
Entity const* BulkData::end_elements(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::ELEMENT_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->end_elements(mesh_idx.bucket_ordinal);
}

inline
ConnectivityOrdinal const* BulkData::end_ordinals(Entity entity, EntityRank rank) const
{
  ThrowAssert(is_valid_connectivity(entity, rank));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->end_ordinals(mesh_idx.bucket_ordinal, rank);
}

inline
ConnectivityOrdinal const* BulkData::end_node_ordinals(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::NODE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->end_node_ordinals(mesh_idx.bucket_ordinal);
}

inline
ConnectivityOrdinal const* BulkData::end_edge_ordinals(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::EDGE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->end_edge_ordinals(mesh_idx.bucket_ordinal);
}

inline
ConnectivityOrdinal const* BulkData::end_face_ordinals(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::FACE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->end_face_ordinals(mesh_idx.bucket_ordinal);
}

inline
ConnectivityOrdinal const* BulkData::end_element_ordinals(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::ELEMENT_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->end_element_ordinals(mesh_idx.bucket_ordinal);
}

inline
Permutation const* BulkData::end_permutations(Entity entity, EntityRank rank) const
{
  ThrowAssert(is_valid_connectivity(entity, rank));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->end_permutations(mesh_idx.bucket_ordinal, rank);
}

inline
Permutation const* BulkData::end_node_permutations(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::NODE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->end_node_permutations(mesh_idx.bucket_ordinal);
}

inline
Permutation const* BulkData::end_edge_permutations(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::EDGE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->end_edge_permutations(mesh_idx.bucket_ordinal);
}

inline
Permutation const* BulkData::end_face_permutations(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::FACE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->end_face_permutations(mesh_idx.bucket_ordinal);
}

inline
Permutation const* BulkData::end_element_permutations(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::ELEMENT_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->end_element_permutations(mesh_idx.bucket_ordinal);
}

inline
bool BulkData::has_permutation(Entity entity, EntityRank rank) const
{
  ThrowAssert(bucket_ptr(entity));
  return bucket(entity).has_permutation(rank);
}

inline
int BulkData::entity_comm_map_owner(const EntityKey & key) const
{
    return internal_entity_comm_map_owner(key);
}
inline
int BulkData::internal_entity_comm_map_owner(const EntityKey & key) const
{
  const int owner_rank = m_entity_comm_map.owner_rank(key);
  ThrowAssertMsg(owner_rank == InvalidProcessRank || owner_rank == parallel_owner_rank(get_entity(key)),
                 "Expected entity " << key.id() << " with rank " << key.rank() << " to have owner " <<
                 parallel_owner_rank(get_entity(key)) << " but in comm map, found " << owner_rank);
  return owner_rank;
}

inline
bool BulkData::in_receive_ghost( EntityKey key ) const
{
  const std::vector<Ghosting*> & ghosts= ghostings();
  for (size_t i=ghosts.size()-1;i>=AURA;--i)
  {
      if ( in_receive_ghost(*ghosts[i], key) )
          return true;
  }
  return false;
}

inline
bool BulkData::in_receive_ghost( const Ghosting & ghost , EntityKey key ) const
{
  const int owner_rank = internal_entity_comm_map_owner(key);
  return in_ghost( ghost , key , owner_rank );
}

inline
bool BulkData::in_send_ghost( EntityKey key) const
{
    const int owner_rank = internal_entity_comm_map_owner(key);
    for ( PairIterEntityComm ec = internal_entity_comm_map(key); ! ec.empty() ; ++ec )
    {
      if ( ec->ghost_id != 0 &&
           ec->proc     != owner_rank)
      {
        return true;
      }
    }
    return false;
}

inline
void BulkData::internal_check_unpopulated_relations(Entity entity, EntityRank rank) const
{
#ifndef NDEBUG
  if (m_check_invalid_rels) {
    const MeshIndex &mesh_idx = mesh_index(entity);
    const Bucket &b = *mesh_idx.bucket;
    Bucket::size_type bucket_ord = mesh_idx.bucket_ordinal;
    ThrowAssert(count_valid_connectivity(entity, rank) == b.num_connectivity(bucket_ord, rank));
  }
#endif
}

struct EntityGhostData
{
    enum DIRECTION {
        INVALID,
        NONE,
        SEND,
        RECEIVE
    };
    enum GHOST_LEVEL {
        LOCALLY_OWNED = -1,
        SHARED = 0,
        AURA = 1
    };
    DIRECTION direction;
    int ghostingLevel;
    int processor;
    Entity entity;
    const BulkData * bulkData;

    EntityGhostData()
    : direction(INVALID)
    , ghostingLevel(-2)
    , processor(-1)
    , entity()
    , bulkData(NULL) { }

    // insert to any ostream-like s
    template<class OStream> friend inline OStream& operator << (OStream& s, const DIRECTION& dir)
    {
        switch (dir) {
            case INVALID:
                s << "INVALID";
                break;
            case NONE:
                s << "NONE";
                break;
            case SEND:
                s << "SEND";
                break;
            case RECEIVE:
                s << "RECEIVE";
                break;
            default:
                s << "INVALID";
        }
        return s;
    }
    template<class OStream> inline OStream& printGhostLevel(OStream& s, int gl) const
    {
        switch (gl) {
            case LOCALLY_OWNED:
                s << "LOCALLY_OWNED";
                break;
            case SHARED:
                s << "SHARED";
                break;
            case AURA:
                s << "AURA";
                break;
            default:
                s << "CUSTOM_" << (gl-1);
        }
        return s;
    }
    template<class OStream> friend inline OStream& operator << (OStream& s, const EntityGhostData& egd)
    {
        if (egd.bulkData != NULL) {
            s << "(Entity_gid=";
            s << egd.bulkData->identifier(egd.entity)
              << ", rank=" << static_cast<unsigned int>(egd.bulkData->entity_rank(egd.entity));
        }
        else {
            s << "(Entity_lid=";
            s << egd.entity;
        }
        s << ", direction=" << egd.direction
          << ", processor=" << egd.processor
          << ", ghosting level=";
        egd.printGhostLevel(s,egd.ghostingLevel);
        s << ")";
        return s;
    }
};

struct StoreEntityKeyInSet
{
    StoreEntityKeyInSet(const BulkData & mesh_in) : mesh(mesh_in) {}
    void operator()(Entity entity) {
       entity_key_set.insert(mesh.entity_key(entity));
    }
    std::set<EntityKey> entity_key_set;
    const BulkData & mesh;
};

struct StoreEntityProcInSet
{
    StoreEntityProcInSet(const BulkData & mesh_in)
    :mesh(mesh_in)
    ,entity_proc_set(EntityLess(mesh_in))
    ,target(-1) {}

    bool operator()(Entity entity) {
        EntityProc ep(entity,target);
        if (entity_proc_set.find(ep) == entity_proc_set.end()) {
            entity_proc_set.insert(ep);
            return true;
        }
        return false;
    }
    const BulkData & mesh;
    std::set<EntityProc,EntityLess> entity_proc_set;
    int target;
};

////////////////

inline void BulkData::copy_entity_fields( Entity src, Entity dst)
{
  //TODO fix const correctness for src
  MeshIndex & src_mesh_idx = mesh_index(src);
  MeshIndex & dst_mesh_idx = mesh_index(dst);

  //// Pre-upgrade stk_mesh did not have this restriction, and it was easy enough to remove.
  //    ThrowAssert(src_mesh_idx.bucket->entity_rank() == dst_mesh_idx.bucket->entity_rank());

  copy_entity_fields_callback(dst_mesh_idx.bucket->entity_rank(),
                       dst_mesh_idx.bucket->bucket_id(),
                       dst_mesh_idx.bucket_ordinal,
                       src_mesh_idx.bucket->bucket_id(),
                       src_mesh_idx.bucket_ordinal);
}

inline Entity BulkData::get_entity( EntityRank ent_rank , EntityId entity_id ) const
{
  if (!is_good_rank_and_id(ent_rank, entity_id)) {
      return Entity();
  }
  return m_entity_repo.get_entity( EntityKey(ent_rank, entity_id));
}

inline bool BulkData::relation_exist( const Entity entity, EntityRank subcell_rank, RelationIdentifier subcell_id )
{
  bool found = false;
  Entity const * rel_entity_it = bucket(entity).begin(bucket_ordinal(entity),subcell_rank);
  const unsigned num_rel = bucket(entity).num_connectivity(bucket_ordinal(entity),subcell_rank);
  ConnectivityOrdinal const * rel_ord_it = bucket(entity).begin_ordinals(bucket_ordinal(entity),subcell_rank);

  for (unsigned i=0 ; i < num_rel ; ++i) {
    if (rel_ord_it[i] == static_cast<ConnectivityOrdinal>(subcell_id) && is_valid(rel_entity_it[i])) {
      found = true;
      break;
    }
  }

  return found;
}

inline bool BulkData::element_side_polarity( const Entity elem ,
      const Entity side , unsigned local_side_id ) const
{
    // 09/14/10:  TODO:  tscoffe:  Will this work in 1D??
    const bool is_side = entity_rank(side) != stk::topology::EDGE_RANK;
    stk::topology elem_top = bucket(elem).topology();

    const unsigned side_count = ! (elem_top != stk::topology::INVALID_TOPOLOGY) ? 0 : (
        is_side ? elem_top.num_sides()
            : elem_top.num_edges() );

    ThrowErrorMsgIf( elem_top == stk::topology::INVALID_TOPOLOGY,
        "For Element[" << identifier(elem) << "], element has no defined topology");

    ThrowErrorMsgIf( static_cast<unsigned>(side_count) <= local_side_id,
        "For Element[" << identifier(elem) << "], " <<
        "side: " << identifier(side) << ", " <<
        "local_side_id = " << local_side_id <<
        " ; unsupported local_side_id");

    stk::topology side_top =
        is_side ? elem_top.side_topology( local_side_id )
            : elem_top.sub_topology( stk::topology::EDGE_RANK, local_side_id );

    std::vector<unsigned> side_map(side_top.num_nodes());
    elem_top.side_node_ordinals( local_side_id, side_map.begin());

    Entity const *elem_nodes = begin_nodes(elem);
    Entity const *side_nodes = begin_nodes(side);
    const unsigned n = side_top.num_nodes();
    bool good = false ;
    for ( unsigned i = 0 ; !good && i < n ; ++i ) {
        good = true;
        for ( unsigned j = 0; good && j < n ; ++j ) {
          good = side_nodes[(j+i)%n] == elem_nodes[ side_map[j] ];
        }
    }
    return good ;
}

inline VolatileFastSharedCommMapOneRank const& BulkData::volatile_fast_shared_comm_map(EntityRank rank) const
{
  ThrowAssert(synchronized_state() == SYNCHRONIZED);
  ThrowAssertMsg(rank < stk::topology::ELEMENT_RANK, "Cannot shared entities of rank: " << rank);
  return m_volatile_fast_shared_comm_map[rank];
}

inline Part& BulkData::ghosting_part(const Ghosting& ghosting) const
{
  ThrowRequireMsg(ghosting.ordinal() < m_ghost_parts.size(), "BulkData::ghosting_part ERROR, no part corresponds to ghosting with name="<<ghosting.name()<<" and ordinal="<<ghosting.ordinal());
  return *m_ghost_parts[ghosting.ordinal()];
}

inline bool BulkData::in_index_range(Entity entity) const
{
  return entity.local_offset() < m_entity_states.size();
}

inline bool BulkData::is_valid(Entity entity) const
{
  return (entity.local_offset() < m_entity_states.size()) && (m_entity_states[entity.local_offset()] != Deleted);
}

inline const MeshIndex& BulkData::mesh_index(Entity entity) const
{
#ifndef NDEBUG
  entity_getter_debug_check(entity);
#endif

  return m_mesh_indexes[entity.local_offset()];
}

inline MeshIndex& BulkData::mesh_index(Entity entity)
{
#ifndef NDEBUG
  entity_setter_debug_check(entity); // setter check due to non-const
#endif

  return m_mesh_indexes[entity.local_offset()];
}

inline EntityId BulkData::identifier(Entity entity) const
{
  entity_getter_debug_check(entity);

  return m_entity_keys[entity.local_offset()].id();
}

inline EntityRank BulkData::entity_rank(Entity entity) const
{
  entity_getter_debug_check(entity);

  return m_entity_keys[entity.local_offset()].rank();
}

inline EntityKey BulkData::entity_key(Entity entity) const
{
  entity_getter_debug_check(entity);

  return m_entity_keys[entity.local_offset()];
}

inline EntityState BulkData::state(Entity entity) const
{
  entity_getter_debug_check(entity);

  return static_cast<EntityState>(m_entity_states[entity.local_offset()]);
}

inline void BulkData::internal_mark_entity(Entity entity, entitySharing sharedType)
{
    m_mark_entity[entity.local_offset()] = static_cast<int>(sharedType);
}

inline BulkData::entitySharing BulkData::internal_is_entity_marked(Entity entity) const
{
    return static_cast<entitySharing>(m_mark_entity[entity.local_offset()]);
}

inline bool BulkData::internal_add_node_sharing_called() const
{
  return m_add_node_sharing_called;
}

inline Bucket & BulkData::bucket(Entity entity) const
{
  entity_getter_debug_check(entity);

  return *mesh_index(entity).bucket;
}

inline Bucket * BulkData::bucket_ptr(Entity entity) const
{
  entity_getter_debug_check(entity);

  return mesh_index(entity).bucket;
}

inline Bucket::size_type BulkData::bucket_ordinal(Entity entity) const
{
  entity_getter_debug_check(entity);

  return mesh_index(entity).bucket_ordinal;
}

inline int BulkData::parallel_owner_rank(Entity entity) const
{
  entity_getter_debug_check(entity);

  return bucket(entity).parallel_owner_rank(bucket_ordinal(entity));
}

inline unsigned BulkData::local_id(Entity entity) const
{
  entity_getter_debug_check(entity);

  return m_local_ids[entity.local_offset()];
}

#ifdef SIERRA_MIGRATION
inline BulkData::FmwkId BulkData::global_id(stk::mesh::Entity entity) const
{
  entity_getter_debug_check(entity);

  return m_fmwk_global_ids[entity.local_offset()];
}

inline const RelationVector& BulkData::aux_relations(Entity entity) const
{
  ThrowAssert(m_add_fmwk_data);
  entity_setter_debug_check(entity); // setter check due to side effects

  if (m_fmwk_aux_relations[entity.local_offset()] == NULL) {
    m_fmwk_aux_relations[entity.local_offset()] = new RelationVector();
  }
  return *m_fmwk_aux_relations[entity.local_offset()];
}

inline RelationVector& BulkData::aux_relations(Entity entity)
{
  ThrowAssert(m_add_fmwk_data);
  entity_setter_debug_check(entity); // setter check due to side effects

  if (m_fmwk_aux_relations[entity.local_offset()] == NULL) {
    m_fmwk_aux_relations[entity.local_offset()] = new RelationVector();
  }
  return *m_fmwk_aux_relations[entity.local_offset()];
}

inline void BulkData::set_global_id(stk::mesh::Entity entity, int id)
{
  entity_setter_debug_check(entity);

  m_fmwk_global_ids[entity.local_offset()] = id;
}

inline RelationIterator BulkData::internal_begin_relation(Entity entity, const Relation::RelationType relation_type) const
{
  ThrowAssert(m_add_fmwk_data);
  if (impl::internal_is_handled_generically(relation_type)) {
    ThrowErrorMsg("stk::Mesh::BulkData::internal_begin_relation(..) requests native stk::mesh relation type");
    return RelationIterator();
  }
  else {
    return aux_relations(entity).begin();
  }
}

inline RelationIterator BulkData::internal_end_relation(Entity entity, const Relation::RelationType relation_type) const
{
  ThrowAssert(m_add_fmwk_data);
  if (impl::internal_is_handled_generically(relation_type)) {
    ThrowErrorMsg("stk::Mesh::BulkData::internal_begin_relation(..) requests native stk::mesh relation type");
    return RelationIterator();
  }
  else {
    return aux_relations(entity).end();
  }
}

inline void BulkData::compress_relation_capacity(Entity entity)
{
  RelationVector &rels = aux_relations(entity);
  RelationVector tmp(rels);
  tmp.swap(rels);
}
#endif

inline void BulkData::set_mesh_index(Entity entity, Bucket * in_bucket, Bucket::size_type ordinal )
{
  // The trace statement forces this method to be defined after Entity
  TraceIfWatching("stk::mesh::BulkData::set_mesh_index", LOG_ENTITY, entity_key(entity));

  entity_setter_debug_check(entity);

  if (in_bucket != NULL) {
    ThrowAssertMsg(in_bucket->size() >= ordinal, "Detected bad bucket/ordinal.");
  }
  MeshIndex &mesh_idx = mesh_index(entity);
  mesh_idx.bucket = in_bucket;
  mesh_idx.bucket_ordinal = ordinal;
}

inline void BulkData::set_entity_key(Entity entity, EntityKey key)
{
  entity_setter_debug_check(entity);

  m_entity_keys[entity.local_offset()] = key;
}

inline void BulkData::set_state(Entity entity, EntityState entity_state)
{
  entity_setter_debug_check(entity);

  m_entity_states[entity.local_offset()] = static_cast<uint16_t>(entity_state);
  m_mark_entity[entity.local_offset()] = NOT_MARKED;
}

inline void BulkData::set_local_id(Entity entity, unsigned id)
{
  entity_setter_debug_check(entity);

  m_local_ids[entity.local_offset()] = id;
}

inline bool BulkData::internal_set_parallel_owner_rank_but_not_comm_lists(Entity entity, int in_owner_rank)
{
  TraceIfWatching("stk::mesh::BulkData::set_entity_owner_rank", LOG_ENTITY, entity_key(entity));
  DiagIfWatching(LOG_ENTITY, entity_key(entity), "new owner: " << in_owner_rank);

  entity_setter_debug_check(entity);

  int & nonconst_processor_rank = bucket(entity).m_owner_ranks[bucket_ordinal(entity)];
  if ( in_owner_rank != nonconst_processor_rank ) {
    nonconst_processor_rank = in_owner_rank;

    mark_entity_and_upward_related_entities_as_modified(entity);
    return true;
  }
  return false;
}

inline void BulkData::log_created_parallel_copy(Entity entity)
{
  if (state(entity) == Created) {
    set_state(entity, Modified);
  }
}

inline bool BulkData::is_valid_connectivity(Entity entity, EntityRank rank) const
{
  if (!is_valid(entity)) return false;
  if (bucket_ptr(entity) == NULL) return false;
  internal_check_unpopulated_relations(entity, rank);
  return true;
}

}
}

#endif
