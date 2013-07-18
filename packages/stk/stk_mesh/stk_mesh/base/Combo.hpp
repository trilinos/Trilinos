/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_mesh_base_Combo_hpp
#define stk_mesh_base_Combo_hpp

#include <stk_mesh/base/EntityKey.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Trace.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/ConnectivityMap.hpp>

#include <stk_util/environment/ReportHandler.hpp>
#include <stk_util/environment/OutputLog.hpp>
#include <stk_util/util/PairIter.hpp>

#include <boost/static_assert.hpp>
#include <boost/range.hpp>
#include <boost/type_traits/is_pod.hpp>

#include <stk_topology/topology.hpp>

#include <utility>
#include <vector>
#include <iosfwd>
#include <string>
#include <algorithm>

#include <boost/functional/hash.hpp>

#ifdef SIERRA_MIGRATION
namespace stk {
namespace mesh {
typedef RelationVector::const_iterator   RelationIterator;
typedef boost::iterator_range<RelationIterator> RelationRange;
}
}
#endif // SIERRA_MIGRATION

#include <stk_mesh/base/Entity.tcc> //only place where this file should be included

#include <stk_mesh/base/Relation.tcc> //only place where this file should be included

#include <stk_mesh/base/BucketConnectivity.tcc> //only place where this file should be included

#include <stk_mesh/base/Bucket.tcc> //only place where this file should be included

///////////////////////////////////////////////////////////////////////////////
// Put methods below that could not otherwise be inlined due to cyclic dependencies between
// Relation/Entity/Bucket
///////////////////////////////////////////////////////////////////////////////

namespace stk {
namespace mesh {

//
// Relation
//

inline
Entity Relation::entity() const
{
  return m_target_entity;
}

#ifdef SIERRA_MIGRATION

inline
Relation::Relation(EntityRank rel_rank, Entity obj, const unsigned relation_type, const unsigned ordinal, const unsigned orient)
  :
      m_raw_relation( Relation::raw_relation_id(rel_rank, ordinal )),
      m_attribute( (relation_type << fmwk_orientation_digits) | orient ),
      m_target_entity(obj)
{
  ThrowAssertMsg( orient <= fmwk_orientation_mask,
      "orientation " << orient << " exceeds maximum allowed value");
}

inline
void Relation::setMeshObj(Entity object, EntityRank object_rank )
{
  m_raw_relation = Relation::raw_relation_id( object_rank, relation_ordinal() );
  m_target_entity = object;
}

#endif

//
// Entity
//

inline
size_t hash_value( Entity entity) {
  return boost::hash_value(entity.local_offset());
}

//
// BucketConnectivity
//

template <EntityRank TargetRank>
template <typename BulkData> // hack to get around dependency
inline
void impl::BucketConnectivity<TargetRank, FIXED_CONNECTIVITY>::end_modification(BulkData* mesh)
{
  //TODO: If bucket is blocked, no longer need to shrink to fit!

  {
    EntityVector temp(m_targets.begin(), m_targets.end());
    m_targets.swap(temp);
  }

  {
    PermutationVector temp(m_permutations.begin(), m_permutations.end());
    m_permutations.swap(temp);
  }

  invariant_check_helper(mesh);
}

template <EntityRank TargetRank>
template <typename BulkData>
inline
void impl::BucketConnectivity<TargetRank, DYNAMIC_CONNECTIVITY>::end_modification(BulkData* mesh)
{
  if (m_active) {
    resize_and_order_by_index();

    {
      UInt32Vector temp(m_indices.begin(), m_indices.end());
      m_indices.swap(temp);
    }

    {
      UInt16Vector temp(m_num_connectivities.begin(), m_num_connectivities.end());
      m_num_connectivities.swap(temp);
    }

    invariant_check_helper(mesh);
  }
}

template <typename BulkData>
inline
bool impl::LowerConnectivitityRankSensitiveCompare<BulkData>::operator()(Entity first_entity, ConnectivityOrdinal first_ordinal,
                                                                         Entity second_entity, ConnectivityOrdinal second_ordinal) const
{
  const EntityRank first_rank = m_mesh.entity_rank(first_entity);
  const EntityRank second_rank = m_mesh.entity_rank(second_entity);

  return (first_rank < second_rank)
         || ((first_rank == second_rank) && (first_ordinal < second_ordinal));
}

template <typename BulkData>
inline
bool impl::HigherConnectivityRankSensitiveCompare<BulkData>::operator()(Entity first_entity, ConnectivityOrdinal first_ordinal, Entity second_entity, ConnectivityOrdinal second_ordinal) const
{
  const EntityRank first_rank = m_mesh.entity_rank(first_entity);
  const EntityRank second_rank = m_mesh.entity_rank(second_entity);

  if (first_rank < second_rank) {
    return true;
  }
  if (first_rank > second_rank) {
    return false;
  }
  // Needs to match LessRelation in BulkData.hpp
  return std::make_pair(first_ordinal,  first_entity.is_local_offset_valid() ?  first_entity.local_offset()  : Entity::MaxEntity) <
         std::make_pair(second_ordinal, second_entity.is_local_offset_valid() ? second_entity.local_offset() : Entity::MaxEntity);
}

} // namespace mesh
} // namespace stk

#endif /* stk_mesh_Combo_hpp */
