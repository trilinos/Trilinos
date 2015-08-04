#ifndef STK_ENTITYLESS_HPP
#define STK_ENTITYLESS_HPP

namespace stk {
namespace mesh {

class EntityLess {
public:
#ifdef SIERRA_MIGRATION
  EntityLess(const BulkData& mesh);
#else
  inline EntityLess(const BulkData& mesh);
#endif
  /** \brief  Comparison operator */
  inline bool operator()(const Entity lhs, const Entity rhs) const;
  inline bool operator()(const Entity lhs, const EntityKey & rhs) const;
  inline bool operator()( const EntityProc & lhs, const EntityProc & rhs) const;
  inline bool operator()( const EntityProc & lhs, const Entity rhs) const;
  inline bool operator()( const EntityProc & lhs, const EntityKey & rhs) const;
  inline EntityLess& operator=(const EntityLess& rhs);
private:
  const BulkData* m_mesh;
#ifdef SIERRA_MIGRATION
  const bool m_shouldSortFacesByNodeIds;
  const EntityRank m_sideRank;
#endif
}; //struct EntityLess

}
}

#endif
