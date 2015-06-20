#ifndef STK_ENTITYLESS_HPP
#define STK_ENTITYLESS_HPP

namespace stk {
namespace mesh {

struct EntityLess {
  inline EntityLess(const BulkData& mesh);
  /** \brief  Comparison operator */
  inline bool operator()(const Entity lhs, const Entity rhs) const;
  inline bool operator()(const Entity lhs, const EntityKey & rhs) const;
  inline bool operator()( const EntityProc & lhs, const EntityProc & rhs) const;
  inline bool operator()( const EntityProc & lhs, const Entity rhs) const;
  inline bool operator()( const EntityProc & lhs, const EntityKey & rhs) const;
  inline EntityLess& operator=(const EntityLess& rhs);
  const BulkData* m_mesh;
}; //struct EntityLess

}
}

#endif
