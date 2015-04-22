#ifndef STK_ENTITYCOMMLIST_INFO_HPP
#define STK_ENTITYCOMMLIST_INFO_HPP

namespace stk {
namespace mesh {

struct EntityCommListInfo
{
  EntityKey key;
  Entity    entity; // Might be invalid if entity has been deleted.
  int  owner;
  const EntityComm* entity_comm; // Might be NULL if entity has been deleted.
};

inline
bool operator<(const EntityKey& key, const EntityCommListInfo& comm) { return key < comm.key; }

inline
bool operator<(const EntityCommListInfo& comm, const EntityKey& key) { return comm.key < key; }

inline
bool operator<(const EntityCommListInfo& lhs, const EntityCommListInfo& rhs) { return lhs.key < rhs.key; }

inline
bool operator==(const EntityCommListInfo& lhs, const EntityCommListInfo& rhs) { return lhs.key == rhs.key; }

inline
bool operator!=(const EntityCommListInfo& lhs, const EntityCommListInfo& rhs) { return !(lhs == rhs); }

typedef TrackedVectorMetaFunc<EntityCommListInfo, EntityCommTag>::type EntityCommListInfoVector;

struct IsInvalid
{
  bool operator()(const EntityCommListInfo& comm) const
  {
    return comm.key == EntityKey();
  }
};

}
}

#endif
