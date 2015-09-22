#ifndef __ENTITY_KEY_HASH__
#define __ENTITY_KEY_HASH__

#include <unordered_map>
#include <unordered_set>
#include <stk_mesh/base/EntityKey.hpp>

namespace std
{
    template<>
    struct hash<stk::mesh::EntityKey>
    {
        const size_t operator()(const stk::mesh::EntityKey& entityKey) const
        {
            return hash<unsigned long long>()(entityKey.m_value);
        }
        const size_t operator()(const stk::mesh::EntityKey& entityKey, const stk::mesh::Entity& entity) const
        {
            return hash<unsigned long long>()(entityKey.m_value)^hash<unsigned long long>()(entity.m_value);
        }
        const size_t operator()(const stk::mesh::EntityKey& entityKey1, const stk::mesh::EntityKey& entityKey2) const
        {
            return hash<unsigned long long>()(entityKey1.m_value)^hash<unsigned long long>()(entityKey2.m_value);
        }
    };

    template<>
    struct hash<stk::mesh::Entity>
    {
        const size_t operator()(const stk::mesh::Entity& entity) const
        {
            return hash<unsigned long long>()(entity.m_value);
        }
        const size_t operator()(const stk::mesh::Entity& entity, const stk::mesh::EntityKey& entityKey) const
        {
            return hash<unsigned long long>()(entity.m_value)^hash<unsigned long long>()(entityKey.m_value);
        }
        const size_t operator()(const stk::mesh::Entity& entity1, const stk::mesh::Entity& entity2) const
        {
            return hash<unsigned long long>()(entity1.m_value)^hash<unsigned long long>()(entity2.m_value);
        }
    };
}

#endif // __ENTITY_KEY_HASH__
