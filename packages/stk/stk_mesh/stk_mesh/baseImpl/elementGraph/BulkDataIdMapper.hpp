#ifndef BULK_DATA_ID_MAPPER_HPP
#define BULK_DATA_ID_MAPPER_HPP

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include "../../base/GetEntities.hpp"
#include "ElemElemGraphImpl.hpp"
#include "ElemGraphCoincidentElems.hpp"

namespace stk
{
namespace mesh
{
namespace impl
{

class ElementLocalIdMapper
{
public:
    ElementLocalIdMapper()
    : localIdToElement(), elementToLocalId()
    {
    }
    void initialize(const stk::mesh::BulkData &bulk)
    {
        unsigned numElems = stk::mesh::count_selected_entities(bulk.mesh_meta_data().locally_owned_part(), bulk.buckets(stk::topology::ELEM_RANK));
        localIdToElement.resize(numElems, Entity());
        elementToLocalId.resize(stk::mesh::get_num_entities(bulk)+1, impl::INVALID_LOCAL_ID);

        const stk::mesh::BucketVector & elemBuckets = bulk.get_buckets(stk::topology::ELEM_RANK, bulk.mesh_meta_data().locally_owned_part());
        impl::LocalId localId = 0;
        for(size_t i=0; i<elemBuckets.size(); ++i)
        {
            const stk::mesh::Bucket& bucket = *elemBuckets[i];
            for(size_t j=0; j<bucket.size(); ++j)
            {
                add_new_entity_with_local_id(bucket[j], localId);
                localId++;
            }
        }
    }
    void clear()
    {
        localIdToElement.clear();
        elementToLocalId.clear();
    }
    stk::mesh::Entity local_to_entity(stk::mesh::impl::LocalId local) const
    {
        return localIdToElement[local];
    }

    bool does_entity_have_local_id(stk::mesh::Entity entity) const
    {
        return (entity.local_offset() < elementToLocalId.size() && elementToLocalId[entity.local_offset()] != stk::mesh::impl::INVALID_LOCAL_ID);
    }

    stk::mesh::impl::LocalId entity_to_local(stk::mesh::Entity entity) const
    {
        ThrowAssert(elementToLocalId.size() > entity.local_offset());

        return elementToLocalId[entity.local_offset()];
    }

    void make_space_for_local_id(stk::mesh::impl::LocalId localId)
    {
        if(static_cast<size_t>(localId) >= localIdToElement.size())
            localIdToElement.resize(localId + 1);
    }

    void make_space_for_new_elements(const stk::mesh::EntityVector & elements)
    {
        stk::mesh::Entity maxEntity;
        maxEntity.set_local_offset(0);
        for(stk::mesh::Entity elem : elements)
            maxEntity.set_local_offset(std::max(maxEntity.local_offset(), elem.local_offset()));

        if(maxEntity.local_offset() >= elementToLocalId.size())
            elementToLocalId.resize(maxEntity.local_offset() + 1);
    }

    void create_local_ids_for_elements(const stk::mesh::EntityVector & elements)
    {
        make_space_for_new_elements(elements);
        localIdToElement.resize(elements.size());
        for(size_t localId = 0; localId < elements.size(); localId++)
            add_new_entity_with_local_id(elements[localId], localId);
    }

    void add_new_entity_with_local_id(stk::mesh::Entity elem, impl::LocalId id)
    {
        localIdToElement[id] = elem;
        elementToLocalId[elem.local_offset()] = id;
    }

    void clear_local_id_for_elem(stk::mesh::Entity elem)
    {
        localIdToElement[entity_to_local(elem)] = stk::mesh::Entity::InvalidEntity;
        elementToLocalId[elem.local_offset()] = impl::INVALID_LOCAL_ID;
    }
private:
    stk::mesh::EntityVector localIdToElement;
    std::vector<stk::mesh::impl::LocalId> elementToLocalId;
};

class BulkDataIdMapper : public IdMapper
{
public:
    BulkDataIdMapper(const stk::mesh::BulkData &b,
                     const ElementLocalIdMapper & l)
    : bulk(b), localIdMaps(l)
    {
    }

    virtual stk::mesh::EntityId localToGlobal(stk::mesh::impl::LocalId localId) const
    {
        if(localId < 0)
            return -localId;
        else
            return bulk.identifier(localIdMaps.local_to_entity(localId));
    }
    virtual stk::mesh::impl::LocalId globalToLocal(stk::mesh::EntityId global) const
    {
        stk::mesh::Entity elem = bulk.get_entity(stk::mesh::EntityKey(stk::topology::ELEM_RANK, global));
        return localIdMaps.entity_to_local(elem);
    }

private:
    const stk::mesh::BulkData &bulk;
    const ElementLocalIdMapper & localIdMaps;
};

}
}
}

#endif
