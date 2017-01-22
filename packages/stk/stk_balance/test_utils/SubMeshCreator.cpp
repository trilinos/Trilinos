#include <test_utils/SubMeshCreator.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/FieldTraits.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/SkinMesh.hpp>
#include <stk_io/IossBridge.hpp>
#include <stk_mesh/baseImpl/MeshImplUtils.hpp>

namespace stk {
namespace balance {
namespace internal {

class SubMeshCreator
{
public:

    SubMeshCreator(stk::mesh::BulkData &bulk, stk::mesh::Selector selector) :
            oldMeta(bulk.mesh_meta_data()),
            oldBulk(bulk),
            subMeshSelector(selector),
            all_parts_ordinal_mapping()
    {
    }

    void create_new_sub_mesh(stk::mesh::BulkData &newBulk)
    {
        stk::mesh::MetaData &newMeta = newBulk.mesh_meta_data();
        copy_parts_to_submesh(newMeta);
        all_parts_ordinal_mapping = get_part_ordinal_mapping(newMeta);
        copy_fields_to_submesh(newMeta);

        newMeta.commit();

        std::map<stk::mesh::Entity, stk::mesh::Entity> oldToNewEntityMap;
        populate_new_bulk_data(newBulk, oldToNewEntityMap);

        copy_field_data_to_submesh(oldToNewEntityMap, newMeta);
    }

private:
    void copy_parts_to_submesh(stk::mesh::MetaData& newMeta)
    {
        const stk::mesh::PartVector &all_parts = oldMeta.get_parts();
        for(size_t i = 0; i < all_parts.size(); i++)
        {
            if(!stk::mesh::is_auto_declared_part(*all_parts[i]))
            {
                stk::topology topology = all_parts[i]->topology();
                if(all_parts[i]->topology()==stk::topology::INVALID_TOPOLOGY && all_parts[i]->primary_entity_rank()==stk::topology::NODE_RANK)
                    topology = stk::topology::NODE;
                create_part(newMeta, all_parts[i]->name(), topology);
            }
        }
    }

    void create_part(stk::mesh::MetaData& newMeta, const std::string& part_name, stk::topology topology)
    {
        if(topology != stk::topology::INVALID_TOPOLOGY)
            create_part_with_valid_topology(newMeta, part_name, topology);
        else
            create_part_with_invalid_topology(newMeta, part_name);
    }

    void create_part_with_valid_topology(stk::mesh::MetaData& newMeta, const std::string& part_name, stk::topology topology)
    {
        stk::mesh::Part *part = &newMeta.declare_part_with_topology(part_name, topology);
        stk::io::put_io_part_attribute(*part);
    }

    void create_part_with_invalid_topology(stk::mesh::MetaData& newMeta, const std::string& part_name)
    {
        stk::mesh::Part *part = &newMeta.declare_part(part_name);
        stk::io::put_io_part_attribute(*part);
    }

    std::vector<unsigned> get_part_ordinal_mapping(const stk::mesh::MetaData& newMeta)
    {
        const stk::mesh::PartVector &all_parts = oldMeta.get_parts();
        return map_parts_to_new_meta_data(all_parts, newMeta);
    }

    std::vector<unsigned> map_parts_to_new_meta_data(const stk::mesh::PartVector &all_parts, const stk::mesh::MetaData& newMeta)
    {
        std::vector<unsigned> allPartsMapping(all_parts.size(), 0);
        for(size_t i = 0; i < all_parts.size(); i++)
            insert_part_ordinal_mapping_if_part_exists(*all_parts[i], newMeta, allPartsMapping);
        return allPartsMapping;
    }

    void insert_part_ordinal_mapping_if_part_exists(const stk::mesh::Part& old_part,
                                                    const stk::mesh::MetaData& newMeta,
                                                    std::vector<unsigned> &allPartsMapping)
    {
        stk::mesh::Part* part = newMeta.get_part(old_part.name());
        if(part != nullptr)
            allPartsMapping[old_part.mesh_meta_data_ordinal()] = part->mesh_meta_data_ordinal();
    }

    stk::mesh::PartVector get_parts_from_field(stk::mesh::FieldBase* field)
    {
        stk::mesh::Selector selectFieldParts = stk::mesh::selectField(*field);
        stk::mesh::PartVector oldParts;
        selectFieldParts.get_parts(oldParts);
        return oldParts;
    }
    
    stk::mesh::PartVector map_parts(const stk::mesh::MetaData& newMetaData, const stk::mesh::PartVector& oldParts)
    {
        stk::mesh::PartVector newParts;
        for(stk::mesh::Part* oldPart : oldParts )
            newParts.push_back(&newMetaData.get_part(all_parts_ordinal_mapping[oldPart->mesh_meta_data_ordinal()]));
        return newParts;
    }
    
    void put_field_on_new_parts(stk::mesh::MetaData& newMetaData, const stk::mesh::PartVector& oldParts, stk::mesh::FieldBase& newField)
    {
        if(!oldParts.empty())
        {
            stk::mesh::PartVector newParts = map_parts(newMetaData, oldParts);
            stk::mesh::Selector selectNewParts = stk::mesh::selectUnion(newParts);
            stk::mesh::put_field(newField, selectNewParts, newMetaData.spatial_dimension());
        }
    }

    stk::mesh::FieldBase* declare_new_field(stk::mesh::FieldBase &oldField, stk::mesh::MetaData& newMetaData)
    {
        return newMetaData.declare_field_base(oldField.name(), oldField.entity_rank(), oldField.data_traits(),
            oldField.field_array_rank(), oldField.dimension_tags(), oldField.number_of_states());
    }

    void copy_fields_to_submesh(stk::mesh::MetaData& newMetaData)
    {
        for(stk::mesh::FieldBase* oldField : oldMeta.get_fields() )
        {
            stk::mesh::FieldBase* newField = declare_new_field(*oldField, newMetaData);
            stk::mesh::PartVector oldParts = get_parts_from_field(oldField);
            put_field_on_new_parts(newMetaData, oldParts, *newField);
        }
    }

    void populate_new_bulk_data(stk::mesh::BulkData& new_bulk_data,
                                std::map<stk::mesh::Entity, stk::mesh::Entity>& oldToNewEntityMap)
    {
        new_bulk_data.modification_begin();
        copy_selected_entities_to_new_mesh(oldToNewEntityMap, new_bulk_data);

        copy_relations_to_new_mesh(oldToNewEntityMap, new_bulk_data);
        copy_selected_faces_to_new_mesh(oldToNewEntityMap, new_bulk_data);

        new_bulk_data.modification_end();
    }


    void copy_selected_faces_to_new_mesh(std::map<stk::mesh::Entity, stk::mesh::Entity> &oldToNewEntityMap,
                                         stk::mesh::BulkData &newBulkData)
    {
        const stk::mesh::BucketVector &buckets = oldBulk.get_buckets(stk::topology::FACE_RANK, subMeshSelector);
        for(size_t i = 0; i < buckets.size(); i++)
        {
            stk::mesh::Bucket &bucket = *buckets[i];
            for(size_t j = 0; j < bucket.size(); j++)
            {
                const stk::mesh::Entity oldSide = bucket[j];
                const stk::mesh::PartVector &oldParts = bucket.supersets();
                stk::mesh::PartVector newParts;
                for(size_t k = 0; k < oldParts.size(); k++)
                {
                    if (!stk::mesh::is_auto_declared_part(*oldParts[k]))
                    {
                        newParts.push_back( &newBulkData.mesh_meta_data().get_part(oldParts[k]->mesh_meta_data_ordinal()) );
                    }
                }

                const stk::mesh::Entity * connected_elements = oldBulk.begin_elements(oldSide);
                const stk::mesh::Entity oldFirstElement = connected_elements[0];
                ThrowRequire(oldBulk.is_valid(oldFirstElement));

                stk::mesh::ConnectivityOrdinal faceOrdinal = stk::mesh::impl::get_ordinal_for_element_side_pair(oldBulk, oldFirstElement, oldSide);
                stk::mesh::Entity newFirstElement = newBulkData.get_entity(oldBulk.entity_key(oldFirstElement));

                stk::mesh::Entity newSide = newBulkData.declare_element_side(newFirstElement, faceOrdinal, newParts);
                oldToNewEntityMap[oldSide] = newSide;
            }
        }
    }

    void copy_selected_entities_to_new_mesh(std::map<stk::mesh::Entity, stk::mesh::Entity> &oldToNewEntityMap,
                                            stk::mesh::BulkData &newBulkData)
    {
        for(stk::mesh::EntityRank rank = stk::topology::NODE_RANK; rank <= stk::topology::ELEMENT_RANK; rank++)
        {
            // Skip over all face relations since those will be managed by declare_element_side() later on
            if (rank != stk::topology::FACE_RANK)
            {
                const stk::mesh::BucketVector &buckets = oldBulk.get_buckets(rank, subMeshSelector);
                for(stk::mesh::Bucket *bucket : buckets)
                {
                    if(bucket->owned() || bucket->shared())
                    {
                        for(stk::mesh::Entity oldEntity : *bucket)
                        {
                            copy_entity_to_submesh(oldEntity, bucket->supersets(), oldToNewEntityMap, newBulkData);
                        }
                    }
                }
            }
        }
    }

    void copy_relations_to_new_mesh(std::map<stk::mesh::Entity, stk::mesh::Entity> &oldToNewEntityMap,
                                    stk::mesh::BulkData &newBulkData)
    {
        for(stk::mesh::EntityRank rank = stk::topology::NODE_RANK; rank <= stk::topology::ELEMENT_RANK; rank++)
        {
            // Skip over all face relations since those will be managed by declare_element_side() later on
            if (rank != stk::topology::FACE_RANK)
            {
                const stk::mesh::BucketVector &buckets = oldBulk.get_buckets(rank, subMeshSelector);
                for(size_t i = 0; i < buckets.size(); i++)
                {
                    stk::mesh::Bucket &bucket = *buckets[i];
                    for(size_t j = 0; j < bucket.size(); j++)
                    {
                        stk::mesh::Entity newFromEntity = oldToNewEntityMap.find(bucket[j])->second;
                        for(stk::mesh::EntityRank downRank = stk::topology::NODE_RANK; downRank < rank; downRank++)
                        {
                            if (downRank != stk::topology::FACE_RANK)
                            {
                                unsigned numConnectedEntities = oldBulk.num_connectivity(bucket[j], downRank);
                                const stk::mesh::Entity *connectedEntities = oldBulk.begin(bucket[j], downRank);
                                for(unsigned connectOrder = 0; connectOrder < numConnectedEntities; connectOrder++)
                                {
                                    stk::mesh::Entity newToEntity = oldToNewEntityMap.find(connectedEntities[connectOrder])->second;
                                    newBulkData.declare_relation(newFromEntity, newToEntity, connectOrder);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    void copy_entity_to_submesh(stk::mesh::Entity oldEntity,
                                const stk::mesh::PartVector &oldParts,
                                std::map<stk::mesh::Entity, stk::mesh::Entity> &oldToNewEntityMap,
                                stk::mesh::BulkData &newBulkData)
    {
        stk::mesh::PartVector newParts = get_new_parts_for_entity(oldParts, newBulkData.mesh_meta_data());
        stk::mesh::Entity newEntity = newBulkData.declare_entity(oldBulk.entity_rank(oldEntity),
                                                                 oldBulk.identifier(oldEntity),
                                                                 newParts);
        oldToNewEntityMap[oldEntity] = newEntity;
        //copy_relations(oldEntity, oldBulk.entity_rank(oldEntity), oldToNewEntityMap, newBulkData, newEntity);
    }

    stk::mesh::PartVector get_new_parts_for_entity(const stk::mesh::PartVector &oldParts, const stk::mesh::MetaData &newMeta)
    {
        stk::mesh::PartVector newParts(oldParts.size(), 0);
        for(size_t k = 0; k < oldParts.size(); k++)
            newParts[k] = &newMeta.get_part(all_parts_ordinal_mapping[oldParts[k]->mesh_meta_data_ordinal()]);
        return newParts;
    }

    void copy_relations(stk::mesh::Entity oldEntity,
                        stk::mesh::EntityRank rank,
                        const std::map<stk::mesh::Entity, stk::mesh::Entity> &oldToNewEntityMap,
                        stk::mesh::BulkData &newBulkData,
                        stk::mesh::Entity newEntity)
    {
        for(stk::mesh::EntityRank downRank = stk::topology::NODE_RANK; downRank < rank; downRank++)
        {
            unsigned numConnectedEntities = oldBulk.num_connectivity(oldEntity, downRank);
            const stk::mesh::Entity *connectedEntities = oldBulk.begin(oldEntity, downRank);
            for(unsigned connectOrder = 0; connectOrder < numConnectedEntities; connectOrder++)
            {
                stk::mesh::Entity newToEntity = oldToNewEntityMap.find(connectedEntities[connectOrder])->second;
                newBulkData.declare_relation(newEntity, newToEntity, connectOrder);
            }
        }
    }

    void copy_field_data_to_submesh(const std::map<stk::mesh::Entity, stk::mesh::Entity> &oldToNewEntityMap,
                                    stk::mesh::MetaData &newMeta)
    {
        const stk::mesh::FieldVector &fields = oldMeta.get_fields();
        const stk::mesh::FieldVector &newFields = newMeta.get_fields();
        for(stk::mesh::EntityRank rank = stk::topology::NODE_RANK; rank <= stk::topology::ELEMENT_RANK; rank++)
        {
            const stk::mesh::BucketVector &buckets = oldBulk.get_buckets(rank, subMeshSelector);
            for(size_t i = 0; i < buckets.size(); i++)
            {
                stk::mesh::Bucket &bucket = *buckets[i];
                if(bucket.owned() || bucket.shared())
                {
                    for(size_t k = 0; k < fields.size(); k++)
                    {
                        if(bucket.field_data_is_allocated(*fields[k]))
                        {
                            for(size_t j = 0; j < bucket.size(); j++)
                            {
                                stk::mesh::Entity oldEntity = bucket[j];
                                stk::mesh::Entity newEntity = oldToNewEntityMap.find(oldEntity)->second;
                                void *oldData = stk::mesh::field_data(*fields[k], oldEntity);
                                void *newData = stk::mesh::field_data(*newFields[k], newEntity);
                                memcpy(newData, oldData, stk::mesh::field_bytes_per_entity(*fields[k], oldEntity));
                            }
                        }
                    }
                }
            }
        }
    }

private:
    stk::mesh::MetaData &oldMeta;
    stk::mesh::BulkData &oldBulk;
    stk::mesh::Selector subMeshSelector;
    std::vector<unsigned> all_parts_ordinal_mapping;
};

void createNewSubMesh(stk::mesh::MetaData &oldMeta,
                      stk::mesh::BulkData &oldBulk,
                      stk::mesh::Selector subMeshSelector,
                      stk::mesh::MetaData &newMeta,
                      stk::mesh::BulkData &newBulk)
{
    SubMeshCreator subMeshCreator(oldBulk, subMeshSelector);
    subMeshCreator.create_new_sub_mesh(newBulk);
}

}
}
}
