#include <gtest/gtest.h>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp>

#include <stk_unit_test_utils/ioUtils.hpp>

#include "TetFixture.hpp"

class DisconnectMesh : public DGTetFixture {};

void createSerialSubMesh(const stk::mesh::MetaData &oldMeta,
                         stk::mesh::BulkData& oldBulkData,
                         stk::mesh::Selector subMeshSelector,
                         stk::mesh::MetaData &newMeta,
                         stk::mesh::BulkData &newBulkData);

TEST_F(DisconnectMesh, tet)
{
    std::vector<stk::mesh::EntityIdVector> tet_conn = {
            {1, 2, 3, 4}, // id 1
            {2, 3, 4, 5}  // id 2
    };

    std::vector< std::vector<double> > node_coords= {
            {0, 0, 0}, // 1
            {1, 0, 0}, // 2
            {0, 1, 0}, // 3
            {0.5, 0.5, 1.0}, // 6...just kidding, it's 4
            {1.0, 1.0, 1.0}
    };

    setup_mesh(tet_conn, node_coords);

    //////////////////////////////////////////////////////////////////////////////////////

    stk::mesh::MetaData newMetaData(get_meta().spatial_dimension());
    stk::mesh::BulkData newBulkData(newMetaData, get_bulk().parallel());
    createSerialSubMesh(get_meta(), get_bulk(), get_meta().universal_part(), newMetaData, newBulkData);

    stk::unit_test_util::write_mesh_using_stk_io("mike_new.g", newBulkData, newBulkData.parallel());
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void copyPartsToNewMesh(const stk::mesh::MetaData &oldMeta, stk::mesh::MetaData &newMeta)
{
    const stk::mesh::PartVector &allparts = oldMeta.get_mesh_parts();
    for(size_t i = 0; i < allparts.size(); i++)
    {
        stk::mesh::Part *part = NULL;
        if(allparts[i]->topology() != stk::topology::INVALID_TOPOLOGY)
        {
            part = &newMeta.declare_part_with_topology(allparts[i]->name(), allparts[i]->topology());
        }
        else
        {
            part = &newMeta.declare_part(allparts[i]->name());
        }
        if(stk::io::is_part_io_part(*allparts[i]))
        {
            stk::io::put_io_part_attribute(*part);
        }
    }
}

void copyFieldsToNewMesh(const stk::mesh::MetaData &oldMeta, stk::mesh::MetaData &newMeta)
{
    const stk::mesh::FieldVector &fields = oldMeta.get_fields();
    for(size_t i = 0; i < fields.size(); i++)
    {
        stk::mesh::FieldBase* newField = newMeta.declare_field_base(fields[i]->name(),
                                                                    fields[i]->entity_rank(),
                                                                    fields[i]->data_traits(),
                                                                    fields[i]->field_array_rank(),
                                                                    fields[i]->dimension_tags(),
                                                                    fields[i]->number_of_states());

        stk::mesh::Selector selectFieldParts = stk::mesh::selectField(*fields[i]);
        stk::mesh::PartVector oldParts;
        selectFieldParts.get_parts(oldParts);
        if(!oldParts.empty())
        {
            stk::mesh::PartVector newParts(oldParts.size());
            for(size_t k = 0; k < oldParts.size(); k++)
            {
                newParts[k] = &newMeta.get_part(oldParts[k]->mesh_meta_data_ordinal());
            }
            stk::mesh::Selector selectNewParts = stk::mesh::selectUnion(newParts);
            stk::mesh::put_field(*newField, selectNewParts, fields[i]->max_size(fields[i]->entity_rank()));
        }
    }
}

void copySelectedEntitiesToNewMesh(const stk::mesh::BulkData& oldBulkData,
                                   const stk::mesh::Selector subMeshSelector,
                                   std::map<stk::mesh::Entity, stk::mesh::Entity> &oldToNewEntityMap,
                                   stk::mesh::MetaData &newMeta,
                                   stk::mesh::BulkData &newBulkData)
{
    for(stk::mesh::EntityRank rank = stk::topology::NODE_RANK; rank <= stk::topology::ELEMENT_RANK; rank++)
    {
        const stk::mesh::BucketVector &buckets = oldBulkData.get_buckets(rank, subMeshSelector);
        for(size_t i = 0; i < buckets.size(); i++)
        {
            stk::mesh::Bucket &bucket = *buckets[i];
            for(size_t j = 0; j < bucket.size(); j++)
            {
                const stk::mesh::PartVector &oldParts = bucket.supersets();
                stk::mesh::PartVector newParts(oldParts.size(), 0);
                for(size_t k = 0; k < oldParts.size(); k++)
                {
                    newParts[k] = &newMeta.get_part(oldParts[k]->mesh_meta_data_ordinal());
                }
                stk::mesh::Entity newEntity = newBulkData.declare_entity(rank, oldBulkData.identifier(bucket[j]), newParts);
                oldToNewEntityMap[bucket[j]] = newEntity;
            }
        }
    }
}

void copyRelationsToNewMesh(const stk::mesh::BulkData& oldBulkData,
                            const stk::mesh::Selector subMeshSelector,
                            const std::map<stk::mesh::Entity, stk::mesh::Entity> &oldToNewEntityMap,
                            stk::mesh::BulkData &newBulkData)
{
    for(stk::mesh::EntityRank rank = stk::topology::NODE_RANK; rank <= stk::topology::ELEMENT_RANK; rank++)
    {
        const stk::mesh::BucketVector &buckets = oldBulkData.get_buckets(rank, subMeshSelector);
        for(size_t i = 0; i < buckets.size(); i++)
        {
            stk::mesh::Bucket &bucket = *buckets[i];
            for(size_t j = 0; j < bucket.size(); j++)
            {
                stk::mesh::Entity newFromEntity = oldToNewEntityMap.find(bucket[j])->second;
                for(stk::mesh::EntityRank downRank = stk::topology::NODE_RANK; downRank < rank; downRank++)
                {
                    unsigned numConnectedEntities = oldBulkData.num_connectivity(bucket[j], downRank);
                    const stk::mesh::Entity *connectedEntities = oldBulkData.begin(bucket[j], downRank);
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

void copyFieldDataToNewMesh(const stk::mesh::MetaData &oldMeta,
                            const stk::mesh::BulkData& oldBulkData,
                            const stk::mesh::Selector subMeshSelector,
                            const std::map<stk::mesh::Entity, stk::mesh::Entity> &oldToNewEntityMap,
                            stk::mesh::MetaData &newMeta)
{
    const stk::mesh::FieldVector &oldFields = oldMeta.get_fields();
    const stk::mesh::FieldVector &newFields = newMeta.get_fields();
    for(stk::mesh::EntityRank rank = stk::topology::NODE_RANK; rank <= stk::topology::ELEMENT_RANK; rank++)
    {
        const stk::mesh::BucketVector &buckets = oldBulkData.get_buckets(rank, subMeshSelector);
        for(size_t i = 0; i < buckets.size(); i++)
        {
            stk::mesh::Bucket &bucket = *buckets[i];
            for(size_t k = 0; k < oldFields.size(); k++)
            {
                if(bucket.field_data_is_allocated(*oldFields[k]))
                {
                    for(size_t j = 0; j < bucket.size(); j++)
                    {
                        stk::mesh::Entity oldEntity = bucket[j];
                        stk::mesh::Entity newEntity = oldToNewEntityMap.find(oldEntity)->second;
                        void *oldData = stk::mesh::field_data(*oldFields[k], oldEntity);
                        void *newData = stk::mesh::field_data(*newFields[k], newEntity);
                        memcpy(newData, oldData, stk::mesh::field_bytes_per_entity(*oldFields[k], oldEntity));
                    }
                }
            }
        }
    }
}

void createSerialSubMesh(const stk::mesh::MetaData &oldMeta,
                         stk::mesh::BulkData& oldBulkData,
                         stk::mesh::Selector subMeshSelector,
                         stk::mesh::MetaData &newMeta,
                         stk::mesh::BulkData &newBulkData)
{
    copyPartsToNewMesh(oldMeta, newMeta);
    copyFieldsToNewMesh(oldMeta, newMeta);

    newMeta.commit();

    std::map<stk::mesh::Entity, stk::mesh::Entity> oldToNewEntityMap;

    newBulkData.modification_begin();
    copySelectedEntitiesToNewMesh(oldBulkData, subMeshSelector, oldToNewEntityMap, newMeta, newBulkData);
    copyRelationsToNewMesh(oldBulkData, subMeshSelector, oldToNewEntityMap, newBulkData);
    newBulkData.modification_end();

    copyFieldDataToNewMesh(oldMeta, oldBulkData, subMeshSelector, oldToNewEntityMap, newMeta);
}
