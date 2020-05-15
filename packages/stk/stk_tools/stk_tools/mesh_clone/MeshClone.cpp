#include "MeshClone.hpp"
#include "MeshCloneIo.hpp"
#include "MeshCloneUtils.hpp"
#include <stk_mesh/base/DestroyElements.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include "stk_mesh/base/Part.hpp"
#include "stk_mesh/base/Field.hpp"
#include "stk_mesh/base/FieldBase.hpp"
#include "stk_topology/topology.hpp"
#include "stk_util/util/ReportHandler.hpp"

namespace stk {
namespace mesh {

//Ugly comment admitting failure to write good code: this BulkData friend
//hopefully goes away after everyone uses declare_element_side that doesn't take an id
stk::mesh::Entity clone_element_side(stk::mesh::BulkData &bulk,
                                     stk::mesh::EntityId id,
                                     stk::mesh::Entity elem,
                                     stk::mesh::ConnectivityOrdinal ord,
                                     const stk::mesh::PartVector &parts)
{
    return bulk.declare_element_side_with_id(id, elem, ord, parts);
}

}}

namespace stk {
namespace tools {

stk::mesh::Part *create_new_part(const stk::mesh::Part &oldPart, stk::mesh::MetaData &newMeta)
{
    if(oldPart.topology() != stk::topology::INVALID_TOPOLOGY)
        return &newMeta.declare_part_with_topology(oldPart.name(), oldPart.topology());
    else
    {
        if(oldPart.primary_entity_rank() == stk::topology::INVALID_RANK)
            return &newMeta.declare_part(oldPart.name());
        else
            return &newMeta.declare_part(oldPart.name(), oldPart.primary_entity_rank(), oldPart.force_no_induce());
    }
}

void clone_part_to_other_meta(const stk::mesh::Part &oldPart, stk::mesh::MetaData &newMeta)
{
    stk::mesh::Part *newPart = create_new_part(oldPart, newMeta);
    newMeta.set_part_id(*newPart, oldPart.id());
}

void copy_parts(const stk::mesh::MetaData &oldMeta, stk::mesh::MetaData &newMeta)
{
    const stk::mesh::PartVector &allParts = oldMeta.get_mesh_parts();
    for(size_t i = 0; i < allParts.size(); i++)
        clone_part_to_other_meta(*allParts[i], newMeta);
    for(size_t i = 0; i < allParts.size(); i++)
        copy_part_supersets(*allParts[i], *get_corresponding_part(newMeta, allParts[i]), oldMeta, newMeta);
}


stk::mesh::FieldBase* clone_field(const stk::mesh::FieldBase& field, stk::mesh::MetaData& newMeta)
{
    return newMeta.declare_field_base(field.name(),
                                      field.entity_rank(),
                                      field.data_traits(),
                                      field.field_array_rank(),
                                      field.dimension_tags(),
                                      field.number_of_states());
}

void copy_field_restrictions(const stk::mesh::FieldBase& field, stk::mesh::MetaData& newMeta, stk::mesh::FieldBase* newField)
{
    const stk::mesh::FieldRestrictionVector& oldRestrictions = field.restrictions();
    for(const stk::mesh::FieldRestriction& res : oldRestrictions)
    {
        stk::mesh::Selector selectNewParts = res.selector().clone_for_different_mesh(newMeta);
        newMeta.declare_field_restriction(*newField,
                                          selectNewParts,
                                          res.num_scalars_per_entity(),
                                          res.dimension(),
                                          field.get_initial_value());
    }
}

stk::mesh::FieldBase * get_corresponding_field(stk::mesh::MetaData &newMeta, const stk::mesh::FieldBase &oldField)
{
    return newMeta.get_fields()[oldField.mesh_meta_data_ordinal()];
}



void copy_fields(const stk::mesh::MetaData &oldMeta, stk::mesh::MetaData &newMeta)
{
    const stk::mesh::FieldVector &fields = oldMeta.get_fields();
    for(size_t i = 0; i < fields.size(); i++)
    {
        if(fields[i]->state() == stk::mesh::StateNone)
        {
            stk::mesh::FieldBase* newField = clone_field(*fields[i], newMeta);
            copy_field_restrictions(*fields[i], newMeta, newField);
        }
    }
    newMeta.set_coordinate_field(get_corresponding_field(newMeta, *oldMeta.coordinate_field()));
}

void copy_surface_to_block_mapping(const stk::mesh::MetaData &oldMeta, stk::mesh::MetaData &newMeta)
{
    std::vector<const stk::mesh::Part *> oldSurfacesInMap = oldMeta.get_surfaces_in_surface_to_block_map();
    for(size_t i=0;i<oldSurfacesInMap.size();++i)
    {
        stk::mesh::Part* surfaceNew = get_corresponding_part(newMeta, oldSurfacesInMap[i]);
        ThrowRequireWithSierraHelpMsg(surfaceNew != nullptr);

        std::vector<const stk::mesh::Part*> oldBlocks = oldMeta.get_blocks_touching_surface(oldSurfacesInMap[i]);
        std::vector<const stk::mesh::Part*> newBlocks(oldBlocks.size());
        for(size_t ii=0;ii<oldBlocks.size();++ii)
        {
            stk::mesh::Part* newBlock = get_corresponding_part(newMeta, oldBlocks[ii]);
            ThrowRequireWithSierraHelpMsg(newBlock != nullptr);
            newBlocks[ii] = newBlock;
        }
        newMeta.set_surface_to_block_mapping(surfaceNew, newBlocks);
    }
}

void copy_meta(const stk::mesh::MetaData &inputMeta, stk::mesh::MetaData &outputMeta)
{
    // Query the coordinate field, to figure out the final name (if none set by the user)
    inputMeta.coordinate_field();

    outputMeta.initialize(inputMeta.spatial_dimension(),
                          inputMeta.entity_rank_names(),
                          inputMeta.coordinate_field_name());
    copy_parts(inputMeta, outputMeta);
    copy_fields(inputMeta, outputMeta);
    copy_surface_to_block_mapping(inputMeta, outputMeta);
}

void copy_relations(const stk::mesh::BulkData& oldBulk,
                    stk::mesh::Selector oldSelector,
                    stk::mesh::EntityRank rank,
                    stk::mesh::Entity oldEntity,
                    stk::mesh::Entity newEntity,
                    stk::mesh::BulkData& outputBulk)
{
    for(stk::mesh::EntityRank relationRank = stk::topology::NODE_RANK; relationRank < outputBulk.mesh_meta_data().entity_rank_count(); relationRank++)
    {
        unsigned numConnected = oldBulk.num_connectivity(oldEntity, relationRank);
        const stk::mesh::Entity* connected = oldBulk.begin(oldEntity, relationRank);
        const stk::mesh::ConnectivityOrdinal *oldOrds = oldBulk.begin_ordinals(oldEntity, relationRank);
        const stk::mesh::Permutation *oldPerms = oldBulk.begin_permutations(oldEntity, relationRank);
        for(unsigned conIndex = 0; conIndex < numConnected; conIndex++)
        {
            if(oldSelector(oldBulk.bucket(connected[conIndex])))
            {
                stk::mesh::Entity to = outputBulk.get_entity(oldBulk.entity_key(connected[conIndex]));
                if(outputBulk.is_valid(to))
                {
                    if(oldPerms != nullptr)
                    {
                        if(rank > relationRank)
                            outputBulk.declare_relation(newEntity, to, oldOrds[conIndex], oldPerms[conIndex]);
                        else
                            outputBulk.declare_relation(to, newEntity, oldOrds[conIndex], oldPerms[conIndex]);
                    }
                    else
                    {
                        if(rank > relationRank)
                            outputBulk.declare_relation(newEntity, to, oldOrds[conIndex]);
                        else
                            outputBulk.declare_relation(to, newEntity, oldOrds[conIndex]);
                    }
                }
            }
        }
    }
}

bool is_comm_self(stk::mesh::BulkData &bulk)
{
    return bulk.parallel() == MPI_COMM_SELF;
}

bool should_copy_bucket(const stk::mesh::Bucket* bucket, stk::mesh::BulkData& outputBulk)
{
    return bucket->owned() || bucket->shared() || is_comm_self(outputBulk);
}

void copy_relations(const stk::mesh::BulkData& oldBulk,
                    stk::mesh::Selector oldSelector,
                    stk::mesh::EntityRank rank,
                    stk::mesh::EntityRank relationRank,
                    stk::mesh::BulkData& outputBulk)
{

    for(const stk::mesh::Bucket* bucket : oldBulk.get_buckets(rank, oldSelector))
    {
        if(should_copy_bucket(bucket, outputBulk))
        {
            for(stk::mesh::Entity oldEntity : *bucket)
            {
                unsigned numConnected = oldBulk.num_connectivity(oldEntity, relationRank);
                const stk::mesh::Entity* connected = oldBulk.begin(oldEntity, relationRank);
                const stk::mesh::ConnectivityOrdinal *oldOrds = oldBulk.begin_ordinals(oldEntity, relationRank);
                const stk::mesh::Permutation *oldPerms = oldBulk.begin_permutations(oldEntity, relationRank);
                for(unsigned conIndex = 0; conIndex < numConnected; conIndex++)
                {
                    if(oldSelector(oldBulk.bucket(connected[conIndex])))
                    {
                        stk::mesh::Entity to = outputBulk.get_entity(oldBulk.entity_key(connected[conIndex]));
                        stk::mesh::Entity newEntity = outputBulk.get_entity(oldBulk.entity_key(oldEntity));
                        if(outputBulk.is_valid(to) && outputBulk.is_valid(newEntity))
                        {
                            if(oldPerms != nullptr)
                                outputBulk.declare_relation(newEntity, to, oldOrds[conIndex], oldPerms[conIndex]);
                            else
                                outputBulk.declare_relation(newEntity, to, oldOrds[conIndex]);
                        }
                    }
                }
            }
        }
    }
}


void copy_field_data(const stk::mesh::BulkData& oldBulk,
                     stk::mesh::EntityRank rank,
                     stk::mesh::Entity oldEntity,
                     stk::mesh::Entity newEntity,
                     stk::mesh::BulkData& newBulk)
{
    const stk::mesh::FieldVector &oldFields = oldBulk.mesh_meta_data().get_fields(rank);
    const stk::mesh::FieldVector &newFields = newBulk.mesh_meta_data().get_fields(rank);

    for(unsigned i=0; i<oldFields.size(); ++i)
    {
        unsigned oldBytesPerEntity = stk::mesh::field_bytes_per_entity(*oldFields[i], oldEntity);

        unsigned char* oldData = static_cast<unsigned char*>(stk::mesh::field_data(*oldFields[i], oldEntity));
        unsigned char* newData = static_cast<unsigned char*>(stk::mesh::field_data(*newFields[i], newEntity));

        for(unsigned j=0; j<oldBytesPerEntity; ++j)
            newData[j] = oldData[j];
    }
}

stk::mesh::PartVector get_corresponding_part_vector(stk::mesh::MetaData &newMeta, const stk::mesh::PartVector &oldParts)
{
    stk::mesh::PartVector newParts;
    newParts.reserve(oldParts.size());

    for(size_t p=0; p < oldParts.size(); p++)
    {
        stk::mesh::Part* part = get_corresponding_part(newMeta, oldParts[p]);
        if(part!=nullptr)
            newParts.push_back(part);
    }

    return newParts;
}

void remove_shared_part(stk::mesh::PartVector &oldParts, stk::mesh::Part &oldSharedPart)
{
    std::remove(oldParts.begin(), oldParts.end(), &oldSharedPart);
}

stk::mesh::PartVector get_new_parts(const stk::mesh::Bucket *bucket, stk::mesh::BulkData &outputBulk)
{
    const stk::mesh::PartVector &oldParts = bucket->supersets();
    stk::mesh::PartVector newParts = get_corresponding_part_vector(outputBulk.mesh_meta_data(), oldParts);
    if(is_comm_self(outputBulk) && bucket->shared())
        remove_shared_part(newParts, outputBulk.mesh_meta_data().globally_shared_part());
    std::remove(newParts.begin(), newParts.end(), &outputBulk.mesh_meta_data().universal_part());
    std::remove(newParts.begin(), newParts.end(), &outputBulk.mesh_meta_data().globally_shared_part());
    std::remove(newParts.begin(), newParts.end(), &outputBulk.mesh_meta_data().locally_owned_part());
    std::remove(newParts.begin(), newParts.end(), &outputBulk.mesh_meta_data().aura_part());
    return newParts;
}

void make_nodes_shared(const stk::mesh::BulkData& inputBulk,
                       stk::mesh::Entity oldEntity,
                       stk::mesh::BulkData& outputBulk,
                       stk::mesh::Entity newEntity)
{
    std::vector<int> commShared;
    inputBulk.comm_shared_procs(inputBulk.entity_key(oldEntity), commShared);
    for(int sharedProc : commShared)
        outputBulk.add_node_sharing(newEntity, sharedProc);
}

void copy_bucket_entities(const stk::mesh::BulkData &inputBulk, stk::mesh::Selector inputSelector, const stk::mesh::Bucket *bucket, stk::mesh::BulkData &outputBulk)
{
    stk::mesh::EntityRank rank = bucket->entity_rank();
    stk::mesh::PartVector newParts = get_new_parts(bucket, outputBulk);

    for(stk::mesh::Entity oldEntity : *bucket)
    {
        stk::mesh::Entity newEntity = outputBulk.declare_entity(rank, inputBulk.identifier(oldEntity), newParts);

        if(!is_comm_self(outputBulk) && rank == stk::topology::NODE_RANK && bucket->shared())
            make_nodes_shared(inputBulk, oldEntity, outputBulk, newEntity);
        copy_field_data(inputBulk, rank, oldEntity, newEntity, outputBulk);
    }
}

void copy_side_entities(const stk::mesh::BulkData &inputBulk, stk::mesh::Selector inputSelector, stk::mesh::BulkData &outputBulk)
{
    for(const stk::mesh::Bucket* bucket : inputBulk.get_buckets(inputBulk.mesh_meta_data().side_rank(), inputSelector))
    {
        if(should_copy_bucket(bucket, outputBulk))
        {
            stk::mesh::EntityRank rank = bucket->entity_rank();
            stk::mesh::PartVector newParts = get_new_parts(bucket, outputBulk);

            for(stk::mesh::Entity oldEntity : *bucket)
            {
                stk::mesh::Entity newEntity;
                unsigned numElems = inputBulk.num_elements(oldEntity);
                if(numElems > 0)
                {
                    const stk::mesh::Entity* elems = inputBulk.begin_elements(oldEntity);
                    const stk::mesh::ConnectivityOrdinal* ords = inputBulk.begin_element_ordinals(oldEntity);
                    for(unsigned i = 0; i < numElems; i++)
                    {
                        if(inputBulk.bucket(elems[i]).owned() && inputSelector(inputBulk.bucket(elems[i])))
                        {
                            stk::mesh::Entity newElem = outputBulk.get_entity(inputBulk.entity_key(elems[i]));
                            newEntity = stk::mesh::clone_element_side(outputBulk, inputBulk.identifier(oldEntity), newElem, ords[i], newParts);
                            break;
                        }
                    }
                }

                if(outputBulk.is_valid(newEntity))
                    copy_field_data(inputBulk, rank, oldEntity, newEntity, outputBulk);
            }
        }
    }
}

void copy_sidesets(const stk::mesh::BulkData & inputBulk, stk::mesh::Selector inputSelector, stk::mesh::BulkData & outputBulk)
{
  stk::mesh::MetaData & outputMeta = outputBulk.mesh_meta_data();

  std::vector<const stk::mesh::SideSet *> inputSideSets = inputBulk.get_sidesets();
  for (const stk::mesh::SideSet * inputSideSet : inputSideSets) {
    const std::string & partName = inputSideSet->get_name();
    const stk::mesh::Part * outputPart = outputMeta.get_part(partName);
    ThrowRequire(outputPart != nullptr);

    stk::mesh::SideSet & outputSideSet = outputBulk.create_sideset(*outputPart, inputSideSet->is_from_input());

    for (const stk::mesh::SideSetEntry & inputEntry : *inputSideSet) {
        const stk::mesh::Entity inputElement = inputEntry.element;
        const stk::mesh::ConnectivityOrdinal inputOrdinal = inputEntry.side;

        if (inputSelector(inputBulk.bucket(inputElement))) {
            stk::mesh::Entity outputEntity = outputBulk.get_entity(inputBulk.entity_key(inputElement));
            ThrowRequire(outputBulk.is_valid(outputEntity));
            outputSideSet.add(outputEntity, inputOrdinal);
        }
    }
  }
}

void create_entities_of_rank(const stk::mesh::BulkData& inputBulk, const stk::mesh::Selector& inputSelector, stk::mesh::EntityRank rank, stk::mesh::BulkData& outputBulk)
{
    for(const stk::mesh::Bucket* bucket : inputBulk.get_buckets(rank, inputSelector))
        if(should_copy_bucket(bucket, outputBulk))
            copy_bucket_entities(inputBulk, inputSelector, bucket, outputBulk);
}

void create_entities_for_remaining_ranks(const stk::mesh::BulkData& inputBulk, const stk::mesh::Selector& inputSelector, stk::mesh::BulkData& outputBulk)
{
    stk::mesh::EntityRank beginRank = static_cast<stk::mesh::EntityRank>(stk::topology::ELEM_RANK+1);
    stk::mesh::EntityRank endRank = static_cast<stk::mesh::EntityRank>(inputBulk.mesh_meta_data().entity_rank_count());
    for(stk::mesh::EntityRank rank = beginRank; rank < endRank; rank++)
        create_entities_of_rank(inputBulk, inputSelector, rank, outputBulk);
}

void copy_relations_for_remaining_ranks(const stk::mesh::BulkData& inputBulk, const stk::mesh::Selector& inputSelector, stk::mesh::BulkData& outputBulk)
{
    stk::mesh::EntityRank beginRank = static_cast<stk::mesh::EntityRank>(stk::topology::ELEM_RANK+1);
    stk::mesh::EntityRank endRank = static_cast<stk::mesh::EntityRank>(inputBulk.mesh_meta_data().entity_rank_count());
    for(stk::mesh::EntityRank rank = beginRank; rank < endRank; rank++)
        for(stk::mesh::EntityRank toRank = stk::topology::NODE_RANK; toRank < rank; toRank++)
            copy_relations(inputBulk, inputSelector, rank, toRank, outputBulk);
}

void copy_selected(const stk::mesh::BulkData& inputBulk, const stk::mesh::Selector& inputSelector, stk::mesh::BulkData& outputBulk)
{
    outputBulk.modification_begin();
    create_entities_of_rank(inputBulk, inputSelector, stk::topology::NODE_RANK, outputBulk);
    create_entities_of_rank(inputBulk, inputSelector, stk::topology::ELEM_RANK, outputBulk);
    copy_relations(inputBulk, inputSelector, stk::topology::ELEM_RANK, stk::topology::NODE_RANK, outputBulk);
    outputBulk.modification_end();

    if(inputBulk.has_face_adjacent_element_graph()) {
        outputBulk.initialize_face_adjacent_element_graph();
    }

    outputBulk.modification_begin();
    create_entities_of_rank(inputBulk, inputSelector, stk::topology::EDGE_RANK, outputBulk);

    if(inputBulk.mesh_meta_data().side_rank() != stk::topology::EDGE_RANK) {
      copy_side_entities(inputBulk, inputSelector, outputBulk);
    }

    copy_sidesets(inputBulk, inputSelector, outputBulk);
    create_entities_for_remaining_ranks(inputBulk, inputSelector, outputBulk);

    copy_relations(inputBulk, inputSelector, stk::topology::ELEM_RANK, stk::topology::FACE_RANK, outputBulk);
    copy_relations(inputBulk, inputSelector, stk::topology::ELEM_RANK, stk::topology::EDGE_RANK, outputBulk);
    copy_relations(inputBulk, inputSelector, stk::topology::FACE_RANK, stk::topology::EDGE_RANK, outputBulk);
    copy_relations(inputBulk, inputSelector, stk::topology::FACE_RANK, stk::topology::NODE_RANK, outputBulk);
    copy_relations(inputBulk, inputSelector, stk::topology::EDGE_RANK, stk::topology::NODE_RANK, outputBulk);
    copy_relations_for_remaining_ranks(inputBulk, inputSelector, outputBulk);
    outputBulk.modification_end();
}

void copy_bulk(const stk::mesh::BulkData &inputBulk, stk::mesh::Selector inputSelector, stk::mesh::BulkData &outputBulk)
{
    outputBulk.modification_begin();
    copy_selected(inputBulk, inputSelector, outputBulk);
    outputBulk.modification_end();

    outputBulk.set_large_ids_flag(inputBulk.supports_large_ids());
}

void copy_meta_with_io_attributes(const stk::mesh::MetaData &inputMeta, stk::mesh::MetaData &outputMeta)
{
    copy_meta(inputMeta, outputMeta);
    copy_io_attributes(inputMeta, outputMeta);
}

bool is_connected_to_element(const stk::mesh::BulkData &bulk, stk::mesh::Entity entity)
{
    unsigned numElements = bulk.num_elements(entity);
    const stk::mesh::Entity *elems = bulk.begin_elements(entity);
    for(unsigned i = 0; i < numElements; i++)
        if(bulk.bucket(elems[i]).owned())
            return true;
    return false;
}

bool is_preexisting_orphan(const stk::mesh::BulkData &inBulk, const stk::mesh::BulkData &outputBulk, stk::mesh::EntityRank rank, stk::mesh::Entity entity)
{
    stk::mesh::Entity originalEntity = inBulk.get_entity(rank, outputBulk.identifier(entity));
    return !is_connected_to_element(inBulk, originalEntity);
}

void destroy_all_orphans(const stk::mesh::BulkData &inBulk, stk::mesh::BulkData &outputBulk)
{
    stk::mesh::EntityVector entitiesToDelete;
    for(stk::mesh::EntityRank rank : {stk::topology::FACE_RANK, stk::topology::EDGE_RANK, stk::topology::NODE_RANK})
    {
        for(const stk::mesh::Bucket* bucket : outputBulk.buckets(rank))
        {
            if(bucket->owned())
            {
                for(stk::mesh::Entity entity : *bucket)
                {
                    bool shouldDeleteEntity = !is_connected_to_element(outputBulk, entity);
                    if(shouldDeleteEntity)
                    {
                        if(!is_preexisting_orphan(inBulk, outputBulk, rank, entity))
                            entitiesToDelete.push_back(entity);
                    }
                }
            }
        }
    }

    outputBulk.modification_begin();
    for(stk::mesh::Entity entity : entitiesToDelete)
    {
        destroy_upward_connected_aura_entities(outputBulk, entity, outputBulk.entity_rank(entity));
        bool didDestroy = outputBulk.destroy_entity(entity);
        ThrowRequireMsg(didDestroy, "entity key: " << outputBulk.entity_key(entity));
    }
    outputBulk.modification_end();
}


void copy_mesh(const stk::mesh::BulkData &inputBulk, stk::mesh::Selector inputSelector, stk::mesh::BulkData &outputBulk)
{
    ThrowRequireMsg(&inputBulk != &outputBulk, "Can't copy to same mesh.");
    ThrowRequireMsg(inputBulk.in_modifiable_state() == false, "Can't copy mesh during modification.");
    copy_meta_with_io_attributes(inputBulk.mesh_meta_data(), outputBulk.mesh_meta_data());
    copy_bulk(inputBulk, inputSelector, outputBulk);
    destroy_all_orphans(inputBulk, outputBulk);
}

}
}
