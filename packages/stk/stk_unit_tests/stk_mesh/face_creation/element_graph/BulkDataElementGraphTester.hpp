#ifndef bulkdatalementgraphtester_hpp
#define bulkdatalementgraphtester_hpp

#include <vector>
#include <algorithm>
#include <stdlib.h>

#include <stk_topology/topology.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/baseImpl/MeshImplUtils.hpp>
#include <stk_mesh/base/FieldTraits.hpp>
#include <stk_mesh/base/Field.hpp>

#include <stk_mesh/base/SkinMesh.hpp>
#include <stk_mesh/base/CreateFaces.hpp>
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp>
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraphImpl.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelVectorConcat.hpp>
#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/environment/memory_util.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <stk_util/environment/ReportHandler.hpp>

#include <stk_io/IossBridge.hpp>
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker

#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_unit_test_utils/getOption.h>

#include "stk_unit_test_utils/ElemGraphTestUtils.hpp"
#include "stk_unit_test_utils/unittestMeshUtils.hpp"
#include <stk_unit_tests/stk_mesh/SetupKeyholeMesh.hpp>

#include <stk_mesh/fixtures/QuadFixture.hpp>  // for QuadFixture
#include <stk_mesh/fixtures/heterogeneous_mesh.hpp>
#include <stk_mesh/fixtures/degenerate_mesh.hpp>

#include <stk_mesh/baseImpl/DeletedElementInfo.hpp>

class BulkDataElementGraphTester : public stk::mesh::BulkData
{
public:
    BulkDataElementGraphTester(stk::mesh::MetaData &mesh_meta_data, MPI_Comm comm)
    : stk::mesh::BulkData(mesh_meta_data, comm)
    {
    }

    BulkDataElementGraphTester(stk::mesh::MetaData &meta, MPI_Comm comm, enum AutomaticAuraOption autoAuraOption)
    : stk::mesh::BulkData(meta, comm, autoAuraOption)
    {
    }

    ~BulkDataElementGraphTester()
    {
    }

    bool my_internal_modification_end_for_skin_mesh(stk::mesh::EntityRank entity_rank, stk::mesh::impl::MeshModification::modification_optimization opt, stk::mesh::Selector selectedToSkin,
            const stk::mesh::Selector * only_consider_second_element_from_this_selector = 0)
    {
        return this->internal_modification_end_for_skin_mesh(entity_rank, opt, selectedToSkin, only_consider_second_element_from_this_selector);
    }

    bool my_modification_end_for_entity_creation(const std::vector<stk::mesh::sharing_info>& shared_modified, stk::mesh::impl::MeshModification::modification_optimization opt = stk::mesh::impl::MeshModification::MOD_END_SORT)
    {
        if ( this->in_synchronized_state() ) { return false ; }

        ThrowAssertMsg(stk::mesh::impl::check_for_connected_nodes(*this)==0, "BulkData::modification_end ERROR, all entities with rank higher than node are required to have connected nodes.");

        if (parallel_size() > 1)
        {
            stk::mesh::OrdinalVector shared_part, owned_part, empty;
            shared_part.push_back(m_mesh_meta_data.globally_shared_part().mesh_meta_data_ordinal());
            owned_part.push_back(m_mesh_meta_data.locally_owned_part().mesh_meta_data_ordinal());

            stk::mesh::EntityVector modified_entities(shared_modified.size());
            for(size_t i = 0; i < shared_modified.size(); ++i)
            {
                stk::mesh::Entity entity = shared_modified[i].m_entity;
                int sharing_proc = shared_modified[i].m_sharing_proc;
                entity_comm_map_insert(entity, stk::mesh::EntityCommInfo(stk::mesh::BulkData::SHARED, sharing_proc));
                int owning_proc = shared_modified[i].m_owner;
                const bool am_not_owner = this->internal_set_parallel_owner_rank_but_not_comm_lists(entity, owning_proc);
                if (am_not_owner)
                {
                    stk::mesh::EntityKey key = this->entity_key(entity);
                    internal_change_owner_in_comm_data(key, owning_proc);
                    internal_change_entity_parts(entity, shared_part /*add*/, owned_part /*remove*/);
                }
                else
                {
                    internal_change_entity_parts(entity, shared_part /*add*/, empty /*remove*/);
                }
                modified_entities[i] = entity;
            }

            std::sort(modified_entities.begin(), modified_entities.end(), stk::mesh::EntityLess(*this));
            stk::mesh::EntityVector::iterator iter = std::unique(modified_entities.begin(), modified_entities.end());
            modified_entities.resize(iter-modified_entities.begin());

            add_comm_list_entries_for_entities( modified_entities );

            internal_resolve_shared_membership();

            if ( is_automatic_aura_on() )
            {
                bool connectFacesToPreexistingGhosts = true;
                resolve_incremental_ghosting_for_entity_creation_or_skin_mesh(mesh_meta_data().side_rank(), mesh_meta_data().universal_part(), connectFacesToPreexistingGhosts);
            }

            check_mesh_consistency();
        }
        else
        {
            std::vector<stk::mesh::Entity> modified_entities ;
            internal_update_sharing_comm_map_and_fill_list_modified_shared_entities_of_rank(mesh_meta_data().side_rank(), modified_entities);
        }

        this->internal_finish_modification_end(opt);

        return true ;
    }

    size_t num_entity_keys() const
    {
        return m_entity_keys.size();
    }

    std::vector<uint64_t> my_internal_get_ids_in_use_this_proc_for_locally_owned(stk::topology::rank_t rank) const
    {
        return internal_get_ids_in_use(rank);
    }
};

#endif
