#include <gtest/gtest.h>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include "stk_util/parallel/ParallelReduce.hpp"
#include "stk_mesh/base/Comm.hpp"
#include "stk_mesh/base/FieldBase.hpp"
#include "stk_io/IossBridge.hpp"
#include "stk_unit_test_utils/getOption.h"
#include "stk_balance/internal/privateDeclarations.hpp"
#include "stk_balance/balance.hpp"
#include "stk_mesh/base/Selector.hpp"
#include "stk_balance/internal/Zoltan2ParallelGraph.hpp"

namespace
{

class SpiderElementMesh : public stk::unit_test_util::MeshFixture
{
protected:
    void setup_spider_mesh(const std::string & meshSpec, stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        setup_empty_mesh(auraOption);

        stk::mesh::Part& spherePart = get_meta().declare_part_with_topology("spheres", stk::topology::PARTICLE);
        stk::mesh::Part& beamPart = get_meta().declare_part_with_topology("beams", stk::topology::BEAM_2);

        stk::io::put_io_part_attribute(spherePart);
        stk::io::put_io_part_attribute(beamPart);

        const int initValue = 0;
        stk::mesh::Field<int> & beamField = get_meta().declare_field<stk::mesh::Field<int>>(stk::topology::NODE_RANK,
                                                                                            m_balanceSettings.getSpiderBeamConnectivityCountFieldName());
        stk::mesh::put_field_on_mesh(beamField, get_meta().universal_part(), &initValue);

        stk::mesh::Field<int> & volumeField = get_meta().declare_field<stk::mesh::Field<int>>(stk::topology::ELEM_RANK,
                                                                                              m_balanceSettings.getSpiderVolumeConnectivityCountFieldName());
        stk::mesh::put_field_on_mesh(volumeField, get_meta().universal_part(), &initValue);

        stk::io::fill_mesh(meshSpec, get_bulk());

        std::vector<size_t> entity_counts;
        stk::mesh::comm_mesh_counts(get_bulk(), entity_counts);
        const size_t numNodes = entity_counts[stk::topology::NODE_RANK];
        const size_t numElems = entity_counts[stk::topology::ELEM_RANK];

        const stk::mesh::EntityId sphereNodeID = numNodes+1;

        get_bulk().modification_begin();

        if (get_bulk().parallel_rank() == 0) {
            stk::mesh::EntityId elemID = numElems+1;
            stk::mesh::EntityIdVector nodes = {sphereNodeID};
            stk::mesh::declare_element(get_bulk(), spherePart, elemID, nodes);
        }
        else {
            get_bulk().declare_node(sphereNodeID);
        }

        stk::mesh::Entity sphereNode = get_bulk().get_entity(stk::topology::NODE_RANK, sphereNodeID);
        stk::mesh::FieldBase * coordsField = get_meta().get_field(stk::topology::NODE_RANK,"coordinates");
        EXPECT_TRUE(coordsField != nullptr);
        double * coords = reinterpret_cast<double*>(stk::mesh::field_data(*coordsField, sphereNode));
        coords[0] = 0.0;
        coords[1] = -1.0;
        coords[2] = 0.0;

        for (int otherProc = 0; otherProc < get_bulk().parallel_size(); ++otherProc) {
            if (otherProc != get_bulk().parallel_rank()) {
                get_bulk().add_node_sharing(sphereNode, otherProc);
            }
        }

        std::vector<stk::mesh::EntityId> newBeamIds;
        get_bulk().generate_new_ids(stk::topology::ELEM_RANK, numNodes, newBeamIds);

        size_t idIndex = 0;
        stk::mesh::EntityVector nodes;
        stk::mesh::get_selected_entities(get_meta().locally_owned_part(), get_bulk().buckets(stk::topology::NODE_RANK), nodes);
        for (stk::mesh::Entity node : nodes) {
            double * nodeCoords = reinterpret_cast<double*>(stk::mesh::field_data(*coordsField, node));
            if (nodeCoords[1] < 0.5 && nodeCoords[1] > -0.5) {
                stk::mesh::EntityIdVector beamNodes = {sphereNodeID, get_bulk().identifier(node)};
                stk::mesh::declare_element(get_bulk(), beamPart, newBeamIds[idIndex++], beamNodes);
            }
        }

        get_bulk().modification_end();

        stk::mesh::EntityProcVec beamsToMove;
        if (get_bulk().parallel_rank() != 0) {
            stk::mesh::EntityVector beams;
            stk::mesh::get_selected_entities(beamPart & get_meta().locally_owned_part(), get_bulk().buckets(stk::topology::ELEM_RANK), beams);
            for (stk::mesh::Entity beam : beams) {
                beamsToMove.push_back(std::make_pair(beam, 0));
            }
        }
        get_bulk().change_entity_owner(beamsToMove);
    }

    stk::balance::GraphCreationSettings m_balanceSettings;
};

TEST_F(SpiderElementMesh, move_spider_legs_to_volume_elem_proc)
{
    if (get_parallel_size() > 4) return;

    m_balanceSettings.setShouldFixSpiders(true);
    std::string meshSpec = stk::unit_test_util::get_option("--mesh-spec", "generated:30x3x30");
    setup_spider_mesh(meshSpec, stk::mesh::BulkData::NO_AUTO_AURA);

    stk::balance::balanceStkMesh(m_balanceSettings, get_bulk());

    stk::mesh::EntityVector beams;
    stk::mesh::get_selected_entities(get_meta().get_topology_root_part(stk::topology::BEAM_2) & get_meta().locally_owned_part(),
                                     get_bulk().buckets(stk::topology::ELEM_RANK), beams);
    const int localNumBeams = beams.size();
    int globalNumBeams = 0;

    stk::all_reduce_sum(get_bulk().parallel(), &localNumBeams, &globalNumBeams, 1);
    const int targetNumBeams = globalNumBeams / get_bulk().parallel_size();
    const double beamNumError = std::abs(localNumBeams - targetNumBeams) / (double)targetNumBeams;
    EXPECT_LE(beamNumError, 0.1);
}

} // namespace
