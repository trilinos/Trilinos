#include <gtest/gtest.h>
#include <unistd.h>                     // for unlink
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_io/DatabasePurpose.hpp>   // for DatabasePurpose::READ_MESH, etc
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, put_field
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_topology/topology.hpp>    // for topology, etc
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/meshCreationHelpers.hpp>
#include <stk_util/parallel/Parallel.hpp>

namespace
{

class MeshWithSide : public stk::unit_test_util::MeshFixture { };

//-BEGIN
TEST_F(MeshWithSide, whenCheckingSideEquivalency_returnsCorrectPermutation)
{
  if (stk::parallel_machine_size(get_comm()) == 1) {
    setup_mesh("generated:1x1x4|sideset:x", stk::mesh::BulkData::NO_AUTO_AURA);
    stk::mesh::Entity elem1 = get_bulk().get_entity(stk::topology::ELEM_RANK, 1);
    ASSERT_EQ(1u, get_bulk().num_faces(elem1));
    const stk::mesh::Entity side = *get_bulk().begin_faces(elem1);
    const stk::mesh::Permutation perm = *get_bulk().begin_face_permutations(elem1);
    const stk::mesh::ConnectivityOrdinal ordinal = *get_bulk().begin_face_ordinals(elem1);
    const stk::mesh::Entity* sideNodes = get_bulk().begin_nodes(side);
    unsigned numNodes = get_bulk().num_nodes(side);

    stk::EquivalentPermutation equivAndPermutation = stk::mesh::side_equivalent(get_bulk(), elem1, ordinal, sideNodes);
    EXPECT_TRUE(equivAndPermutation.is_equivalent);
    EXPECT_EQ(perm, static_cast<stk::mesh::Permutation>(equivAndPermutation.permutation_number));

    EXPECT_TRUE(stk::mesh::is_side_equivalent(get_bulk(), elem1, ordinal, sideNodes));

    stk::mesh::EquivAndPositive result = stk::mesh::is_side_equivalent_and_positive(get_bulk(), elem1, ordinal, sideNodes, numNodes);
    EXPECT_TRUE(result.is_equiv);
    EXPECT_TRUE(result.is_positive);
  }
}
//-END

}
