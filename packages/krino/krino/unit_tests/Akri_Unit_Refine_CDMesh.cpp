#include <Akri_MeshSpecs.hpp>
#include <Akri_Unit_RefinementFixture.hpp>
#include <Akri_CDMesh.hpp>
#include <Akri_Interface_Name_Generator.hpp>
#include <Akri_IO_Helpers.hpp>
#include <Akri_Phase_Support.hpp>
#include <Akri_PhaseTag.hpp>
#include <Akri_RefinementSupport.hpp>

namespace krino {

class RegularTriRefinementWithCDMesh : public RefinementFixture<RegularTri>
{
public:
  RegularTriRefinementWithCDMesh()
  {
    StkMeshTriFixture::build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1});
    setup_phases();
    cdmesh = std::make_unique<CDMesh>(mMesh);
    mMesh.mesh_meta_data().enable_late_fields();
    RefinementSupport::get(mMesh.mesh_meta_data()).activate_nonconformal_adaptivity(2);
    RefinementSupport::get(mMesh.mesh_meta_data()).set_non_interface_conforming_refinement(KrinoRefinement::create(mMesh.mesh_meta_data()));
  }

  void setup_phases()
  {
    stk::mesh::Part & block_1 = get_aux_meta().get_part("block_1");
    std::vector<stk::mesh::Part *> blocks{&block_1};

    phaseA.add(surfaceIdentifier, 1);
    phaseB.add(surfaceIdentifier, -1);
    PhaseVec named_phases{{"A", phaseA}, {"B", phaseB}};

    Phase_Support & phase_support = Phase_Support::get(mMesh.mesh_meta_data());
    Block_Surface_Connectivity block_surface_info;
    phase_support.set_input_block_surface_connectivity(block_surface_info);
    phase_support.register_blocks_for_level_set(surfaceIdentifier, blocks);
    std::vector<std::tuple<stk::mesh::PartVector,
      std::shared_ptr<Interface_Name_Generator>, PhaseVec>> ls_sets;
    auto interface_name_gen = std::shared_ptr<Interface_Name_Generator>(new LS_Name_Generator());
    ls_sets.push_back(std::make_tuple(blocks, interface_name_gen, named_phases));
    phase_support.decompose_blocks(ls_sets);
  }

  size_t get_num_elements_in_block(const stk::mesh::Part & block) const
  {
    const stk::mesh::Selector locallyOwnedBlockSelector = block & mMesh.mesh_meta_data().locally_owned_part();
    size_t numInBlock = stk::mesh::count_entities(mMesh, stk::topology::ELEM_RANK, locallyOwnedBlockSelector);
    const size_t localNumInBlock= numInBlock;
    stk::all_reduce_sum(mMesh.parallel(), &localNumInBlock, &numInBlock, 1);
    return numInBlock;
  }

  void move_elements_to_phase(const std::vector<stk::mesh::EntityId> & elemIds, const PhaseTag & phase)
  {
    stk::mesh::PartVector add_parts;
    stk::mesh::PartVector remove_parts;

    mMesh.modification_begin();
    for (auto elemId : elemIds)
    {
      const stk::mesh::Entity elem = mMesh.get_entity(stk::topology::ELEMENT_RANK, elemId);
      if (mMesh.is_valid(elem) && mMesh.bucket(elem).owned())
      {
        cdmesh->determine_conformal_parts(elem, phase, add_parts, remove_parts);
        remove_parts.push_back(&cdmesh->get_parent_part());
        remove_parts.push_back(&cdmesh->get_child_part());

        mMesh.change_entity_parts(elem, add_parts, remove_parts);
      }
    }
    mMesh.modification_end();
  }

  void update_adaptivity_parent_entities()
  {
    mMesh.modification_begin();
    cdmesh->update_adaptivity_parent_entities();
    mMesh.modification_end();
  }

protected:
  std::unique_ptr<CDMesh> cdmesh;
  const Surface_Identifier surfaceIdentifier{0};
  PhaseTag phaseA;
  PhaseTag phaseB;
};


TEST_F(RegularTriRefinementWithCDMesh, meshWithRebalancedElements_phaseChangedOnLeafChildren_adaptivityParentHaveCorrectPhase)
{
  move_elements_to_phase({1001}, phaseA);

  perform_iterations_of_uniform_refinement_with_uniform_marker(2);

  if (stk::parallel_machine_size(mComm) > 1)
    move_owned_elements_with_given_ids_and_owned_attached_entities_to_processor({1003, 1010, 1011, 1012, 1013}, 1);

  move_elements_to_phase({1010, 1011, 1012, 1013}, phaseB);
  update_adaptivity_parent_entities();

  stk::mesh::Part & block_1_A = get_aux_meta().get_part("block_1_A");
  stk::mesh::Part & block_1_B = get_aux_meta().get_part("block_1_B");
  stk::mesh::Part & block_1_nonconformal = get_aux_meta().get_part("block_1_nonconformal");

  EXPECT_EQ(15u, get_num_elements_in_block(block_1_A));
  EXPECT_EQ(5u, get_num_elements_in_block(block_1_B));
  EXPECT_EQ(1u, get_num_elements_in_block(block_1_nonconformal));
}

}
