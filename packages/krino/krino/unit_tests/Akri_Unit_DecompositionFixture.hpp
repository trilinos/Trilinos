/*
 * Akri_Unit_DecompositionFixture.hpp
 *
 *  Created on: Apr 18, 2023
 *      Author: drnoble
 */

#ifndef KRINO_KRINO_UNIT_TESTS_AKRI_UNIT_DECOMPOSITIONFIXTURE_HPP_
#define KRINO_KRINO_UNIT_TESTS_AKRI_UNIT_DECOMPOSITIONFIXTURE_HPP_
#include <Akri_StkMeshFixture.hpp>
#include <stk_io/IossBridge.hpp>
#include <Akri_CDFEM_Support.hpp>
#include <Akri_CDMesh.hpp>
#include <Akri_CreateInterfaceGeometry.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_Unit_LogRedirecter.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_LevelSetPolicy.hpp>

namespace krino {

template <typename MESHSPEC, typename LS_FIELD_POLICY, unsigned NUM_LS>
class DecompositionFixture : public StkMeshFixture<MESHSPEC::TOPOLOGY>
{
public:
using StkMeshFixture<MESHSPEC::TOPOLOGY>::mMesh;
using StkMeshFixture<MESHSPEC::TOPOLOGY>::mBuilder;
using StkMeshFixture<MESHSPEC::TOPOLOGY>::mComm;

DecompositionFixture()
{
  setup_cdfem_support();

  cdmesh = std::make_unique<CDMesh>(mMesh);
}

CDFEM_Support & cdfem_support() { return CDFEM_Support::get(mMesh.mesh_meta_data()); }

FieldRef get_coordinates_field() { return this->get_aux_meta().get_current_coordinates(); }

void setup_ls_fields_with_options(const bool isDeath, const bool doRegisterField)
{
  mMesh.mesh_meta_data().enable_late_fields();

  Block_Surface_Connectivity blockSurfaceConnectivity(mMesh.mesh_meta_data());

  if (isDeath)
  {
    deathSpec.reset(new CDFEM_Inequality_Spec("death_spec"));
    myLSFields = LS_FIELD_POLICY::setup_levelsets_on_blocks(mMesh.mesh_meta_data(), NUM_LS, mBuilder.get_block_parts(), blockSurfaceConnectivity, doRegisterField, deathSpec.get());
  }
  else
  {
    myLSFields = LS_FIELD_POLICY::setup_levelsets_on_blocks(mMesh.mesh_meta_data(), NUM_LS, mBuilder.get_block_parts(), blockSurfaceConnectivity, doRegisterField);
  }
}

void setup_ls_fields_for_death()
{
  setup_ls_fields_with_options(true, true);
}

void setup_ls_fields_without_registering_fields()
{
  setup_ls_fields_with_options(false, false);
}

void setup_ls_fields_for_death_without_registering_fields()
{
  setup_ls_fields_with_options(true, false);
}

void setup_ls_fields()
{
  setup_ls_fields_with_options(false, true);
}

std::vector<LS_Field> & levelset_fields()
{
  return myLSFields;
}

FieldRef get_ls_field(const unsigned lsIndex=0)
{
  STK_ThrowRequire(lsIndex < myLSFields.size());
  return myLSFields[lsIndex].isovar;
}

stk::mesh::Entity get_node(const unsigned nodeIndex)
{
  return this->get_assigned_node_for_index(nodeIndex);
}

double & node_ls_field(FieldRef lsField, const stk::mesh::Entity node)
{
  return *field_data<double>(lsField, node);
}

double & node_ls_field(const stk::mesh::Entity node)
{
  STK_ThrowRequire(1 == myLSFields.size());
  return *field_data<double>(get_ls_field(0), node);
}

void set_level_set(const std::vector<unsigned> & nodeIndices,
    const std::vector<double> & nodeLS)
{
  STK_ThrowRequire(1 == myLSFields.size());
  STK_ThrowRequire(nodeIndices.size() == nodeLS.size());
  for (size_t i = 0; i < nodeIndices.size(); ++i)
    node_ls_field(get_node(nodeIndices[i])) = nodeLS[i];
}

void attempt_decompose_mesh()
{
  try
  {
    decompose_mesh();
  }
  catch (const std::exception & exception)
  {
    std::cout << "Decomposing mesh failed with exception:\n";
    std::cout << exception.what() << "\n";
    ASSERT_TRUE(false);
  }
}

void check_mesh_consistency()
{
  EXPECT_TRUE(check_induced_parts(mMesh));
  EXPECT_TRUE(check_face_and_edge_ownership(mMesh));
  EXPECT_TRUE(check_face_and_edge_relations(mMesh));
  EXPECT_TRUE(check_shared_entity_nodes(mMesh));
  EXPECT_TRUE(cdmesh->check_element_side_parts());
}

void check_nonfatal_error(const std::string & baseName, const int iteration)
{
  if (this->HasNonfatalFailure())
  {
    std::cout << "Failure on iteration " << iteration << std::endl;
    cdmesh->debug_output();
    std::cout << log.get_log() << std::endl;
    std::ostringstream fname;
    fname << baseName << iteration << ".e";
    this->write_mesh(fname.str());
    ASSERT_TRUE(false);
  }
}

void decompose_mesh()
{
  NodeToCapturedDomainsMap nodesToSnappedDomains;
  std::unique_ptr<InterfaceGeometry> interfaceGeometry = create_levelset_geometry(mMesh.mesh_meta_data().spatial_dimension(), cdmesh->get_active_part(), cdfem_support(), Phase_Support::get(mMesh.mesh_meta_data()), levelset_fields());
  if (cdfem_support().get_cdfem_edge_degeneracy_handling() == SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE)
  {
    const double minIntPtWeightForEstimatingCutQuality = cdfem_support().get_snapper().get_edge_tolerance();
    nodesToSnappedDomains = snap_as_much_as_possible_while_maintaining_quality(cdmesh->stk_bulk(),
        cdmesh->get_active_part(),
        cdfem_support().get_snap_fields(),
        *interfaceGeometry,
        cdfem_support().get_global_ids_are_parallel_consistent(),
        cdfem_support().get_snapping_sharp_feature_angle_in_degrees(),
        minIntPtWeightForEstimatingCutQuality,
        cdfem_support().get_max_edge_snap());
  }
  interfaceGeometry->prepare_to_decompose_elements(mMesh, nodesToSnappedDomains);

  cdmesh->generate_nonconformal_elements();
  if (cdfem_support().get_cdfem_edge_degeneracy_handling() == SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE)
    cdmesh->snap_nearby_intersections_to_nodes(*interfaceGeometry, nodesToSnappedDomains);
  cdmesh->set_phase_of_uncut_elements(*interfaceGeometry);
  cdmesh->triangulate(*interfaceGeometry);
  cdmesh->decompose(*interfaceGeometry);
  cdmesh->stash_field_data(-1);
  cdmesh->modify_mesh();
  cdmesh->prolongation();

  if (krinolog.shouldPrint(LOG_DEBUG))
  {
    cdmesh->debug_output();
  }
}

void expect_num_entities(const unsigned numGold, const stk::mesh::Selector & select, const stk::mesh::EntityRank entityRank)
{
  EXPECT_EQ(numGold, stk::mesh::count_selected_entities(select, mMesh.buckets(entityRank))) << "Mismatch for rank " << entityRank << " with selector " << select;
}

void expect_num_elements(const unsigned numGold, const stk::mesh::Selector & select)
{
  expect_num_entities(numGold, select, stk::topology::ELEMENT_RANK);
}

void expect_num_elements(const unsigned numGold)
{
  expect_num_elements(numGold, mMesh.mesh_meta_data().universal_part());
}

void expect_num_sides(const unsigned numGold, const stk::mesh::Selector & select)
{
  expect_num_entities(numGold, select, mMesh.mesh_meta_data().side_rank());
}

void expect_num_sides(const unsigned numGold)
{
  expect_num_sides(numGold, mMesh.mesh_meta_data().universal_part());
}

void expect_num_nodes(const unsigned numGold, const stk::mesh::Selector & select)
{
  expect_num_entities(numGold, select, stk::topology::NODE_RANK);
}

void expect_num_nodes(const unsigned numGold)
{
  expect_num_nodes(numGold, mMesh.mesh_meta_data().universal_part());
}

protected:
void setup_cdfem_support()
{
  FieldRef coordsField = mMesh.mesh_meta_data().coordinate_field();
  cdfem_support().set_coords_field(coordsField);
  cdfem_support().add_edge_interpolation_field(coordsField);
  cdfem_support().register_parent_node_ids_field();

  cdfem_support().set_prolongation_model(INTERPOLATION);
}

protected:
MESHSPEC meshSpec;
std::unique_ptr<CDMesh> cdmesh;
std::unique_ptr<CDFEM_Inequality_Spec> deathSpec;
std::vector<LS_Field> myLSFields;
LogRedirecter log;

};

}

#endif /* KRINO_KRINO_UNIT_TESTS_AKRI_UNIT_DECOMPOSITIONFIXTURE_HPP_ */
