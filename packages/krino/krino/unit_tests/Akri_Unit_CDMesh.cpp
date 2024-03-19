// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <string>
#include <memory>
#include <random>

#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_io/IossBridge.hpp>

#include <Akri_BoundingBox.hpp>
#include <Akri_BoundingBoxMesh.hpp>
#include <Akri_CDFEM_Support.hpp>
#include <Akri_CDMesh.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_Interface_Name_Generator.hpp>
#include <Akri_LevelSetPolicy.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_MeshClone.hpp>
#include <Akri_Phase_Support.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_CreateInterfaceGeometry.hpp>
#include <Akri_QualityMetric.hpp>
#include <Akri_RebalanceUtils.hpp>
#include <Akri_Snap.hpp>

#include <Akri_Unit_BoundingBoxMesh.hpp>
#include <Akri_Unit_Single_Element_Fixtures.hpp>
#include <Akri_Unit_LogRedirecter.hpp>
#include <Akri_Quality.hpp>
#include <Akri_RefinementInterface.hpp>
#include <Akri_RefinementSupport.hpp>
#include <Akri_Unit_DecompositionFixture.hpp>
#include <Akri_MeshSpecs.hpp>

namespace krino
{

class LevelSet;

namespace
{
template <class DECOMP_FIXTURE>
void build_one_tet4_mesh(DECOMP_FIXTURE & fixture,
    stk::mesh::Part & elem_part,
    const int parallel_rank,
    const int parallel_size)
{
  stk::mesh::BulkData & mesh = fixture.fixture.bulk_data();
  stk::mesh::MetaData & meta = fixture.fixture.meta_data();
  AuxMetaData & aux_meta = AuxMetaData::get(meta);

  ASSERT_TRUE(parallel_size == 1 || parallel_size == 2);
  STK_ThrowRequire(elem_part.topology() == stk::topology::TETRAHEDRON_4);
  ASSERT_EQ(stk::topology::ELEMENT_RANK, elem_part.primary_entity_rank());

  mesh.modification_begin();
  {
    stk::mesh::PartVector elem_parts(2);
    elem_parts[1] = &aux_meta.active_part();

    std::vector<stk::mesh::EntityId> elem_nodes = {1, 2, 3, 4};

    stk::mesh::Entity element;
    if (parallel_rank == 0)
    {
      elem_parts[0] = &elem_part;
      element = fixture.create_element(elem_parts, 1, elem_nodes);
    }
  }
  mesh.modification_end();

  if (parallel_rank == 0)
  {
    EXPECT_EQ(1u,
        stk::mesh::count_selected_entities(
            meta.universal_part(), mesh.buckets(stk::topology::ELEMENT_RANK)));

    const auto node1 = mesh.get_entity(stk::topology::NODE_RANK, 1);
    const auto node2 = mesh.get_entity(stk::topology::NODE_RANK, 2);
    const auto node3 = mesh.get_entity(stk::topology::NODE_RANK, 3);
    const auto node4 = mesh.get_entity(stk::topology::NODE_RANK, 4);

    ASSERT_TRUE(mesh.is_valid(node1));
    ASSERT_TRUE(mesh.is_valid(node2));
    ASSERT_TRUE(mesh.is_valid(node3));
    ASSERT_TRUE(mesh.is_valid(node4));

    double * node1_coords = field_data<double>(fixture.coord_field, node1);
    double * node2_coords = field_data<double>(fixture.coord_field, node2);
    double * node3_coords = field_data<double>(fixture.coord_field, node3);
    double * node4_coords = field_data<double>(fixture.coord_field, node4);

    node1_coords[0] = 1.;
    node1_coords[1] = 1.;
    node1_coords[2] = 1.;

    node2_coords[0] = -1.;
    node2_coords[1] = 1.;
    node2_coords[2] = -1.;

    node3_coords[0] = 1.;
    node3_coords[1] = -1.;
    node3_coords[2] = -1.;

    node4_coords[0] = -1.;
    node4_coords[1] = -1.;
    node4_coords[2] = 1.;
  }
}

template <class DECOMP_FIXTURE>
void build_two_tet4_mesh_np2(DECOMP_FIXTURE & fixture,
    stk::mesh::Part & elem1_part,
    stk::mesh::Part & elem2_part,
    stk::mesh::Part & surface_part,
    const int parallel_rank,
    const int parallel_size,
    const bool add_side = true,
    const bool build_all_on_P0 = false)
{
  stk::mesh::BulkData & mesh = fixture.fixture.bulk_data();
  stk::mesh::MetaData & meta = fixture.fixture.meta_data();
  AuxMetaData & aux_meta = AuxMetaData::get(meta);

  ASSERT_TRUE(parallel_size == 1 || parallel_size == 2);
  STK_ThrowRequire(elem1_part.topology() == stk::topology::TETRAHEDRON_4);
  ASSERT_EQ(meta.side_rank(), surface_part.primary_entity_rank());
  ASSERT_EQ(stk::topology::ELEMENT_RANK, elem1_part.primary_entity_rank());
  ASSERT_EQ(stk::topology::ELEMENT_RANK, elem2_part.primary_entity_rank());

  stk::mesh::Entity sideEntity;
  mesh.modification_begin();
  {
    stk::mesh::PartVector elem_parts(2);
    elem_parts[1] = &aux_meta.active_part();

    std::vector<stk::mesh::EntityId> elem1_nodes = {1, 3, 2, 4};
    std::vector<stk::mesh::EntityId> elem2_nodes = {1, 2, 3, 5};

    stk::mesh::Entity element;
    if (parallel_rank == 0)
    {
      elem_parts[0] = &elem1_part;
      element = fixture.create_element(elem_parts, 1, elem1_nodes);
    }
    if ((parallel_rank == 1 && !build_all_on_P0) || (parallel_rank == 0 && build_all_on_P0) ||
        parallel_size == 1)
    {
      elem_parts[0] = &elem2_part;
      element = fixture.create_element(elem_parts, 2, elem2_nodes);
    }

    if (parallel_size > 1 && !build_all_on_P0)
    {
      const int opp_rank = parallel_rank == 0 ? 1 : 0;
      for (auto i = 1; i <= 3; ++i)
      {
        const auto node = mesh.get_entity(stk::topology::NODE_RANK, i);
        mesh.add_node_sharing(node, opp_rank);
      }
    }

    if (add_side && (parallel_rank == 0 || !build_all_on_P0))
    {
      sideEntity = mesh.declare_solo_side(7, {&surface_part});
      for (auto i = 1; i <= 3; ++i)
      {
        auto side_node = mesh.get_entity(stk::topology::NODE_RANK, i);
        mesh.declare_relation(sideEntity, side_node, i - 1);
      }
      attach_entity_to_elements(mesh, sideEntity);
    }
  }
  mesh.modification_end();

  if (parallel_rank == 0 || !build_all_on_P0)
  {
    EXPECT_EQ(2u,
        stk::mesh::count_selected_entities(
            meta.universal_part(), mesh.buckets(stk::topology::ELEMENT_RANK)));
    if (add_side)
    {
      EXPECT_EQ(1u,
          stk::mesh::count_selected_entities(
              meta.universal_part(), mesh.buckets(meta.side_rank())));
      EXPECT_EQ(
          1u, stk::mesh::count_selected_entities(surface_part, mesh.buckets(meta.side_rank())));

      EXPECT_EQ(2u, mesh.num_connectivity(sideEntity, stk::topology::ELEMENT_RANK));
    }

    const auto node1 = mesh.get_entity(stk::topology::NODE_RANK, 1);
    const auto node2 = mesh.get_entity(stk::topology::NODE_RANK, 2);
    const auto node3 = mesh.get_entity(stk::topology::NODE_RANK, 3);
    const auto node4 = mesh.get_entity(stk::topology::NODE_RANK, 4);
    const auto node5 = mesh.get_entity(stk::topology::NODE_RANK, 5);

    ASSERT_TRUE(mesh.is_valid(node1));
    ASSERT_TRUE(mesh.is_valid(node2));
    ASSERT_TRUE(mesh.is_valid(node3));
    ASSERT_TRUE(mesh.is_valid(node4));
    ASSERT_TRUE(mesh.is_valid(node5));

    double * node1_coords = field_data<double>(fixture.coord_field, node1);
    double * node2_coords = field_data<double>(fixture.coord_field, node2);
    double * node3_coords = field_data<double>(fixture.coord_field, node3);
    double * node4_coords = field_data<double>(fixture.coord_field, node4);
    double * node5_coords = field_data<double>(fixture.coord_field, node5);

    node1_coords[0] = 0.;
    node1_coords[1] = 0.;
    node1_coords[2] = 0.;

    node2_coords[0] = 0.;
    node2_coords[1] = 0.;
    node2_coords[2] = 1.;

    node3_coords[0] = 0.;
    node3_coords[1] = 1.;
    node3_coords[2] = 0.;

    node4_coords[0] = -1.;
    node4_coords[1] = 0.;
    node4_coords[2] = 0.;

    node5_coords[0] = 1.;
    node5_coords[1] = 0.;
    node5_coords[2] = 0.;
  }
}
} // namespace

template <class MESH_FIXTURE, class LS_FIELD_POLICY, unsigned NUM_LS>
class CompleteDecompositionFixture : public ::testing::Test
{
public:
  CompleteDecompositionFixture() : fixture(), cdfemSupport(CDFEM_Support::get(fixture.meta_data()))
  {
    AuxMetaData & aux_meta = AuxMetaData::get(fixture.meta_data());
    auto & vec_type =
        fixture.meta_data().spatial_dimension() == 2 ? FieldType::VECTOR_2D : FieldType::VECTOR_3D;
    coord_field = aux_meta.register_field("coordinates",
        vec_type,
        stk::topology::NODE_RANK,
        1u,
        1u,
        fixture.meta_data().universal_part());
    cdfemSupport.set_coords_field(coord_field);
    cdfemSupport.add_edge_interpolation_field(coord_field);
    cdfemSupport.register_parent_node_ids_field();

    cdfemSupport.set_prolongation_model(INTERPOLATION);
  }

  void
  register_ls_on_blocks(const stk::mesh::PartVector & blocks, const bool doRegisterField = true)
  {
    myLSFields = LS_FIELD_POLICY::setup_levelsets_on_blocks(fixture.meta_data(), NUM_LS, blocks, block_surface_info, doRegisterField);
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

  stk::mesh::Part & declare_input_block(const std::string & name, const stk::topology topo)
  {
    auto & block_part = fixture.meta_data().declare_part_with_topology(name, topo);
    stk::io::put_io_part_attribute(block_part);
    return block_part;
  }

  stk::mesh::Part & declare_input_surface(const std::string & name,
      const stk::topology topo,
      const std::set<stk::mesh::PartOrdinal> & touching_blocks)
  {
    auto & surface_part = fixture.meta_data().declare_part_with_topology(name, topo);
    stk::io::put_io_part_attribute(surface_part);

    block_surface_info.add_surface(surface_part.mesh_meta_data_ordinal(), touching_blocks);
    return surface_part;
  }

  stk::mesh::PartVector declare_input_surfaces_touching_block(
      const unsigned numSurfaces, const stk::mesh::Part & touchingBlock)
  {
    const stk::topology topo = touchingBlock.topology().side_topology();
    stk::mesh::PartVector surfaces;
    for (unsigned i = 0; i < numSurfaces; ++i)
    {
      const std::string name = "InputSurface" + std::to_string(i + 1);
      auto & surfacePart = fixture.meta_data().declare_part_with_topology(name, topo);
      stk::io::put_io_part_attribute(surfacePart);

      block_surface_info.add_surface(
          surfacePart.mesh_meta_data_ordinal(), {touchingBlock.mesh_meta_data_ordinal()});
      surfaces.push_back(&surfacePart);
    }

    return surfaces;
  }

  void decompose_mesh()
  {
    NodeToCapturedDomainsMap nodesToSnappedDomains;
    std::unique_ptr<InterfaceGeometry> interfaceGeometry = create_levelset_geometry(fixture.meta_data().spatial_dimension(), krino_mesh->get_active_part(), cdfemSupport, Phase_Support::get(fixture.meta_data()), levelset_fields());
    if (cdfemSupport.get_cdfem_edge_degeneracy_handling() == SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE)
    {
      const double minIntPtWeightForEstimatingCutQuality = cdfemSupport.get_snapper().get_edge_tolerance();
      nodesToSnappedDomains = snap_as_much_as_possible_while_maintaining_quality(krino_mesh->stk_bulk(),
          krino_mesh->get_active_part(),
          cdfemSupport.get_snap_fields(),
          *interfaceGeometry,
          cdfemSupport.get_global_ids_are_parallel_consistent(),
          cdfemSupport.get_snapping_sharp_feature_angle_in_degrees(),
          minIntPtWeightForEstimatingCutQuality,
          cdfemSupport.get_max_edge_snap());
    }
    interfaceGeometry->prepare_to_decompose_elements(krino_mesh->stk_bulk(), nodesToSnappedDomains);

    krino_mesh->generate_nonconformal_elements();
    if (cdfemSupport.get_cdfem_edge_degeneracy_handling() ==
        SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE)
      krino_mesh->snap_nearby_intersections_to_nodes(*interfaceGeometry, nodesToSnappedDomains);
    krino_mesh->set_phase_of_uncut_elements(*interfaceGeometry);
    krino_mesh->triangulate(*interfaceGeometry);
    krino_mesh->decompose(*interfaceGeometry);
    krino_mesh->stash_field_data(-1);
    krino_mesh->modify_mesh();
    krino_mesh->prolongation();

    if (krinolog.shouldPrint(LOG_DEBUG))
    {
      krino_mesh->debug_output();
    }
  }

  void debug_output() { krino_mesh->debug_output(); }

  void commit()
  {
    krino_mesh = std::make_unique<CDMesh>(fixture.bulk_data());
  }

  stk::mesh::Entity create_element(stk::mesh::PartVector & elem_parts,
      stk::mesh::EntityId elem_id,
      std::vector<stk::mesh::EntityId> elem_nodes)
  {
    auto elem = stk::mesh::declare_element(fixture.bulk_data(), elem_parts, elem_id, elem_nodes);
    {
      const stk::mesh::Entity * const nodes = fixture.bulk_data().begin_nodes(elem);
      for (unsigned i = 0; i < elem_nodes.size(); ++i)
      {
        EXPECT_EQ(elem_nodes[i], fixture.bulk_data().identifier(nodes[i]));
        if (!fixture.bulk_data().bucket(nodes[i]).member(cdfemSupport.get_active_part()))
          fixture.bulk_data().change_entity_parts(
              nodes[i], stk::mesh::ConstPartVector{&cdfemSupport.get_active_part()}, {});
      }
    }
    return elem;
  }

  void run_rebalance_with(const std::string & decomp_method)
  {
    /* This is a 2 processor test to confirm that we can rebalance a mesh with CDFEM cut elements,
     * and then successfully cut the mesh again. We create an initial mesh with 2 tets both owned
     * by P0, do a decomposition, rebalance so that 1 parent element should end up on each of P0 and
     * P1, and then do a second decomposition.
     */

    stk::ParallelMachine pm = MPI_COMM_WORLD;
    const int parallel_size = stk::parallel_machine_size(pm);
    const int parallel_rank = stk::parallel_machine_rank(pm);

    stk::mesh::BulkData & mesh = fixture.bulk_data();
    stk::mesh::MetaData & meta = fixture.meta_data();
    AuxMetaData & aux_meta = AuxMetaData::get(meta);
    RefinementInterface * refinement = nullptr;

    if (parallel_size != 2) return;

    cdfemSupport.set_cdfem_edge_tol(0.1);
    cdfemSupport.set_simplex_generation_method(CUT_QUADS_BY_GLOBAL_IDENTIFIER);

    const stk::topology tet4 = stk::topology::TETRAHEDRON_4;
    auto & block1_part = declare_input_block("block_1", tet4);
    auto & block2_part = declare_input_block("block_2", tet4);
    auto & surface_part = declare_input_surface("surface_1",
        tet4.side_topology(),
        {block1_part.mesh_meta_data_ordinal(), block2_part.mesh_meta_data_ordinal()});

    register_ls_on_blocks({&block1_part, &block2_part});

    FieldRef elem_weight_field = aux_meta.register_field("element_weight",
        FieldType::REAL,
        stk::topology::ELEMENT_RANK,
        1u,
        1,
        fixture.meta_data().universal_part());

    commit();

    const bool build_all_on_P0 = true;
    build_two_tet4_mesh_np2(*this,
        block1_part,
        block2_part,
        surface_part,
        parallel_rank,
        parallel_size,
        true,
        build_all_on_P0);

    stk::mesh::create_exposed_block_boundary_sides(
        mesh, meta.universal_part(), {&aux_meta.exposed_boundary_part()});

    if (parallel_rank == 0)
    {
      auto & node_buckets = mesh.buckets(stk::topology::NODE_RANK);
      for (auto && b_ptr : node_buckets)
      {
        for (auto && node : *b_ptr)
        {
          const double * coords = field_data<double>(coord_field, node);
          double * ls_data = field_data<double>(get_ls_field(), node);
          if (ls_data) *ls_data = coords[1] - 0.5;
        }
      }
    }
    stk::mesh::communicate_field_data(mesh, {&get_ls_field().field(), &coord_field.field()});

    ASSERT_NO_THROW(decompose_mesh());

    if (parallel_rank == 1)
    {
      EXPECT_EQ(0u, stk::mesh::get_num_entities(mesh));
    }

    const stk::mesh::Selector parent_selector = krino_mesh->get_parent_part();
    if (parallel_rank == 0)
    {
      auto & elem_buckets = mesh.buckets(stk::topology::ELEMENT_RANK);
      for (auto && b_ptr : elem_buckets)
      {
        const bool is_parent = parent_selector(*b_ptr);
        for (auto && elem : *b_ptr)
        {
          double * weight = field_data<double>(elem_weight_field, elem);
          *weight = is_parent ? 1. : 0.;
        }
      }
    }

    rebalance_utils::rebalance_mesh(mesh,
        refinement,
        krino_mesh.get(),
        elem_weight_field.name(),
        coord_field.name(),
        {fixture.meta_data().universal_part()},
        10,
        decomp_method);

    // Both procs should now own 1 parent element and 4 children
    EXPECT_EQ(1u,
        stk::mesh::count_selected_entities(
            fixture.meta_data().locally_owned_part() & parent_selector,
            mesh.buckets(stk::topology::ELEMENT_RANK)));
    EXPECT_EQ(4u,
        stk::mesh::count_selected_entities(
            fixture.meta_data().locally_owned_part() & krino_mesh->get_child_part(),
            mesh.buckets(stk::topology::ELEMENT_RANK)));

    krino_mesh = std::make_unique<CDMesh>(mesh);

    auto & node_buckets = mesh.buckets(stk::topology::NODE_RANK);
    for (auto && b_ptr : node_buckets)
    {
      for (auto && node : *b_ptr)
      {
        const double * coords = field_data<double>(coord_field, node);
        double * ls_data = field_data<double>(get_ls_field(), node);
        if (ls_data) *ls_data = coords[2] - 0.5;
      }
    }

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

    EXPECT_TRUE(check_induced_parts(mesh));
    EXPECT_TRUE(check_face_and_edge_ownership(mesh));
    EXPECT_TRUE(check_face_and_edge_relations(mesh));
    EXPECT_TRUE(check_shared_entity_nodes(mesh));
    EXPECT_TRUE(krino_mesh->check_element_side_parts());
  }

  MESH_FIXTURE fixture;
  FieldRef coord_field;
  CDFEM_Support & cdfemSupport;
  std::unique_ptr<CDMesh> krino_mesh;
  Block_Surface_Connectivity block_surface_info;
  std::vector<LS_Field> myLSFields;
  LogRedirecter log;
};

typedef DecompositionFixture<RegularTri, LSPerInterfacePolicy, 1> RegularTriSingleLSDecompositionFixture;
TEST_F(RegularTriSingleLSDecompositionFixture, decompose)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_size = stk::parallel_machine_size(pm);
  if (parallel_size > 1) return;

  StkMeshTriFixture::build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1});

  setup_ls_fields();

  set_level_set({0,1,2} ,{-1., 1., -1.});

  attempt_decompose_mesh();

  std::vector<stk::mesh::Entity> entities;

  // Should have added 2 nodes at the cutting locations
  mMesh.get_entities(stk::topology::NODE_RANK, mMesh.mesh_meta_data().universal_part(), entities);
  EXPECT_EQ(5u, entities.size());
  mMesh.get_entities(stk::topology::NODE_RANK, get_aux_meta().active_part(), entities);
  EXPECT_EQ(5u, entities.size());

  // Should be 1 interface edge
  mMesh.get_entities(stk::topology::EDGE_RANK,
      get_aux_meta().active_part() & get_aux_meta().get_part("surface_block_1_P-_P+") & get_aux_meta().get_part("surface_block_1_P+_P-"),
      entities);
  EXPECT_EQ(1u, entities.size());

  // Should be 3 conformal elements plus the parent element
  mMesh.get_entities(stk::topology::ELEMENT_RANK, mMesh.mesh_meta_data().universal_part(), entities);
  EXPECT_EQ(4u, entities.size());
  mMesh.get_entities(stk::topology::ELEMENT_RANK, get_aux_meta().active_part(), entities);
  EXPECT_EQ(3u, entities.size());
  mMesh.get_entities(stk::topology::ELEMENT_RANK, !get_aux_meta().active_part(), entities);
  EXPECT_EQ(1u, entities.size());

  mMesh.get_entities(stk::topology::ELEMENT_RANK, get_aux_meta().get_part("block_1_P-"), entities);
  EXPECT_EQ(2u, entities.size());
  mMesh.get_entities(stk::topology::ELEMENT_RANK, get_aux_meta().get_part("block_1_P+"), entities);
  EXPECT_EQ(1u, entities.size());
}

class DecompositionFixtureWithoutRegisteringFields : public RegularTriSingleLSDecompositionFixture
{
public:
  DecompositionFixtureWithoutRegisteringFields()
  {
    StkMeshTriFixture::build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1});
    setup_ls_fields_without_registering_fields();
  }
};

TEST_F(DecompositionFixtureWithoutRegisteringFields, IsovariableNotDefinedOnAnyBlock)
{
  EXPECT_ANY_THROW(Phase_Support::check_isovariable_field_existence_on_decomposed_blocks(
      mMesh.mesh_meta_data(), levelset_fields(), true));
}

TEST_F(DecompositionFixtureWithoutRegisteringFields, IsovariableNotDefinedOnBlock1)
{
  auto & A_part = get_aux_meta().get_part("block_1_P+");
  auto & B_part = get_aux_meta().get_part("block_1_P-");

  // Catch the case where the field exists on both conformal parts, but not
  // the initial un-decomposed part so we can't do the initial decomposition.
  get_aux_meta().register_field("LS", FieldType::REAL, stk::topology::NODE_RANK, 1u, 1u, A_part);
  get_aux_meta().register_field("LS", FieldType::REAL, stk::topology::NODE_RANK, 1u, 1u, B_part);

  EXPECT_ANY_THROW(Phase_Support::check_isovariable_field_existence_on_decomposed_blocks(
      mMesh.mesh_meta_data(), levelset_fields(), true));
}

TEST_F(DecompositionFixtureWithoutRegisteringFields, IsovariableOnlyOnBlock1SteadyState)
{
  auto & blockPart = get_aux_meta().get_part("block_1");

  // Catch the case where the field on the initial un-decomposed part but not on the conformal parts
  get_aux_meta().register_field("LS", FieldType::REAL, stk::topology::NODE_RANK, 1u, 1u, blockPart);

  EXPECT_ANY_THROW(Phase_Support::check_isovariable_field_existence_on_decomposed_blocks(
      mMesh.mesh_meta_data(), levelset_fields(), true));
}

class DecompositionFixtureForDeathWithoutRegisteringFields : public RegularTriSingleLSDecompositionFixture
{
public:
  DecompositionFixtureForDeathWithoutRegisteringFields()
  {
    StkMeshTriFixture::build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1});
    setup_ls_fields_for_death_without_registering_fields();
  }
};

TEST_F(DecompositionFixtureForDeathWithoutRegisteringFields, DeathIsovariableNotDefinedOnDecomposedBlock)
{
  EXPECT_ANY_THROW(Phase_Support::check_isovariable_field_existence_on_decomposed_blocks(
      mMesh.mesh_meta_data(), levelset_fields(), true));
}

TEST_F(DecompositionFixtureForDeathWithoutRegisteringFields, DeathIsovariableNotDefinedOnDeadBlock)
{
  auto & block_part = get_aux_meta().get_part("block_1");
  auto & alive_part = get_aux_meta().get_part("block_1_P-");

  get_aux_meta().register_field("LS", FieldType::REAL, stk::topology::NODE_RANK, 1u, 1u, block_part);
  get_aux_meta().register_field("LS", FieldType::REAL, stk::topology::NODE_RANK, 1u, 1u, alive_part);

  EXPECT_NO_THROW(Phase_Support::check_isovariable_field_existence_on_decomposed_blocks(
      mMesh.mesh_meta_data(), levelset_fields(), true));
}

class TwoRightTrisOn1Or2ProcsSingleLSDecompositionFixture : public DecompositionFixture<TwoRightTrisSharingDiagonal, LSPerInterfacePolicy, 1>
{
public:
TwoRightTrisOn1Or2ProcsSingleLSDecompositionFixture()
{
  if(stk::parallel_machine_size(mComm) == 1)
    this->build_mesh(meshSpec.nodeLocs, {meshSpec.allElementConn});
  else if(stk::parallel_machine_size(mComm) == 2)
    this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,1}, {0,1});
}

void check_snapped_mesh_with_one_element_in_each_phase()
{
  const unsigned parallelSize = stk::parallel_machine_size(mComm);
  const unsigned parallelRank = stk::parallel_machine_rank(mComm);

  std::vector<stk::mesh::Entity> entities;

  // Should be no new nodes because of snapped interface, 3 nodes owned by P0, 1 by P1
  mMesh.get_entities(stk::topology::NODE_RANK, mMesh.mesh_meta_data().universal_part(), entities);
  EXPECT_EQ(4u, entities.size());
  mMesh.get_entities(stk::topology::NODE_RANK, get_aux_meta().active_part(), entities);
  EXPECT_EQ(4u, entities.size());
  mMesh.get_entities(
      stk::topology::NODE_RANK, get_aux_meta().active_part() & mMesh.mesh_meta_data().locally_owned_part(), entities);

  const std::array<std::array<unsigned,2>,2> goldNumOwnedNodesPerProc = {{ {{4,0}}, {{3,1}} }};
  const unsigned goldNumOwnedNodes = goldNumOwnedNodesPerProc[parallelSize-1][parallelRank];
  EXPECT_EQ(goldNumOwnedNodes, entities.size());

  // Should be 1 interface edge
  mMesh.get_entities(stk::topology::EDGE_RANK,
      get_aux_meta().active_part() & get_aux_meta().get_part("surface_block_1_P+_P-") &
          get_aux_meta().get_part("surface_block_1_P-_P+"),
      entities);
  EXPECT_EQ(1u, entities.size());

  // Should be 2 coincident subelements, no parents
  mMesh.get_entities(stk::topology::ELEMENT_RANK, mMesh.mesh_meta_data().universal_part(), entities);
  EXPECT_EQ(2u, entities.size());
  mMesh.get_entities(stk::topology::ELEMENT_RANK, get_aux_meta().active_part(), entities);
  EXPECT_EQ(2u, entities.size());
  mMesh.get_entities(stk::topology::ELEMENT_RANK, !get_aux_meta().active_part(), entities);
  EXPECT_EQ(0u, entities.size());
  mMesh.get_entities(stk::topology::ELEMENT_RANK, get_aux_meta().get_part("block_1_P+"), entities);
  EXPECT_EQ(1u, entities.size());
  mMesh.get_entities(stk::topology::ELEMENT_RANK, get_aux_meta().get_part("block_1_P-"), entities);
  EXPECT_EQ(1u, entities.size());
}
};

TEST_F(TwoRightTrisOn1Or2ProcsSingleLSDecompositionFixture, decompose)
{
  if(stk::parallel_machine_size(mComm) > 2) return;
  setup_ls_fields();

  {
    set_level_set({0,1,2,3} ,{-1., 0., 1., 0.});

    attempt_decompose_mesh();

    check_snapped_mesh_with_one_element_in_each_phase();

    cdmesh = std::make_unique<CDMesh>(mMesh);
  }

  {
    // Swap the signs

    set_level_set({0,1,2,3} ,{1., 0., -1., 0.});

    attempt_decompose_mesh();

    check_snapped_mesh_with_one_element_in_each_phase();
  }
}

TEST_F(TwoRightTrisOn1Or2ProcsSingleLSDecompositionFixture, death_decompose)
{
  if(stk::parallel_machine_size(mComm) > 2) return;
  setup_ls_fields_for_death();

  set_level_set({0,1,2,3} ,{-1., 0., 1., 0.});

  attempt_decompose_mesh();

  check_snapped_mesh_with_one_element_in_each_phase();
}

TEST_F(TwoRightTrisOn1Or2ProcsSingleLSDecompositionFixture, death_swapped_signs_decompose)
{
  if(stk::parallel_machine_size(mComm) > 2) return;
  setup_ls_fields_for_death();

  set_level_set({0,1,2,3} ,{1., 0., -1., 0.});

  attempt_decompose_mesh();

  check_snapped_mesh_with_one_element_in_each_phase();
}

TEST_F(TwoRightTrisOn1Or2ProcsSingleLSDecompositionFixture, periodic_decompose)
{
  if(stk::parallel_machine_size(mComm) > 2) return;
  setup_ls_fields();

  cdfem_support().set_cdfem_edge_tol(0.1);

  set_level_set({0,1,2,3} ,{-0.99, -0.8, 0.2, 0.01});

  cdmesh->add_periodic_node_pair(get_node(2), get_node(3));

  attempt_decompose_mesh();

  // Periodic constraint should cause interface to snap to both 3 and 4 so no elements get cut.
  std::vector<stk::mesh::Entity> entities;
  mMesh.get_entities(stk::topology::ELEMENT_RANK, get_aux_meta().get_part("block_1_P+"), entities);
  EXPECT_EQ(0u, entities.size());
  mMesh.get_entities(stk::topology::ELEMENT_RANK, get_aux_meta().get_part("block_1_P-"), entities);
  EXPECT_EQ(2u, entities.size());
}

TEST_F(TwoRightTrisOn1Or2ProcsSingleLSDecompositionFixture, Check_Compatibility_When_Snapping)
{
  if(stk::parallel_machine_size(mComm) > 2) return;
  setup_ls_fields();

  set_level_set({0,1,2,3}, {-1., -1., 1.e20, 1.});

  attempt_decompose_mesh();

  std::vector<stk::mesh::Entity> entities;
  mMesh.get_entities(stk::topology::ELEMENT_RANK, get_aux_meta().get_part("block_1_P+"), entities);
  EXPECT_EQ(2u, entities.size());
  mMesh.get_entities(stk::topology::ELEMENT_RANK, get_aux_meta().get_part("block_1_P-"), entities);
  EXPECT_EQ(1u, entities.size());

  // fixture.write_results("Two_Tri3_Check_Compatibility_When_Snapping.e");
}

class FourDisconnectedTrisOn2Or4ProcsDecompositionFixture : public DecompositionFixture<FourDisconnectedTris, LSPerInterfacePolicy, 1>
{
public:
FourDisconnectedTrisOn2Or4ProcsDecompositionFixture()
{
  if(stk::parallel_machine_size(mComm) == 2)
    this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,1,1,1}, {0,1,0,1});
  else if(stk::parallel_machine_size(mComm) == 4)
    this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,1,1,1}, {0,1,2,3});

  setup_ls_fields();
}

void setup_ghosting()
{
  mMesh.modification_begin();
  auto & ghosting = mMesh.create_ghosting("test_ghosting");
  mMesh.modification_end();
  stk::mesh::EntityProcVec add_send;
  const int parallel_rank = stk::parallel_machine_rank(mComm);
  const int mod = parallel_rank % 2;
  const stk::mesh::Entity node1 = get_node(3 * parallel_rank);
  if (mod == 0)
  {
    add_send.push_back(std::make_pair(node1, parallel_rank + 1));
  }
  else
  {
    add_send.push_back(std::make_pair(node1, parallel_rank - 1));
  }
  mMesh.modification_begin();
  mMesh.change_ghosting(ghosting, add_send);
  mMesh.modification_end();

  const unsigned goldNumElems = 4/stk::parallel_machine_size(mComm);
  ASSERT_EQ(3*goldNumElems+1, stk::mesh::count_selected_entities(mMesh.mesh_meta_data().universal_part(), mMesh.buckets(stk::topology::NODE_RANK)));

  auto other_node1 = get_node((mod == 0) ? 3*(parallel_rank+1) : 3*(parallel_rank-1));
  cdmesh->add_periodic_node_pair(node1, other_node1);
}

};

TEST_F(FourDisconnectedTrisOn2Or4ProcsDecompositionFixture, PeriodicParallelNonSharedNode)
{
  if(stk::parallel_machine_size(mComm) != 2 && stk::parallel_machine_size(mComm) != 4) return;

  const unsigned goldNumElems = 4/stk::parallel_machine_size(mComm);
  EXPECT_EQ(goldNumElems, stk::mesh::count_selected_entities(mMesh.mesh_meta_data().universal_part(), mMesh.buckets(stk::topology::ELEMENT_RANK)));

  cdfem_support().set_cdfem_edge_tol(0.1);

  setup_ghosting();

  const unsigned parallelRank = stk::parallel_machine_rank(mComm);
  set_level_set({3*parallelRank, 3*parallelRank+1, 3*parallelRank+2}, {(parallelRank % 2 == 0) ? -0.01 : -1., 1., 1.});

  attempt_decompose_mesh();

  // Periodic constraint should cause interface to snap to node1 so the element(s) are entirely in the P+ phase
  std::vector<stk::mesh::Entity> entities;
  mMesh.get_entities(stk::topology::ELEMENT_RANK, get_aux_meta().get_part("block_1_P+"), entities);
  EXPECT_EQ(goldNumElems, entities.size());
  mMesh.get_entities(stk::topology::ELEMENT_RANK, get_aux_meta().get_part("block_1_P-"), entities);
  EXPECT_EQ(0u, entities.size());
}

void set_level_sets(const stk::mesh::BulkData & mesh,
    const std::vector<LS_Field> & lsFields,
    const std::vector<stk::mesh::EntityId> & nodeIds,
    const std::vector<std::vector<double>> & nodeLS)
{
  STK_ThrowRequire(nodeIds.size() == nodeLS.size());
  const size_t numFields = lsFields.size();
  for (size_t i = 0; i < nodeIds.size(); ++i)
  {
    const auto node = mesh.get_entity(stk::topology::NODE_RANK, nodeIds[i]);
    STK_ThrowRequire(mesh.is_valid(node));
    STK_ThrowRequire(numFields == nodeLS[i].size());
    for (size_t j = 0; j < numFields; ++j)
      *field_data<double>(lsFields[j].isovar, node) = nodeLS[i][j];
  }
}

class CDMeshTestsTwoRightTris3LSPerPhase : public DecompositionFixture<TwoRightTrisSharingDiagonal, LSPerPhasePolicy, 3>
{
public:
CDMeshTestsTwoRightTris3LSPerPhase()
{
  if(stk::parallel_machine_size(mComm) <= 2)
  {
    const SideIdAndTriSideNodes sideset1{1, { {{0,1}}, {{2,3}} }};
    const SideIdAndTriSideNodes sideset2{2, { {{1,2}}, {{3,0}} }};
    const std::vector<SideIdAndTriSideNodes> sidesets{ sideset1, sideset2 };
    mBuilder.create_sideset_parts(sidesets);

    if(stk::parallel_machine_size(mComm) == 1)
      this->build_mesh(meshSpec.nodeLocs, {meshSpec.allElementConn});
    else if(stk::parallel_machine_size(mComm) == 2)
      this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,1}, {0,1});

    mBuilder.add_sides_to_sidesets(sidesets);

    setup_ls_fields();
  }
}
};

TEST_F(CDMeshTestsTwoRightTris3LSPerPhase, Two_Tri3_Check_Compatibility_When_Snapping)
{
  if(stk::parallel_machine_size(mComm) > 2) return;

  cdfem_support().set_cdfem_edge_degeneracy_handling(SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE);

  get_aux_meta().set_assert_32bit_flag();
  get_aux_meta().clear_force_64bit_flag();

  const double eps = 1.e-13;
  set_level_sets(mMesh, levelset_fields(),
      get_ids_of_nodes_with_given_indices({0, 1, 2, 3}),
      {{-1., 1., 1. + eps}, {-1., 1., 1. + eps}, {1.e2, 1., -1.}, {1., -1., -1. + eps}});

  attempt_decompose_mesh();

  check_mesh_consistency();

  std::vector<stk::mesh::Entity> entities;
  mMesh.get_entities(stk::topology::NODE_RANK, mMesh.mesh_meta_data().universal_part(), entities);
  EXPECT_EQ(7u, entities.size());
  mMesh.get_entities(stk::topology::ELEMENT_RANK, get_aux_meta().get_part("block_1_P0"), entities);
  EXPECT_EQ(3u, entities.size());
  mMesh.get_entities(stk::topology::ELEMENT_RANK, get_aux_meta().get_part("block_1_P1"), entities);
  EXPECT_EQ(1u, entities.size());
  mMesh.get_entities(stk::topology::ELEMENT_RANK, get_aux_meta().get_part("block_1_P2"), entities);
  EXPECT_EQ(2u, entities.size());

  //write_mesh("Two_Tri3_Check_Compatibility_When_Snapping_LSPerPhase.e");
}



class TwoRegularTetsSharingNodeAtOriginOn1or2ProcsDecompositionFixture : public DecompositionFixture<TwoRegularTetsSharingNodeAtOrigin, LSPerInterfacePolicy, 1>
{
public:
TwoRegularTetsSharingNodeAtOriginOn1or2ProcsDecompositionFixture()
{
  if(stk::parallel_machine_size(mComm) == 1)
    this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,1}, {0,0});
  else if(stk::parallel_machine_size(mComm) == 2)
    this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,1}, {0,1});

  setup_ls_fields();
}
};

TEST_F(TwoRegularTetsSharingNodeAtOriginOn1or2ProcsDecompositionFixture, onlyCutFacesImpactNodeScoringToImproveQuality)
{
  if (stk::parallel_machine_size(mComm) > 2) return;

  set_level_set({0,1,2,3,4,5,6}, {-0.49, 1.5, 1.5, 1.5, -0.06, 0.18, 0.18});

  attempt_decompose_mesh();

  const ScaledJacobianQualityMetric qualityMetric;
  const double qualityAfterCut = compute_mesh_quality(mMesh, get_aux_meta().active_part(), qualityMetric);
  EXPECT_LT(0.07, qualityAfterCut);

  //write_mesh("test.e");
}

class TwoRightTetsWith2BlocksOn1or2ProcsDecompositionFixture : public DecompositionFixture<TwoRightTets, LSPerInterfacePolicy, 1>
{
public:
TwoRightTetsWith2BlocksOn1or2ProcsDecompositionFixture() {}

void build_mesh_with_optional_sideset(const bool addSideset)
{
  const std::array<unsigned, 3> side1Nodes{{0,2,4}};
  const std::vector<SideIdAndTetSideNodes> sideIdsAndNodes{{1, {side1Nodes}}};

  if (addSideset)
    mBuilder.create_sideset_parts(sideIdsAndNodes);

  if(stk::parallel_machine_size(mComm) == 1)
    this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,2}, {0,0});
  else if(stk::parallel_machine_size(mComm) == 2)
    this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,2}, {0,1});

  if (addSideset)
    mBuilder.add_sides_to_sidesets(sideIdsAndNodes);

  setup_ls_fields();
}
};

TEST_F(TwoRightTetsWith2BlocksOn1or2ProcsDecompositionFixture, Write_Results_No_Side)
{
  if(stk::parallel_machine_size(mComm) > 2) return;

  build_mesh_with_optional_sideset(false);

  write_mesh("Write_Results_No_Side.e");
}

TEST_F(TwoRightTetsWith2BlocksOn1or2ProcsDecompositionFixture, Write_Results_With_Side)
{
  if(stk::parallel_machine_size(mComm) > 2) return;

  build_mesh_with_optional_sideset(true);

  write_mesh("Write_Results_With_Side.e");
}

namespace
{
void randomize_ls_field(const stk::mesh::BulkData & mesh,
    const FieldRef & field,
    std::mt19937 & mt,
    std::uniform_real_distribution<double> & dist)
{
  auto & node_buckets = mesh.buckets(stk::topology::NODE_RANK);
  for (auto && b_ptr : node_buckets)
  {
    for (auto && node : *b_ptr)
    {
      double * ls_data = field_data<double>(field, node);
      if (ls_data) *ls_data = dist(mt);
    }
  }
}
void set_ls_field_on_part(const stk::mesh::BulkData & mesh,
    const stk::mesh::Part & part,
    const FieldRef & ls_field,
    const double ls_value)
{
  auto & node_buckets = mesh.get_buckets(stk::topology::NODE_RANK, part);
  for (auto && b_ptr : node_buckets)
  {
    for (auto && node : *b_ptr)
    {
      double * ls_data = field_data<double>(ls_field, node);
      if (ls_data) *ls_data = ls_value;
    }
  }
}

} // namespace

TEST_F(TwoRightTetsWith2BlocksOn1or2ProcsDecompositionFixture, Random_TwoTet4_InternalSideset)
{
  if(stk::parallel_machine_size(mComm) > 2) return;

  build_mesh_with_optional_sideset(true);

  // Use a large snap tolerance to make snapping more common since it is a frequent source of parallel bugs
  cdfem_support().set_cdfem_edge_tol(0.1);
  cdfem_support().set_simplex_generation_method(CUT_QUADS_BY_GLOBAL_IDENTIFIER);

  std::mt19937 mt(std::mt19937::default_seed + stk::parallel_machine_rank(mComm));
  std::uniform_real_distribution<double> dist(-1., 1.);
#ifdef NDEBUG
  const int num_cases = 5000;
#else
  const int num_cases = 1000;
#endif
  for (int i = 0; i < num_cases; ++i)
  {
    if (i % 1000 == 0) std::cout << "Testing random configuration " << i << std::endl;

    randomize_ls_field(mMesh, get_ls_field(), mt, dist);
    stk::mesh::communicate_field_data(mMesh, {&get_ls_field().field(), &get_coordinates_field().field()});

    attempt_decompose_mesh();

    check_mesh_consistency();
    check_nonfatal_error("Random_TwoTet4_InternalSideset_iter", i);

    cdmesh = std::make_unique<CDMesh>(mMesh);
  }
}

typedef CompleteDecompositionFixture<BoundingBoxMeshTri3, LSPerInterfacePolicy, 1> CDMeshTestsBboxMesh2D;
TEST_F(CDMeshTestsBboxMesh2D, Random_SnapMesh)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_size = stk::parallel_machine_size(pm);
  const int parallel_rank = stk::parallel_machine_rank(pm);

  stk::mesh::MetaData & meta = fixture.meta_data();
  AuxMetaData & aux_meta = AuxMetaData::get(meta);

  if (parallel_size > 2) return;

  cdfemSupport.set_cdfem_edge_degeneracy_handling(
      SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE);
  const double approxMinRelativeSize = 0.25;

  auto & block1_part = aux_meta.get_part("block_1");

  register_ls_on_blocks({&block1_part});

  typename BoundingBoxMesh::BoundingBoxType domain(stk::math::Vector3d::ZERO, stk::math::Vector3d(1., 1., 0.));
  const double mesh_size = 1. / 3.;
  fixture.set_domain(domain, mesh_size);
  fixture.populate_mesh();
  stk::mesh::BulkData & mesh = fixture.bulk_data();
  stk::mesh::create_exposed_block_boundary_sides(
      mesh, meta.universal_part(), {&aux_meta.exposed_boundary_part(), &aux_meta.active_part()});

  double overallMinEdgeLength = std::numeric_limits<double>::max();
  double overallMaxEdgeLength = -std::numeric_limits<double>::max();
  double overallMinVolume = std::numeric_limits<double>::max();
  double overallMaxVolume = -std::numeric_limits<double>::max();
  const double expectedMinLength = mesh_size * approxMinRelativeSize;
  const double expectedMinVol = std::pow(mesh_size * approxMinRelativeSize, 2.) / 2.;

  std::mt19937 mt(std::mt19937::default_seed + parallel_rank);
  std::uniform_real_distribution<double> dist(-1., 1.);
#ifdef NDEBUG
  const int num_cases = 5000;
#else
  const int num_cases = 1000;
#endif
  for (int i = 0; i < num_cases; ++i)
  {
    if (i % 1000 == 0) std::cout << "Testing random configuration " << i << std::endl;

    MeshClone::stash_or_restore_mesh(mesh, 0); // restore original uncut mesh
    commit();                                  // new krino_mesh each time

    randomize_ls_field(mesh, get_ls_field(), mt, dist);
    set_ls_field_on_part(mesh, aux_meta.exposed_boundary_part(), get_ls_field(), 1.);
    stk::mesh::communicate_field_data(mesh, {&get_ls_field().field(), &coord_field.field()});

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

    EXPECT_TRUE(check_induced_parts(mesh));
    EXPECT_TRUE(check_face_and_edge_ownership(mesh));
    EXPECT_TRUE(check_face_and_edge_relations(mesh));
    EXPECT_TRUE(check_shared_entity_nodes(mesh));
    EXPECT_TRUE(krino_mesh->check_element_side_parts());

    double minEdgeLength, maxEdgeLength, minVolume, maxVolume;
    compute_element_quality(mesh, minEdgeLength, maxEdgeLength, minVolume, maxVolume);
    overallMinEdgeLength = std::min(overallMinEdgeLength, minEdgeLength);
    overallMaxEdgeLength = std::max(overallMaxEdgeLength, maxEdgeLength);
    overallMinVolume = std::min(overallMinVolume, minVolume);
    overallMaxVolume = std::max(overallMaxVolume, maxVolume);

    bool failedQuality = false;
    if (minVolume < 0.5 * expectedMinVol)
    {
      failedQuality = true;
      std::cout << "Failed quality requirements: minEdgeLength=" << minEdgeLength
                << ", maxEdgeLength=" << maxEdgeLength << ", minVolume=" << minVolume
                << ", maxVolume=" << maxVolume << std::endl;
    }

    if (HasNonfatalFailure() || failedQuality)
    {
      std::cout << "Failure on iteration " << i << std::endl;
      krino_mesh->debug_output();
      std::cout << log.get_log() << std::endl;
      std::ostringstream fname;
      fname << "Random_SnapTri3_iter_" << i << ".e";
      SimpleStkFixture::write_results(fname.str(), mesh);
      ASSERT_TRUE(false);
    }
  }
  std::cout << "Expected quality: minEdgeLength~=" << expectedMinLength
            << ", minVolume~=" << expectedMinVol << std::endl;
  std::cout << "Quality results: minEdgeLength=" << overallMinEdgeLength
            << ", maxEdgeLength=" << overallMaxEdgeLength << ", minVolume=" << overallMinVolume
            << ", maxVolume=" << overallMaxVolume << std::endl;
}

typedef CompleteDecompositionFixture<BoundingBoxMeshTri3, LSPerPhasePolicy, 3>
    CDMeshTestsBboxMesh2DLSPerPhase;
TEST_F(CDMeshTestsBboxMesh2DLSPerPhase, Random_SnapMesh)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_size = stk::parallel_machine_size(pm);
  const int parallel_rank = stk::parallel_machine_rank(pm);

  stk::mesh::MetaData & meta = fixture.meta_data();
  AuxMetaData & aux_meta = AuxMetaData::get(meta);

  if (parallel_size > 2) return;

  cdfemSupport.set_cdfem_edge_degeneracy_handling(
      SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE);
  const double approxMinRelativeSize = 0.25;

  auto & block1_part = aux_meta.get_part("block_1");

  register_ls_on_blocks({&block1_part});

  typename BoundingBoxMesh::BoundingBoxType domain(stk::math::Vector3d::ZERO, stk::math::Vector3d(1., 1., 0.));
  const double mesh_size = 1. / 3.;
  fixture.set_domain(domain, mesh_size);
  fixture.populate_mesh();
  stk::mesh::BulkData & mesh = fixture.bulk_data();
  stk::mesh::create_exposed_block_boundary_sides(
      mesh, meta.universal_part(), {&aux_meta.exposed_boundary_part(), &aux_meta.active_part()});

  const auto & lsFields = levelset_fields();
  std::vector<const stk::mesh::FieldBase *> sync_fields = {&coord_field.field()};
  for (auto && lsField : lsFields)
    sync_fields.push_back(&lsField.isovar.field());

  double overallMinEdgeLength = std::numeric_limits<double>::max();
  double overallMaxEdgeLength = -std::numeric_limits<double>::max();
  double overallMinVolume = std::numeric_limits<double>::max();
  double overallMaxVolume = -std::numeric_limits<double>::max();
  const double expectedMinLength = mesh_size * approxMinRelativeSize;
  const double expectedMinVol = std::pow(mesh_size * approxMinRelativeSize, 2.) / 2.;

  std::mt19937 mt(std::mt19937::default_seed + parallel_rank);
  std::uniform_real_distribution<double> dist(-1., 1.);
#ifdef NDEBUG
  const int num_cases = 5000;
#else
  const int num_cases = 1000;
#endif
  for (int i = 0; i < num_cases; ++i)
  {
    if (i % 1000 == 0) std::cout << "Testing random configuration " << i << std::endl;

    MeshClone::stash_or_restore_mesh(mesh, 0); // restore original uncut mesh
    commit();                                  // new krino_mesh each time

    for (auto && lsField : lsFields)
    {
      randomize_ls_field(mesh, lsField.isovar, mt, dist);
    }
    stk::mesh::communicate_field_data(mesh, sync_fields);

    try
    {
      decompose_mesh();
    }
    catch (const std::exception & exception)
    {
      std::cout << log.get_log() << std::endl;
      std::cout << "Decomposing mesh failed with exception:\n";
      std::cout << exception.what() << "\n";
      ASSERT_TRUE(false);
    }

    EXPECT_TRUE(check_induced_parts(mesh));
    EXPECT_TRUE(check_face_and_edge_ownership(mesh));
    EXPECT_TRUE(check_face_and_edge_relations(mesh));
    EXPECT_TRUE(check_shared_entity_nodes(mesh));
    EXPECT_TRUE(krino_mesh->check_element_side_parts());

    double minEdgeLength, maxEdgeLength, minVolume, maxVolume;
    compute_element_quality(mesh, minEdgeLength, maxEdgeLength, minVolume, maxVolume);
    overallMinEdgeLength = std::min(overallMinEdgeLength, minEdgeLength);
    overallMaxEdgeLength = std::max(overallMaxEdgeLength, maxEdgeLength);
    overallMinVolume = std::min(overallMinVolume, minVolume);
    overallMaxVolume = std::max(overallMaxVolume, maxVolume);

    if (HasNonfatalFailure())
    {
      std::cout << "Failure on iteration " << i << std::endl;
      krino_mesh->debug_output();
      std::cout << log.get_log() << std::endl;
      std::ostringstream fname;
      fname << "Random_SnapTri3_iter_" << i << ".e";
      SimpleStkFixture::write_results(fname.str(), mesh);
      ASSERT_TRUE(false);
    }
  }
  std::cout << "Expected quality: minEdgeLength~=" << expectedMinLength
            << ", minVolume~=" << expectedMinVol << std::endl;
  std::cout << "Quality results: minEdgeLength=" << overallMinEdgeLength
            << ", maxEdgeLength=" << overallMaxEdgeLength << ", minVolume=" << overallMinVolume
            << ", maxVolume=" << overallMaxVolume << std::endl;
}

typedef CompleteDecompositionFixture<BoundingBoxMeshTet4, LSPerInterfacePolicy, 1> CDMeshTestsBboxMesh3D;
TEST_F(CDMeshTestsBboxMesh3D, Random_SnapMesh)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_size = stk::parallel_machine_size(pm);
  const int parallel_rank = stk::parallel_machine_rank(pm);

  stk::mesh::MetaData & meta = fixture.meta_data();
  AuxMetaData & aux_meta = AuxMetaData::get(meta);

  if (parallel_size > 2) return;

  cdfemSupport.set_cdfem_edge_degeneracy_handling(
      SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE);
  const double approxMinRelativeSize = 0.25;

  auto & block1_part = aux_meta.get_part("block_1");

  register_ls_on_blocks({&block1_part});

  typename BoundingBoxMesh::BoundingBoxType domain(stk::math::Vector3d::ZERO, stk::math::Vector3d(1., 1., 1.));
  const double mesh_size = 1. / 3.;
  fixture.set_domain(domain, mesh_size);
  fixture.populate_mesh();
  stk::mesh::BulkData & mesh = fixture.bulk_data();
  stk::mesh::create_exposed_block_boundary_sides(
      mesh, meta.universal_part(), {&aux_meta.exposed_boundary_part(), &aux_meta.active_part()});

  double overallMinEdgeLength = std::numeric_limits<double>::max();
  double overallMaxEdgeLength = -std::numeric_limits<double>::max();
  double overallMinVolume = std::numeric_limits<double>::max();
  double overallMaxVolume = -std::numeric_limits<double>::max();
  const double expectedMinLength = mesh_size * approxMinRelativeSize;
  const double expectedMinVol = std::pow(mesh_size * approxMinRelativeSize, 3.) / 6.;

  std::mt19937 mt(std::mt19937::default_seed + parallel_rank);
  std::uniform_real_distribution<double> dist(-1., 1.);
#ifdef NDEBUG
  const int num_cases = 1000;
#else
  const int num_cases = 250;
#endif
  for (int i = 0; i < num_cases; ++i)
  {
    if (i % 250 == 0) std::cout << "Testing random configuration " << i << std::endl;

    MeshClone::stash_or_restore_mesh(mesh, 0); // restore original uncut mesh
    commit();                                  // new krino_mesh each time

    randomize_ls_field(mesh, get_ls_field(), mt, dist);
    set_ls_field_on_part(mesh, aux_meta.exposed_boundary_part(), get_ls_field(), 1.);
    stk::mesh::communicate_field_data(mesh, {&get_ls_field().field(), &coord_field.field()});

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

    EXPECT_TRUE(check_induced_parts(mesh));
    EXPECT_TRUE(check_face_and_edge_ownership(mesh));
    EXPECT_TRUE(check_face_and_edge_relations(mesh));
    EXPECT_TRUE(check_shared_entity_nodes(mesh));
    EXPECT_TRUE(krino_mesh->check_element_side_parts());

    double minEdgeLength, maxEdgeLength, minVolume, maxVolume;
    compute_element_quality(mesh, minEdgeLength, maxEdgeLength, minVolume, maxVolume);
    overallMinEdgeLength = std::min(overallMinEdgeLength, minEdgeLength);
    overallMaxEdgeLength = std::max(overallMaxEdgeLength, maxEdgeLength);
    overallMinVolume = std::min(overallMinVolume, minVolume);
    overallMaxVolume = std::max(overallMaxVolume, maxVolume);

    bool failedQuality = false;
    if (minVolume < 0.5 * expectedMinVol)
    {
      failedQuality = true;
      std::cout << "Failed quality requirements: minEdgeLength=" << minEdgeLength
                << ", maxEdgeLength=" << maxEdgeLength << ", minVolume=" << minVolume
                << ", maxVolume=" << maxVolume << std::endl;
    }

    if (HasNonfatalFailure() || failedQuality)
    {
      std::cout << "Failure on iteration " << i << std::endl;
      krino_mesh->debug_output();
      std::cout << log.get_log() << std::endl;
      std::ostringstream fname;
      fname << "Random_SnapTet4_iter_" << i << ".e";
      SimpleStkFixture::write_results(fname.str(), mesh);
      ASSERT_TRUE(false);
    }
  }
  std::cout << "Expected quality: minEdgeLength~=" << expectedMinLength
            << ", minVolume~=" << expectedMinVol << std::endl;
  std::cout << "Actual quality results: minEdgeLength=" << overallMinEdgeLength
            << ", maxEdgeLength=" << overallMaxEdgeLength << ", minVolume=" << overallMinVolume
            << ", maxVolume=" << overallMaxVolume << std::endl;
}

typedef CompleteDecompositionFixture<BoundingBoxMeshTet4, LSPerInterfacePolicy, 1> CDMeshTestsBboxMesh3DBCC;
TEST_F(CDMeshTestsBboxMesh3DBCC, Random_SnapMesh)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_size = stk::parallel_machine_size(pm);
  const int parallel_rank = stk::parallel_machine_rank(pm);

  stk::mesh::MetaData & meta = fixture.meta_data();
  AuxMetaData & aux_meta = AuxMetaData::get(meta);

  if (parallel_size > 2) return;

  cdfemSupport.set_cdfem_edge_degeneracy_handling(
      SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE);
  const double approxMinRelativeSize = 0.25;

  auto & block1_part = aux_meta.get_part("block_1");

  register_ls_on_blocks({&block1_part});

  typename BoundingBoxMesh::BoundingBoxType domain(stk::math::Vector3d::ZERO, stk::math::Vector3d(1., 1., 1.));
  const double mesh_size = 1. / 3.;
  fixture.set_domain(domain, mesh_size);
  fixture.set_mesh_structure_type(BCC_BOUNDING_BOX_MESH);
  fixture.populate_mesh();
  stk::mesh::BulkData & mesh = fixture.bulk_data();
  stk::mesh::create_exposed_block_boundary_sides(
      mesh, meta.universal_part(), {&aux_meta.exposed_boundary_part(), &aux_meta.active_part()});

  const double BCCsize = std::sqrt(3.) / 2. * mesh_size;
  double overallMinEdgeLength = std::numeric_limits<double>::max();
  double overallMaxEdgeLength = -std::numeric_limits<double>::max();
  double overallMinVolume = std::numeric_limits<double>::max();
  double overallMaxVolume = -std::numeric_limits<double>::max();
  const double expectedMinLength = BCCsize * approxMinRelativeSize;
  const double expectedMinVol =
      std::pow(BCCsize * approxMinRelativeSize, 3.) / (6. * std::sqrt(2.));

  std::mt19937 mt(std::mt19937::default_seed + parallel_rank);
  std::uniform_real_distribution<double> dist(-1., 1.);
#ifdef NDEBUG
  const int num_cases = 1000;
#else
  const int num_cases = 250;
#endif
  for (int i = 0; i < num_cases; ++i)
  {
    if (i % 250 == 0) std::cout << "Testing random configuration " << i << std::endl;

    MeshClone::stash_or_restore_mesh(mesh, 0); // restore original uncut mesh
    commit();                                  // new krino_mesh each time

    randomize_ls_field(mesh, get_ls_field(), mt, dist);
    set_ls_field_on_part(mesh, aux_meta.exposed_boundary_part(), get_ls_field(), 1.);
    stk::mesh::communicate_field_data(mesh, {&get_ls_field().field(), &coord_field.field()});

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

    EXPECT_TRUE(check_induced_parts(mesh));
    EXPECT_TRUE(check_face_and_edge_ownership(mesh));
    EXPECT_TRUE(check_face_and_edge_relations(mesh));
    EXPECT_TRUE(check_shared_entity_nodes(mesh));
    EXPECT_TRUE(krino_mesh->check_element_side_parts());

    double minEdgeLength, maxEdgeLength, minVolume, maxVolume;
    compute_element_quality(mesh, minEdgeLength, maxEdgeLength, minVolume, maxVolume);
    overallMinEdgeLength = std::min(overallMinEdgeLength, minEdgeLength);
    overallMaxEdgeLength = std::max(overallMaxEdgeLength, maxEdgeLength);
    overallMinVolume = std::min(overallMinVolume, minVolume);
    overallMaxVolume = std::max(overallMaxVolume, maxVolume);

    bool failedQuality = false;
    if (minVolume < 0.5 * expectedMinVol)
    {
      failedQuality = true;
      std::cout << "Failed quality requirements: minEdgeLength=" << minEdgeLength
                << ", maxEdgeLength=" << maxEdgeLength << ", minVolume=" << minVolume
                << ", maxVolume=" << maxVolume << std::endl;
    }

    if (HasNonfatalFailure() || failedQuality)
    {
      std::cout << "Failure on iteration " << i << std::endl;
      krino_mesh->debug_output();
      std::cout << log.get_log() << std::endl;
      std::ostringstream fname;
      fname << "Random_SnapTet4_iter_" << i << ".e";
      SimpleStkFixture::write_results(fname.str(), mesh);
      ASSERT_TRUE(false);
    }
  }
  std::cout << "Expected quality: minEdgeLength~=" << expectedMinLength
            << ", minVolume~=" << expectedMinVol << std::endl;
  std::cout << "Actual quality results: minEdgeLength=" << overallMinEdgeLength
            << ", maxEdgeLength=" << overallMaxEdgeLength << ", minVolume=" << overallMinVolume
            << ", maxVolume=" << overallMaxVolume << std::endl;
}

class CDMeshTestsRightTri3LSPerPhase : public DecompositionFixture<RightTri, LSPerPhasePolicy, 3>
{
public:
CDMeshTestsRightTri3LSPerPhase()
{
  if(stk::parallel_machine_size(mComm) == 1)
    this->build_mesh(meshSpec.nodeLocs, {meshSpec.allElementConn});

  setup_ls_fields();
}
};

TEST_F(CDMeshTestsRightTri3LSPerPhase, Tri3_3LS_SnappedTriplePoint)
{
  if (stk::parallel_machine_size(mComm) > 1) return;

  cdfem_support().set_cdfem_edge_tol(0.15);

  set_level_sets(mMesh, levelset_fields(),
      get_ids_of_nodes_with_given_indices({0, 1, 2}),
      {{0.0, 0.1, 0.2}, {0.0, -0.2, -0.25}, {0.0, -0.01, -0.005}});

  attempt_decompose_mesh();

  // This assumes that the element is cut first by the (0,1) interface,
  // then the (0, 2) virtual interface, and then the (1, 2) interface.
  // THere are other valid decompositions if the element is cut in a different order.
  // Should have added 3 nodes at the cutting locations
  expect_num_nodes(5u);
  expect_num_nodes(5u, get_aux_meta().active_part());

  expect_num_sides(1u, get_aux_meta().active_part() & get_aux_meta().get_part("surface_block_1_P0_P1") & get_aux_meta().get_part("surface_block_1_P1_P0"));
  expect_num_sides(1u, get_aux_meta().active_part() & get_aux_meta().get_part("surface_block_1_P1_P2") & get_aux_meta().get_part("surface_block_1_P2_P1"));
  expect_num_sides(0u, get_aux_meta().active_part() & get_aux_meta().get_part("surface_block_1_P0_P2") & get_aux_meta().get_part("surface_block_1_P2_P0"));

  // Should be 3 conformal elements plus the parent element
  expect_num_elements(4u);
  expect_num_elements(3u, get_aux_meta().active_part());
  expect_num_elements(1u, get_aux_meta().get_part("block_1_P0"));
  expect_num_elements(1u, get_aux_meta().get_part("block_1_P1"));
  expect_num_elements(1u, get_aux_meta().get_part("block_1_P2"));
}

TEST_F(CDMeshTestsRightTri3LSPerPhase, Tri3_3LS_TriplePointDebug)
{
  if (stk::parallel_machine_size(mComm) > 1) return;

  cdfem_support().set_cdfem_edge_tol(0.01);

  set_level_sets(mMesh, levelset_fields(),
      get_ids_of_nodes_with_given_indices({0, 1, 2}),
      {{-0.2, 0.7, 0.1}, {0.1, -0.5, 0.3}, {0.6, 0.4, -0.5}});

  attempt_decompose_mesh();

  // Test output, but remove output unless actually debugging.
  write_mesh("debug_2d.e");
  std::remove("debug_2d.e");
}

class TwoRightTrisOn1Or2Procs3LSPerPhaseDecompositionFixture : public DecompositionFixture<TwoRightTrisSharingDiagonal, LSPerPhasePolicy, 3>
{
public:
TwoRightTrisOn1Or2Procs3LSPerPhaseDecompositionFixture()
{
  if(stk::parallel_machine_size(mComm) == 1)
    this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,2}, {0,0});
  else if(stk::parallel_machine_size(mComm) == 2)
    this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,2}, {0,1});

  setup_ls_fields();
}
};

TEST_F(TwoRightTrisOn1Or2Procs3LSPerPhaseDecompositionFixture, Random_TwoTri3_InternalSideset_Snap)
{
  if (stk::parallel_machine_size(mComm) > 2) return;

  cdfem_support().set_cdfem_edge_degeneracy_handling(SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE);

  const auto & lsFields = levelset_fields();
  std::vector<const stk::mesh::FieldBase *> sync_fields = {&get_coordinates_field().field()};
  for (auto && ls_field : lsFields)
    sync_fields.push_back(&ls_field.isovar.field());

  std::mt19937 mt(std::mt19937::default_seed + stk::parallel_machine_rank(mComm));
  std::uniform_real_distribution<double> dist(-1., 1.);
#ifdef NDEBUG
  const int num_cases = 5000;
#else
  const int num_cases = 1000;
#endif
  for (int i = 0; i < num_cases; ++i)
  {
    if (i % 1000 == 0) std::cout << "Testing random configuration " << i << std::endl;

    MeshClone::stash_or_restore_mesh(mMesh, 0); // restore original uncut mesh
    cdmesh = std::make_unique<CDMesh>(mMesh);   // new krino_mesh each time

    for (auto && ls_field : lsFields)
    {
      randomize_ls_field(mMesh, ls_field.isovar, mt, dist);
    }
    stk::mesh::communicate_field_data(mMesh, sync_fields);

    attempt_decompose_mesh();

    check_mesh_consistency();

    check_nonfatal_error("Random_TwoTri3_InternalSideset_Snap_iter", i);

    cdmesh = std::make_unique<CDMesh>(mMesh);
  }
}

typedef CompleteDecompositionFixture<SimpleStkFixture3d, LSPerInterfacePolicy, 3>
    CDMeshTests3DLSPerInterface;
TEST_F(CDMeshTests3DLSPerInterface, Random_OneTet4)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_size = stk::parallel_machine_size(pm);
  const int parallel_rank = stk::parallel_machine_rank(pm);

  if (parallel_size > 2) return;

  stk::mesh::BulkData & mesh = fixture.bulk_data();
  stk::mesh::MetaData & meta = fixture.meta_data();
  AuxMetaData & aux_meta = AuxMetaData::get(meta);

  aux_meta.set_assert_32bit_flag();
  aux_meta.clear_force_64bit_flag();

  // Use a large snap tolerance to make snapping more common since it is a frequent source
  // of parallel bugs
  cdfemSupport.set_cdfem_edge_tol(0.1);
  cdfemSupport.set_simplex_generation_method(CUT_QUADS_BY_GLOBAL_IDENTIFIER);

  const stk::topology tet4 = stk::topology::TETRAHEDRON_4;
  auto & block1_part = declare_input_block("block_1", tet4);

  register_ls_on_blocks({&block1_part});

  commit();

  build_one_tet4_mesh(*this, block1_part, parallel_rank, parallel_size);

  stk::mesh::create_exposed_block_boundary_sides(
      mesh, meta.universal_part(), {&aux_meta.exposed_boundary_part()});

  const auto & ls_fields = levelset_fields();
  std::vector<const stk::mesh::FieldBase *> sync_fields = {&coord_field.field()};
  for (auto && ls_field : ls_fields)
    sync_fields.push_back(&ls_field.isovar.field());

  std::mt19937 mt(std::mt19937::default_seed + parallel_rank);
  std::uniform_real_distribution<double> dist(-1., 1.);
#ifdef NDEBUG
  const int num_cases = 5000;
#else
  const int num_cases = 1000;
#endif
  for (int i = 0; i < num_cases; ++i)
  {
    if (i % 1000 == 0) std::cout << "Testing random configuration " << i << std::endl;

    MeshClone::stash_or_restore_mesh(mesh, 0); // restore original uncut mesh
    commit();                                  // new krino_mesh each time

    for (auto && ls_field : ls_fields)
    {
      randomize_ls_field(mesh, ls_field.isovar, mt, dist);
    }
    stk::mesh::communicate_field_data(mesh, sync_fields);

    try
    {
      decompose_mesh();
    }
    catch (const std::exception & exception)
    {
      std::cout << log.get_log() << std::endl;
      std::cout << "Decomposing mesh failed with exception:\n";
      std::cout << exception.what() << "\n";
      ASSERT_TRUE(false);
    }

    EXPECT_TRUE(check_induced_parts(mesh));
    EXPECT_TRUE(check_face_and_edge_ownership(mesh));
    EXPECT_TRUE(check_face_and_edge_relations(mesh));
    EXPECT_TRUE(check_shared_entity_nodes(mesh));
    EXPECT_TRUE(krino_mesh->check_element_side_parts());

    if (HasNonfatalFailure())
    {
      std::cout << log.get_log() << std::endl;
      std::cout << "Failure on iteration " << i << std::endl;
      std::ostringstream fname;
      fname << "Random_TwoTet4_InternalSideset_iter_" << i << ".e";
      fixture.write_results(fname.str());
      ASSERT_TRUE(false);
    }

    krino_mesh = std::make_unique<CDMesh>(mesh);
  }
}

typedef CompleteDecompositionFixture<SimpleStkFixture3d, LSPerInterfacePolicy, 3>
    CDMeshTests3DLSPerInterface;
TEST_F(CDMeshTests3DLSPerInterface, OneTet4_CutBasedOnNearestEdgeCut)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_size = stk::parallel_machine_size(pm);
  const int parallel_rank = stk::parallel_machine_rank(pm);

  if (parallel_size > 1) return;

  stk::mesh::BulkData & mesh = fixture.bulk_data();
  stk::mesh::MetaData & meta = fixture.meta_data();
  AuxMetaData & aux_meta = AuxMetaData::get(meta);

  aux_meta.set_assert_32bit_flag();
  aux_meta.clear_force_64bit_flag();

  // Use a large snap tolerance to make snapping more common since it is a frequent source
  // of parallel bugs
  cdfemSupport.set_cdfem_edge_tol(0.001);
  cdfemSupport.set_simplex_generation_method(CUT_QUADS_BY_NEAREST_EDGE_CUT);

  const stk::topology tet4 = stk::topology::TETRAHEDRON_4;
  auto & block1_part = declare_input_block("block_1", tet4);

  register_ls_on_blocks({&block1_part});

  commit();

  build_one_tet4_mesh(*this, block1_part, parallel_rank, parallel_size);

  stk::mesh::create_exposed_block_boundary_sides(
      mesh, meta.universal_part(), {&aux_meta.exposed_boundary_part()});

  const std::array<stk::mesh::Entity, 4> nodes = {{mesh.get_entity(stk::topology::NODE_RANK, 1),
      mesh.get_entity(stk::topology::NODE_RANK, 2),
      mesh.get_entity(stk::topology::NODE_RANK, 3),
      mesh.get_entity(stk::topology::NODE_RANK, 4)}};
  const auto & ls_fields = levelset_fields();
  *field_data<double>(ls_fields[0].isovar, nodes[0]) = -1;
  *field_data<double>(ls_fields[0].isovar, nodes[1]) = 1.;
  *field_data<double>(ls_fields[0].isovar, nodes[2]) = 2.;
  *field_data<double>(ls_fields[0].isovar, nodes[3]) = 3.;

  commit();
  decompose_mesh();

  EXPECT_TRUE(check_induced_parts(mesh));
  EXPECT_TRUE(check_face_and_edge_ownership(mesh));
  EXPECT_TRUE(check_face_and_edge_relations(mesh));
  EXPECT_TRUE(check_shared_entity_nodes(mesh));
  EXPECT_TRUE(krino_mesh->check_element_side_parts());

  EXPECT_EQ(1u + 1u, mesh.num_elements(nodes[0]));
  EXPECT_EQ(1u + 1u, mesh.num_elements(nodes[1]));
  EXPECT_EQ(2u + 1u, mesh.num_elements(nodes[2]));
  EXPECT_EQ(3u + 1u, mesh.num_elements(nodes[3]));

  // regression test
  const ScaledJacobianQualityMetric qualityMetric;
  const double quality = compute_mesh_quality(mesh, krino_mesh->get_active_part(), qualityMetric);
  const double goldQuality = 0.21;
  EXPECT_GT(quality, goldQuality);
}

typedef CompleteDecompositionFixture<SimpleStkFixture3d, LSPerPhasePolicy, 3>
    CDMeshTests3DLSPerPhase;
TEST_F(CDMeshTests3DLSPerPhase, Random_TwoTet4_InternalSideset)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_size = stk::parallel_machine_size(pm);
  const int parallel_rank = stk::parallel_machine_rank(pm);

  if (parallel_size > 2) return;

  stk::mesh::BulkData & mesh = fixture.bulk_data();
  stk::mesh::MetaData & meta = fixture.meta_data();
  AuxMetaData & aux_meta = AuxMetaData::get(meta);

  // Use a large snap tolerance to make snapping more common since it is a frequent source
  // of parallel bugs
  cdfemSupport.set_cdfem_edge_tol(0.1);
  cdfemSupport.set_simplex_generation_method(CUT_QUADS_BY_GLOBAL_IDENTIFIER);

  const stk::topology tet4 = stk::topology::TETRAHEDRON_4;
  auto & block1_part = declare_input_block("block_1", tet4);
  auto & block2_part = declare_input_block("block_2", tet4);
  auto & surface_part = declare_input_surface("surface_1",
      tet4.side_topology(),
      {block1_part.mesh_meta_data_ordinal(), block2_part.mesh_meta_data_ordinal()});

  register_ls_on_blocks({&block1_part, &block2_part});

  commit();

  build_two_tet4_mesh_np2(
      *this, block1_part, block2_part, surface_part, parallel_rank, parallel_size);

  stk::mesh::create_exposed_block_boundary_sides(
      mesh, meta.universal_part(), {&aux_meta.exposed_boundary_part()});

  const auto & ls_fields = levelset_fields();
  std::vector<const stk::mesh::FieldBase *> sync_fields = {&coord_field.field()};
  for (auto && ls_field : ls_fields)
    sync_fields.push_back(&ls_field.isovar.field());

  std::mt19937 mt(std::mt19937::default_seed + parallel_rank);
  std::uniform_real_distribution<double> dist(-1., 1.);
#ifdef NDEBUG
  const int num_cases = 5000;
#else
  const int num_cases = 1000;
#endif
  for (int i = 0; i < num_cases; ++i)
  {
    if (i % 1000 == 0) std::cout << "Testing random configuration " << i << std::endl;

    for (auto && ls_field : ls_fields)
    {
      randomize_ls_field(mesh, ls_field.isovar, mt, dist);
    }
    stk::mesh::communicate_field_data(mesh, sync_fields);

    try
    {
      decompose_mesh();
    }
    catch (const std::exception & exception)
    {
      std::cout << log.get_log() << std::endl;
      std::cout << "Decomposing mesh failed with exception:\n";
      std::cout << exception.what() << "\n";
      ASSERT_TRUE(false);
    }

    EXPECT_TRUE(check_induced_parts(mesh));
    EXPECT_TRUE(check_face_and_edge_ownership(mesh));
    EXPECT_TRUE(check_face_and_edge_relations(mesh));
    EXPECT_TRUE(check_shared_entity_nodes(mesh));
    EXPECT_TRUE(krino_mesh->check_element_side_parts());

    if (HasNonfatalFailure())
    {
      std::cout << "Failure on iteration " << i << std::endl;
      std::ostringstream fname;
      fname << "Random_TwoTet4_InternalSideset_iter_" << i << ".e";
      fixture.write_results(fname.str());
      ASSERT_TRUE(false);
    }

    krino_mesh = std::make_unique<CDMesh>(mesh);
  }
}

TEST_F(CDMeshTests3DLSPerPhase, Random_TwoTet4_InternalSideset_Snap)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_size = stk::parallel_machine_size(pm);
  const int parallel_rank = stk::parallel_machine_rank(pm);

  if (parallel_size > 2) return;

  stk::mesh::BulkData & mesh = fixture.bulk_data();
  stk::mesh::MetaData & meta = fixture.meta_data();
  AuxMetaData & aux_meta = AuxMetaData::get(meta);

  cdfemSupport.set_cdfem_edge_degeneracy_handling(
      SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE);

  const stk::topology tet4 = stk::topology::TETRAHEDRON_4;
  auto & block1_part = declare_input_block("block_1", tet4);
  auto & block2_part = declare_input_block("block_2", tet4);
  auto & surface_part = declare_input_surface("surface_1",
      tet4.side_topology(),
      {block1_part.mesh_meta_data_ordinal(), block2_part.mesh_meta_data_ordinal()});

  register_ls_on_blocks({&block1_part, &block2_part});

  commit();

  build_two_tet4_mesh_np2(
      *this, block1_part, block2_part, surface_part, parallel_rank, parallel_size);

  stk::mesh::create_exposed_block_boundary_sides(
      mesh, meta.universal_part(), {&aux_meta.exposed_boundary_part()});

  const auto & ls_fields = levelset_fields();
  std::vector<const stk::mesh::FieldBase *> sync_fields = {&coord_field.field()};
  for (auto && ls_field : ls_fields)
    sync_fields.push_back(&ls_field.isovar.field());

  std::mt19937 mt(std::mt19937::default_seed + parallel_rank);
  std::uniform_real_distribution<double> dist(-1., 1.);
#ifdef NDEBUG
  const int num_cases = 5000;
#else
  const int num_cases = 1000;
#endif
  for (int i = 0; i < num_cases; ++i)
  {
    if (i % 1000 == 0) std::cout << "Testing random configuration " << i << std::endl;

    MeshClone::stash_or_restore_mesh(mesh, 0); // restore original uncut mesh
    commit();                                  // new krino_mesh each time

    for (auto && ls_field : ls_fields)
    {
      randomize_ls_field(mesh, ls_field.isovar, mt, dist);
    }
    stk::mesh::communicate_field_data(mesh, sync_fields);

    try
    {
      decompose_mesh();
    }
    catch (const std::exception & exception)
    {
      std::cout << log.get_log() << std::endl;
      std::cout << "Decomposing mesh failed with exception:\n";
      std::cout << exception.what() << "\n";
      ASSERT_TRUE(false);
    }

    EXPECT_TRUE(check_induced_parts(mesh));
    EXPECT_TRUE(check_face_and_edge_ownership(mesh));
    EXPECT_TRUE(check_face_and_edge_relations(mesh));
    EXPECT_TRUE(check_shared_entity_nodes(mesh));
    EXPECT_TRUE(krino_mesh->check_element_side_parts());

    if (HasNonfatalFailure())
    {
      std::cout << "Failure on iteration " << i << std::endl;
      std::ostringstream fname;
      fname << "Random_TwoTet4_InternalSideset_iter_" << i << ".e";
      fixture.write_results(fname.str());
      ASSERT_TRUE(false);
    }

    krino_mesh = std::make_unique<CDMesh>(mesh);
  }
}

TEST_F(CDMeshTests3DLSPerPhase, RestoreAfterFailedStep)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_size = stk::parallel_machine_size(pm);
  const int parallel_rank = stk::parallel_machine_rank(pm);

  if (parallel_size > 2) return;

  stk::mesh::BulkData & mesh = fixture.bulk_data();
  stk::mesh::MetaData & meta = fixture.meta_data();
  AuxMetaData & aux_meta = AuxMetaData::get(meta);

  // Use a large snap tolerance to make snapping more common since it is a frequent source
  // of parallel bugs
  cdfemSupport.set_cdfem_edge_tol(0.1);
  cdfemSupport.set_simplex_generation_method(CUT_QUADS_BY_GLOBAL_IDENTIFIER);

  const stk::topology tet4 = stk::topology::TETRAHEDRON_4;
  auto & block1_part = declare_input_block("block_1", tet4);
  auto & block2_part = declare_input_block("block_2", tet4);
  auto & surface_part = declare_input_surface("surface_1",
      tet4.side_topology(),
      {block1_part.mesh_meta_data_ordinal(), block2_part.mesh_meta_data_ordinal()});

  register_ls_on_blocks({&block1_part, &block2_part});

  commit();

  build_two_tet4_mesh_np2(
      *this, block1_part, block2_part, surface_part, parallel_rank, parallel_size);

  stk::mesh::create_exposed_block_boundary_sides(
      mesh, meta.universal_part(), {&aux_meta.exposed_boundary_part()});

  const auto & ls_fields = levelset_fields();
  std::vector<const stk::mesh::FieldBase *> sync_fields = {&coord_field.field()};
  for (auto && ls_field : ls_fields)
    sync_fields.push_back(&ls_field.isovar.field());

  std::mt19937 mt(std::mt19937::default_seed + parallel_rank);
  std::uniform_real_distribution<double> dist(-1., 1.);
  for (int i = 0; i < 100; ++i)
  {
    if (i % 1000 == 0) std::cout << "Testing random configuration " << i << std::endl;

    // We want to test restoring the mesh after failed time steps for a variety of
    // random decompositions:
    // 1) Do a random decomposition and stash the mesh to act like a successful time
    //    step.
    // 2) Do a second random decomposition to simulate the changes of a time step that
    //    will fail.
    // 3) Restore the mesh to the stashed mesh, and rebuild the CDMesh, confirm that
    //    does not fail.

    auto randomize_all_fields = [&]()
    {
      for (auto && ls_field : ls_fields)
      {
        randomize_ls_field(mesh, ls_field.isovar, mt, dist);
      }
      stk::mesh::communicate_field_data(mesh, sync_fields);
    };

    auto run_decomp = [&]()
    {
      try
      {
        decompose_mesh();
      }
      catch (const std::exception & exception)
      {
        std::cout << log.get_log() << std::endl;
        std::cout << "Decomposing mesh failed with exception:\n";
        std::cout << exception.what() << "\n";
        ASSERT_TRUE(false);
      }
    };
    // Step 1.
    randomize_all_fields();
    run_decomp();
    krino::MeshClone::stash_or_restore_mesh(mesh, i);
    krino_mesh = std::make_unique<CDMesh>(mesh);

    // Step 2.
    randomize_all_fields();
    run_decomp();
    krino::MeshClone::stash_or_restore_mesh(mesh, i);

    // Step 3.
    ASSERT_NO_THROW(krino_mesh->rebuild_after_rebalance_or_failed_step());
    krino_mesh = std::make_unique<CDMesh>(mesh);
  }
}

typedef CompleteDecompositionFixture<SimpleStkFixture3d, LSPerInterfacePolicy, 1> CDMeshTests3D;

TEST_F(CDMeshTests3D, Rebalance_with_rcb)
{
  run_rebalance_with("rcb");
}

TEST_F(CDMeshTests3D, Rebalance_with_parmetis)
{
  if (rebalance_utils::have_parmetis())
    run_rebalance_with("parmetis");
}

typedef CompleteDecompositionFixture<SimpleStkFixture3d, LSPerInterfacePolicy, 1> MeshCloneTest;
TEST_F(MeshCloneTest, FaceOwnershipAndPartChangeBetweenClones)
{
  /*
   * This test regression tests a bug where a shared face both changed parts and
   * parallel owners between calls to MeshClone::stash_or_restore_mesh(), leading to a
   * throw when copying field data to the clone mesh because the parts of the face were not
   * updated on the clone.
   */
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_size = stk::parallel_machine_size(pm);
  const int parallel_rank = stk::parallel_machine_rank(pm);
  stk::mesh::BulkData & mesh = fixture.bulk_data();

  if (parallel_size != 2) return;

  const stk::topology tet4 = stk::topology::TETRAHEDRON_4;
  auto & block1_part = declare_input_block("block_1", tet4);
  auto & block2_part = declare_input_block("block_2", tet4);
  auto & surface_1_part = declare_input_surface("surface_1",
      tet4.side_topology(),
      {block1_part.mesh_meta_data_ordinal(), block2_part.mesh_meta_data_ordinal()});
  auto & surface_2_part = declare_input_surface("surface_2",
      tet4.side_topology(),
      {block1_part.mesh_meta_data_ordinal(), block2_part.mesh_meta_data_ordinal()});

  auto & aux_meta = AuxMetaData::get(fixture.meta_data());
  aux_meta.register_field(
      "side_field", FieldType::REAL, stk::topology::FACE_RANK, 1u, 1u, surface_1_part);

  commit();

  build_two_tet4_mesh_np2(
      *this, block1_part, block2_part, surface_1_part, parallel_rank, parallel_size);

  // Stash the initial mesh.
  ASSERT_NO_THROW(MeshClone::stash_or_restore_mesh(mesh, 0));

  // Change the parallel owner of the shared face, then change its surface part from surface_1 to
  // surface_2.
  const auto side_1 = mesh.get_entity(stk::topology::FACE_RANK, 7);
  const auto current_owner = mesh.parallel_owner_rank(side_1);
  stk::mesh::EntityProcVec owner_changes;
  if (parallel_rank == current_owner)
  {
    owner_changes.emplace_back(side_1, (current_owner + 1) % 2);
  }
  mesh.change_entity_owner(owner_changes);
  ASSERT_TRUE(mesh.is_valid(side_1));
  const auto new_owner = mesh.parallel_owner_rank(side_1);
  ASSERT_NE(new_owner, current_owner);

  mesh.modification_begin();
  if (new_owner == parallel_rank)
  {
    mesh.change_entity_parts(side_1,
        stk::mesh::ConstPartVector{&surface_2_part},
        stk::mesh::ConstPartVector{&surface_1_part});
  }
  mesh.modification_end();

  // Confirm that updating the stash mesh succeeds
  try
  {
    MeshClone::stash_or_restore_mesh(mesh, 1);
  }
  catch (const std::exception & exception)
  {
    std::cout << "Decomposing mesh failed with exception:\n";
    std::cout << exception.what() << "\n";
    ASSERT_TRUE(false);
  }

  // Change the shared face part back to surface_1, then restore from the stash and
  // confirm that the parts of the face are correct.
  mesh.modification_begin();
  if (new_owner == parallel_rank)
  {
    mesh.change_entity_parts(side_1,
        stk::mesh::ConstPartVector{&surface_1_part},
        stk::mesh::ConstPartVector{&surface_2_part});
  }
  mesh.modification_end();
  ASSERT_TRUE(mesh.bucket(side_1).member(surface_1_part));
  MeshClone::stash_or_restore_mesh(mesh, 1);
  const auto side_1_new = mesh.get_entity(stk::topology::FACE_RANK, 7);
  EXPECT_TRUE(mesh.bucket(side_1_new).member(surface_2_part));
  EXPECT_FALSE(mesh.bucket(side_1_new).member(surface_1_part));
}

} // namespace krino
