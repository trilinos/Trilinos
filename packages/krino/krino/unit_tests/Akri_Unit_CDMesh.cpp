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

#include <Akri_AdaptivityInterface.hpp>
#include <Akri_BoundingBox.hpp>
#include <Akri_BoundingBoxMesh.hpp>
#include <Akri_CDMesh.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_Interface_Name_Generator.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_MeshClone.hpp>
#include <Akri_Phase_Support.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_LevelSetInterfaceGeometry.hpp>
#include <Akri_QualityMetric.hpp>
#include <Akri_RebalanceUtils.hpp>
#include <Akri_Snap.hpp>

#include <Akri_Unit_Single_Element_Fixtures.hpp>
#include <Akri_Unit_LogRedirecter.hpp>

namespace krino {

class LevelSet;

namespace {
  template <class DECOMP_FIXTURE>
void build_one_tet4_mesh(DECOMP_FIXTURE & fixture,
    stk::mesh::Part & elem_part,
    const int parallel_rank, const int parallel_size)
{
  stk::mesh::BulkData & mesh = fixture.fixture.bulk_data();
  stk::mesh::MetaData & meta = fixture.fixture.meta_data();
  AuxMetaData & aux_meta = AuxMetaData::get(meta);

  ASSERT_TRUE(parallel_size == 1 || parallel_size == 2);
  ThrowRequire(elem_part.topology() == stk::topology::TETRAHEDRON_4);
  ASSERT_EQ(stk::topology::ELEMENT_RANK, elem_part.primary_entity_rank());

  mesh.modification_begin();
  {
    stk::mesh::PartVector elem_parts(2);
    elem_parts[1] = &aux_meta.active_part();

    std::vector<stk::mesh::EntityId> elem_nodes = {1, 2, 3, 4};

    stk::mesh::Entity element;
    if(parallel_rank == 0)
    {
      elem_parts[0] = &elem_part;
      element = fixture.create_element(elem_parts, 1, elem_nodes);
    }
  }
  mesh.modification_end();

  if(parallel_rank == 0)
  {
    EXPECT_EQ(1u, stk::mesh::count_selected_entities(meta.universal_part(), mesh.buckets(stk::topology::ELEMENT_RANK)));

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

    node1_coords[0] =  1.;
    node1_coords[1] =  1.;
    node1_coords[2] =  1.;

    node2_coords[0] = -1.;
    node2_coords[1] =  1.;
    node2_coords[2] = -1.;

    node3_coords[0] =  1.;
    node3_coords[1] = -1.;
    node3_coords[2] = -1.;

    node4_coords[0] = -1.;
    node4_coords[1] = -1.;
    node4_coords[2] =  1.;
  }
}

template <class DECOMP_FIXTURE>
void build_two_tet4_mesh_np2(DECOMP_FIXTURE & fixture,
    stk::mesh::Part & elem1_part, stk::mesh::Part & elem2_part, stk::mesh::Part & surface_part,
    const int parallel_rank, const int parallel_size, const bool add_side = true,
    const bool build_all_on_P0 = false)
{
  stk::mesh::BulkData & mesh = fixture.fixture.bulk_data();
  stk::mesh::MetaData & meta = fixture.fixture.meta_data();
  AuxMetaData & aux_meta = AuxMetaData::get(meta);

  ASSERT_TRUE(parallel_size == 1 || parallel_size == 2);
  ThrowRequire(elem1_part.topology() == stk::topology::TETRAHEDRON_4);
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
    if(parallel_rank == 0)
    {
      elem_parts[0] = &elem1_part;
      element = fixture.create_element(elem_parts, 1, elem1_nodes);
    }
    if((parallel_rank == 1 && !build_all_on_P0) || (parallel_rank == 0 && build_all_on_P0) ||
        parallel_size == 1)
    {
      elem_parts[0] = &elem2_part;
      element = fixture.create_element(elem_parts, 2, elem2_nodes);
    }

    if(parallel_size > 1 && !build_all_on_P0)
    {
      const int opp_rank = parallel_rank == 0 ? 1 : 0;
      for(auto i=1; i <= 3; ++i)
      {
        const auto node = mesh.get_entity(stk::topology::NODE_RANK, i);
        mesh.add_node_sharing(node, opp_rank);
      }
    }

    if(add_side && (parallel_rank == 0 || !build_all_on_P0))
    {
      sideEntity = mesh.declare_solo_side(7, {&surface_part});
      for(auto i=1; i <= 3; ++i)
      {
        auto side_node = mesh.get_entity(stk::topology::NODE_RANK, i);
        mesh.declare_relation(sideEntity, side_node, i-1);
      }
      attach_entity_to_elements(mesh, sideEntity);
    }
  }
  mesh.modification_end();

  if(parallel_rank == 0 || !build_all_on_P0)
  {
    EXPECT_EQ(2u, stk::mesh::count_selected_entities(meta.universal_part(), mesh.buckets(stk::topology::ELEMENT_RANK)));
    if(add_side)
    {
      EXPECT_EQ(1u, stk::mesh::count_selected_entities(meta.universal_part(), mesh.buckets(meta.side_rank())));
      EXPECT_EQ(1u, stk::mesh::count_selected_entities(surface_part, mesh.buckets(meta.side_rank())));

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

template <class DECOMP_FIXTURE>
void build_two_tri3_mesh_np2(DECOMP_FIXTURE & fixture,
    stk::mesh::Part & elem1_part, stk::mesh::Part & elem2_part, stk::mesh::Part & surface_part,
    const int parallel_rank, const int parallel_size, const bool add_side = true,
    const bool build_all_on_P0 = false)
{
  stk::mesh::BulkData & mesh = fixture.fixture.bulk_data();
  stk::mesh::MetaData & meta = fixture.fixture.meta_data();
  AuxMetaData & aux_meta = AuxMetaData::get(meta);

  ASSERT_TRUE(parallel_size == 1 || parallel_size == 2);
  ThrowRequire(elem1_part.topology() == stk::topology::TRIANGLE_3_2D);
  ASSERT_EQ(meta.side_rank(), surface_part.primary_entity_rank());
  ASSERT_EQ(stk::topology::ELEMENT_RANK, elem1_part.primary_entity_rank());
  ASSERT_EQ(stk::topology::ELEMENT_RANK, elem2_part.primary_entity_rank());

  stk::mesh::Entity sideEntity;
  mesh.modification_begin();
  {
    stk::mesh::PartVector elem_parts(2);
    elem_parts[1] = &aux_meta.active_part();

    std::vector<stk::mesh::EntityId> elem1_nodes = {2, 1, 3};
    std::vector<stk::mesh::EntityId> elem2_nodes = {1, 2, 4};

    stk::mesh::Entity element;
    if(parallel_rank == 0)
    {
      elem_parts[0] = &elem1_part;
      element = fixture.create_element(elem_parts, 1, elem1_nodes);
    }
    if((parallel_rank == 1 && !build_all_on_P0) || (parallel_rank == 0 && build_all_on_P0) ||
        parallel_size == 1)
    {
      elem_parts[0] = &elem2_part;
      element = fixture.create_element(elem_parts, 2, elem2_nodes);
    }

    if(parallel_size > 1 && !build_all_on_P0)
    {
      const int opp_rank = parallel_rank == 0 ? 1 : 0;
      for(auto i=1; i <= 2; ++i)
      {
        const auto node = mesh.get_entity(stk::topology::NODE_RANK, i);
        mesh.add_node_sharing(node, opp_rank);
      }
    }

    if(add_side && (parallel_rank == 0 || !build_all_on_P0))
    {
      sideEntity = mesh.declare_solo_side(7, {&surface_part});
      for(auto i=1; i <= 2; ++i)
      {
        auto side_node = mesh.get_entity(stk::topology::NODE_RANK, i);
        mesh.declare_relation(sideEntity, side_node, i-1);
      }
      attach_entity_to_elements(mesh, sideEntity);
    }
  }
  mesh.modification_end();

  if(parallel_rank == 0 || !build_all_on_P0)
  {
    EXPECT_EQ(2u, stk::mesh::count_selected_entities(meta.universal_part(), mesh.buckets(stk::topology::ELEMENT_RANK)));
    if(add_side)
    {
      EXPECT_EQ(1u, stk::mesh::count_selected_entities(meta.universal_part(), mesh.buckets(meta.side_rank())));
      EXPECT_EQ(1u, stk::mesh::count_selected_entities(surface_part, mesh.buckets(meta.side_rank())));

      EXPECT_EQ(2u, mesh.num_connectivity(sideEntity, stk::topology::ELEMENT_RANK));
    }

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

    node1_coords[0] = 0.;
    node1_coords[1] = 0.;

    node2_coords[0] = 0.;
    node2_coords[1] = 1.;

    node3_coords[0] = 1.;
    node3_coords[1] = 0.;

    node4_coords[0] = -1.;
    node4_coords[1] = 0.;
  }
}
}

struct SingleLSPolicy
{
  void setup_ls_field(const bool is_death, stk::mesh::MetaData & meta, CDFEM_Support & cdfem_support)
  {
    AuxMetaData & aux_meta = AuxMetaData::get(meta);
    Phase_Support::set_one_levelset_per_phase(false);
    ls_isovar = aux_meta.declare_field("LS", FieldType::REAL, stk::topology::NODE_RANK, 1u);
    cdfem_support.add_interpolation_field(ls_isovar);

    LevelSet * ls_ptr = nullptr;
    ls_field = std::make_shared<LS_Field>("LS", LevelSet_Identifier(0), ls_isovar, 0., ls_ptr);
    if(is_death)
    {
      death_spec = std::unique_ptr<CDFEM_Inequality_Spec>(new CDFEM_Inequality_Spec("death_spec"));
    }
    cdfem_support.add_ls_field(*ls_field, death_spec.get());
  }

  void register_ls_on_blocks(const stk::mesh::PartVector & blocks, stk::mesh::MetaData & meta,
      Block_Surface_Connectivity & block_surface_info, const bool register_fields)
  {
    AuxMetaData & aux_meta = AuxMetaData::get(meta);
    PhaseTag p, n;
    const LevelSet_Identifier id0(0);
    p.add(id0, 1);
    n.add(id0, -1);
    PhaseVec named_phases{{"A", p}, {"B", n}};

    if(death_spec)
    {
      death_spec->set_phases(p, n);
    }

    Phase_Support & phase_support = Phase_Support::get(meta);
    phase_support.set_input_block_surface_connectivity(block_surface_info);
    LevelSet * ls_ptr = nullptr;
    phase_support.register_blocks_for_level_set(ls_ptr, blocks);
    std::vector<std::tuple<stk::mesh::PartVector, 
      std::shared_ptr<Interface_Name_Generator>, PhaseVec>> ls_sets;
    auto interface_name_gen = std::shared_ptr<Interface_Name_Generator>(new LS_Name_Generator());
    ls_sets.push_back(std::make_tuple(blocks,interface_name_gen,named_phases));
    phase_support.decompose_blocks(ls_sets);

    for(auto && b : blocks)
    {
      auto conformal_A = phase_support.find_conformal_io_part(*b, p);
      auto conformal_B = phase_support.find_conformal_io_part(*b, n);
      ThrowRequire(conformal_A != nullptr);
      ThrowRequire(conformal_B != nullptr);
      if(register_fields)
      {
        aux_meta.register_field("LS", FieldType::REAL, stk::topology::NODE_RANK, 1u, 1u, *b);
        aux_meta.register_field("LS", FieldType::REAL, stk::topology::NODE_RANK, 1u, 1u, *conformal_A);
        aux_meta.register_field("LS", FieldType::REAL, stk::topology::NODE_RANK, 1u, 1u, *conformal_B);
      }
    }
  }

  FieldRef ls_isovar;
  std::shared_ptr<LS_Field> ls_field;
  std::unique_ptr<CDFEM_Inequality_Spec> death_spec;
};

template <unsigned NUM_LS>
struct LSPerPhasePolicy
{
  void setup_ls_field(const bool is_death, stk::mesh::MetaData & meta, CDFEM_Support & cdfem_support)
  {
    AuxMetaData & aux_meta = AuxMetaData::get(meta);
    Phase_Support::set_one_levelset_per_phase(true);
    ThrowRequire(!is_death);
    for(unsigned i=0; i < NUM_LS; ++i)
    {
      const std::string isovar_name = "LS" + std::to_string(i);
      FieldRef ls_isovar = aux_meta.declare_field(isovar_name, FieldType::REAL, stk::topology::NODE_RANK, 1u);
      ls_isovars.push_back(ls_isovar);
      cdfem_support.add_interpolation_field(ls_isovar);

      LevelSet * ls_ptr = nullptr;
      auto ls_field = std::make_shared<LS_Field>(isovar_name, LevelSet_Identifier(i), ls_isovar, 0., ls_ptr);
      ls_fields.push_back(ls_field);
      cdfem_support.add_ls_field(*ls_field, nullptr);
    }
  }

  void register_ls_on_blocks(const stk::mesh::PartVector & blocks, stk::mesh::MetaData & meta,
      Block_Surface_Connectivity & block_surface_info, const bool register_fields)
  {
    AuxMetaData & aux_meta = AuxMetaData::get(meta);
    PhaseVec named_phases;
    for(unsigned ls=0; ls < NUM_LS; ++ls)
    {
      PhaseTag tag;
      tag.add(LevelSet_Identifier(ls), -1);
      named_phases.push_back(NamedPhase("P" + std::to_string(ls), tag));
    }

    Phase_Support & phase_support = Phase_Support::get(meta);
    phase_support.set_input_block_surface_connectivity(block_surface_info);
    LevelSet * ls_ptr = nullptr;
    phase_support.register_blocks_for_level_set(ls_ptr, blocks);
    std::vector<std::tuple<stk::mesh::PartVector, 
      std::shared_ptr<Interface_Name_Generator>, PhaseVec>> ls_sets;
    auto interface_name_gen = std::shared_ptr<Interface_Name_Generator>(new LS_Name_Generator());
    ls_sets.push_back(std::make_tuple(blocks,interface_name_gen,named_phases));
    phase_support.decompose_blocks(ls_sets);

    for(auto && b : blocks)
    {
      for( unsigned ls=0; ls < NUM_LS; ++ls )
      {
        auto conformal_part = phase_support.find_conformal_io_part(*b, named_phases[ls].tag());
        ThrowRequire(conformal_part != nullptr);
        if(register_fields)
        {
          // Need to register every LS on every conformal part
          for( unsigned ls2=0; ls2 < NUM_LS; ++ls2 )
          {
            aux_meta.register_field(ls_isovars[ls2].name(), FieldType::REAL, stk::topology::NODE_RANK, 1u, 1u, *conformal_part);
          }
          aux_meta.register_field(ls_isovars[ls].name(), FieldType::REAL, stk::topology::NODE_RANK, 1u, 1u, *b);
        }
      }
    }
  }

  std::vector<FieldRef> ls_isovars;
  std::vector<std::shared_ptr<LS_Field> > ls_fields;
};

template <unsigned NUM_LS>
struct LSPerInterfacePolicy
{
  void setup_ls_field(const bool is_death, stk::mesh::MetaData & meta, CDFEM_Support & cdfem_support)
  {
    AuxMetaData & aux_meta = AuxMetaData::get(meta);
    Phase_Support::set_one_levelset_per_phase(false);
    ThrowRequire(!is_death);
    for(unsigned i=0; i < NUM_LS; ++i)
    {
      const std::string isovar_name = "LS" + std::to_string(i);
      FieldRef ls_isovar = aux_meta.declare_field(isovar_name, FieldType::REAL, stk::topology::NODE_RANK, 1u);
      ls_isovars.push_back(ls_isovar);
      cdfem_support.add_interpolation_field(ls_isovar);

      LevelSet * ls_ptr = nullptr;
      auto ls_field = std::make_shared<LS_Field>(isovar_name, LevelSet_Identifier(i), ls_isovar, 0., ls_ptr);
      ls_fields.push_back(ls_field);
      cdfem_support.add_ls_field(*ls_field, nullptr);
    }
  }

  void register_ls_on_blocks(const stk::mesh::PartVector & blocks, stk::mesh::MetaData & meta,
      Block_Surface_Connectivity & block_surface_info, const bool register_fields)
  {
    AuxMetaData & aux_meta = AuxMetaData::get(meta);
    PhaseVec named_phases;
    const unsigned numPhases = 1<<NUM_LS;
    for (unsigned phase=0; phase<numPhases; ++phase)
    {
      std::string phaseName = "P";
      PhaseTag tag;
      for(unsigned ls=0; ls < NUM_LS; ++ls)
      {
        const bool lsIsNeg = (phase>>ls)%2 == 0;
        const int lsSign = lsIsNeg ? -1 : 1;
        tag.add(LevelSet_Identifier(ls), lsSign);
        phaseName += (lsIsNeg ? "-" : "+");
      }
      named_phases.push_back(NamedPhase(phaseName, tag));
    }

    Phase_Support & phase_support = Phase_Support::get(meta);
    phase_support.set_input_block_surface_connectivity(block_surface_info);
    LevelSet * ls_ptr = nullptr;
    phase_support.register_blocks_for_level_set(ls_ptr, blocks);
    std::vector<std::tuple<stk::mesh::PartVector, 
      std::shared_ptr<Interface_Name_Generator>, PhaseVec>> ls_sets;
    auto interface_name_gen = std::shared_ptr<Interface_Name_Generator>(new LS_Name_Generator());
    ls_sets.push_back(std::make_tuple(blocks,interface_name_gen,named_phases));
    phase_support.decompose_blocks(ls_sets);

    if(register_fields)
    {
      for(auto && b : blocks)
      {
        for( unsigned ls=0; ls < NUM_LS; ++ls )
        {
          aux_meta.register_field(ls_isovars[ls].name(), FieldType::REAL, stk::topology::NODE_RANK, 1u, 1u, *b);

          for (unsigned phase=0; phase<numPhases; ++phase)
          {
            auto conformal_part = phase_support.find_conformal_io_part(*b, named_phases[phase].tag());
            ThrowRequire(conformal_part != nullptr);
            aux_meta.register_field(ls_isovars[ls].name(), FieldType::REAL, stk::topology::NODE_RANK, 1u, 1u, *conformal_part);
          }
        }
      }
    }
  }

  std::vector<FieldRef> ls_isovars;
  std::vector<std::shared_ptr<LS_Field> > ls_fields;
};

template <class MESH_FIXTURE, class LS_FIELD_POLICY>
class CompleteDecompositionFixture : public ::testing::Test
{
public:
  CompleteDecompositionFixture()
  : fixture(), cdfemSupport(CDFEM_Support::get(fixture.meta_data()))
  {
    AuxMetaData & aux_meta = AuxMetaData::get(fixture.meta_data());
    auto & vec_type = fixture.meta_data().spatial_dimension() == 2 ? FieldType::VECTOR_2D : FieldType::VECTOR_3D;
    coord_field = aux_meta.register_field("coordinates", vec_type, stk::topology::NODE_RANK, 1u, 1u, fixture.meta_data().universal_part());
    cdfemSupport.set_coords_field(coord_field);
    cdfemSupport.add_interpolation_field(coord_field);
    cdfemSupport.register_parent_node_ids_field();

    cdfemSupport.set_prolongation_model(INTERPOLATION);
  }

  void setup_ls_field(const bool is_death = false)
  {
    ls_policy.setup_ls_field(is_death, fixture.meta_data(), cdfemSupport);
  }

  void register_ls_on_blocks(const stk::mesh::PartVector & blocks, const bool register_fields = true)
  {
    ls_policy.register_ls_on_blocks(blocks, fixture.meta_data(), block_surface_info, register_fields);
  }

  stk::mesh::Part & declare_input_block(const std::string & name, const stk::topology topo)
  {
    auto & block_part = fixture.meta_data().declare_part_with_topology(name, topo);
    stk::io::put_io_part_attribute(block_part);
    return block_part;
  }

  stk::mesh::Part & declare_input_surface(const std::string & name, const stk::topology topo, const std::set<stk::mesh::PartOrdinal> & touching_blocks)
  {
    auto & surface_part = fixture.meta_data().declare_part_with_topology(name, topo);
    stk::io::put_io_part_attribute(surface_part);

    block_surface_info.add_surface(surface_part.mesh_meta_data_ordinal(), touching_blocks);
    return surface_part;
  }

  stk::mesh::PartVector declare_input_surfaces_touching_block(const unsigned numSurfaces, const stk::mesh::Part & touchingBlock)
  {
    const stk::topology topo = touchingBlock.topology().side_topology();
    stk::mesh::PartVector surfaces;
    for (unsigned i=0; i<numSurfaces; ++i)
    {
      const std::string name = "InputSurface" + std::to_string(i+1);
      auto & surfacePart = fixture.meta_data().declare_part_with_topology(name, topo);
      stk::io::put_io_part_attribute(surfacePart);

      block_surface_info.add_surface(surfacePart.mesh_meta_data_ordinal(), {touchingBlock.mesh_meta_data_ordinal()});
      surfaces.push_back(&surfacePart);
    }

    return surfaces;
  }

  void decompose_mesh()
  {
    NodeToCapturedDomainsMap nodesToSnappedDomains;
    LevelSetInterfaceGeometry interfaceGeometry(krino_mesh->get_active_part(), cdfemSupport, Phase_Support::get(fixture.meta_data()));
    if (cdfemSupport.get_cdfem_edge_degeneracy_handling() == SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE)
      nodesToSnappedDomains = snap_as_much_as_possible_while_maintaining_quality(krino_mesh->stk_bulk(), krino_mesh->get_active_part(), cdfemSupport.get_interpolation_fields(), interfaceGeometry, cdfemSupport.get_global_ids_are_parallel_consistent());
    interfaceGeometry.prepare_to_process_elements(krino_mesh->stk_bulk(), nodesToSnappedDomains);

    if(!krino_mesh->my_old_mesh)
    {
      krino_mesh->my_old_mesh = std::make_shared<CDMesh>(fixture.bulk_data(), std::shared_ptr<CDMesh>());
      krino_mesh->my_old_mesh->generate_nonconformal_elements();
    }

    krino_mesh->generate_nonconformal_elements();
    if (cdfemSupport.get_cdfem_edge_degeneracy_handling() == SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE)
      krino_mesh->snap_nearby_intersections_to_nodes(interfaceGeometry, nodesToSnappedDomains);
    krino_mesh->set_phase_of_uncut_elements(interfaceGeometry);
    krino_mesh->triangulate(interfaceGeometry);
    krino_mesh->my_old_mesh->stash_field_data(-1, *krino_mesh);


    krino_mesh->decompose();
    krino_mesh->modify_mesh();
    krino_mesh->prolongation();

    if(krinolog.shouldPrint(LOG_DEBUG))
    {
      krino_mesh->debug_output();
    }
  }

  void debug_output()
  {
    krino_mesh->debug_output();
  }

  void commit()
  {
    krino_mesh = std::make_shared<CDMesh>(fixture.bulk_data(), std::shared_ptr<CDMesh>());
  }

  stk::mesh::Entity create_element(stk::mesh::PartVector & elem_parts, stk::mesh::EntityId elem_id,
      std::vector<stk::mesh::EntityId> elem_nodes)
  {
    auto elem = stk::mesh::declare_element( fixture.bulk_data(), elem_parts, elem_id, elem_nodes );
    {
      const stk::mesh::Entity * const nodes = fixture.bulk_data().begin_nodes(elem);
      for(unsigned i=0; i < elem_nodes.size(); ++i)
      {
        EXPECT_EQ(elem_nodes[i], fixture.bulk_data().identifier(nodes[i]));
        if (!fixture.bulk_data().bucket(nodes[i]).member(cdfemSupport.get_active_part()))
          fixture.bulk_data().change_entity_parts(nodes[i], stk::mesh::ConstPartVector{&cdfemSupport.get_active_part()}, {});
      }
    }
    return elem;
  }

  void run_rebalance_with(const std::string& decomp_method)
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

    if (parallel_size != 2) return;

    cdfemSupport.set_cdfem_edge_tol(0.1);
    cdfemSupport.set_simplex_generation_method(CUT_QUADS_BY_GLOBAL_IDENTIFIER);

    setup_ls_field();

    const stk::topology tet4 = stk::topology::TETRAHEDRON_4;
    auto & block1_part = declare_input_block("block_1", tet4);
    auto & block2_part = declare_input_block("block_2", tet4);
    auto & surface_part = declare_input_surface("surface_1", tet4.side_topology(), {block1_part.mesh_meta_data_ordinal(), block2_part.mesh_meta_data_ordinal()});

    register_ls_on_blocks({&block1_part, &block2_part});

    FieldRef elem_weight_field =
        aux_meta.register_field("element_weight", FieldType::REAL, stk::topology::ELEMENT_RANK,
            1u, 1, fixture.meta_data().universal_part());

    commit();

    const bool build_all_on_P0 = true;
    build_two_tet4_mesh_np2(*this, block1_part, block2_part, surface_part, parallel_rank,
        parallel_size, true, build_all_on_P0);

    stk::mesh::create_exposed_block_boundary_sides(mesh, meta.universal_part(), {&aux_meta.exposed_boundary_part()});

    if(parallel_rank == 0)
    {
      auto & node_buckets = mesh.buckets(stk::topology::NODE_RANK);
      for(auto && b_ptr : node_buckets)
      {
        for(auto && node : *b_ptr)
        {
          const double * coords = field_data<double>(coord_field, node);
          double * ls_data = field_data<double>(ls_policy.ls_isovar, node);
          if(ls_data) *ls_data = coords[1]-0.5;
        }
      }
    }
    stk::mesh::communicate_field_data(mesh, {&ls_policy.ls_isovar.field(), &coord_field.field()});

    ASSERT_NO_THROW(decompose_mesh());

    if(parallel_rank == 1)
    {
      EXPECT_EQ(0u, stk::mesh::get_num_entities(mesh));
    }

    const stk::mesh::Selector parent_selector = krino_mesh->get_parent_part();
    if(parallel_rank == 0)
    {
      auto & elem_buckets = mesh.buckets(stk::topology::ELEMENT_RANK);
      for(auto && b_ptr : elem_buckets)
      {
        const bool is_parent = parent_selector(*b_ptr);
        for(auto && elem : *b_ptr)
        {
          double * weight = field_data<double>(elem_weight_field, elem);
          *weight = is_parent ? 1. : 0.;
        }
      }
    }

    rebalance_utils::rebalance_mesh(mesh,
        krino_mesh.get(),
        elem_weight_field.name(),
        coord_field.name(),
        {fixture.meta_data().universal_part()},
        decomp_method);

    // Both procs should now own 1 parent element and 4 children
    EXPECT_EQ(1u, stk::mesh::count_selected_entities(
        fixture.meta_data().locally_owned_part() & parent_selector,
        mesh.buckets(stk::topology::ELEMENT_RANK)));
    EXPECT_EQ(4u, stk::mesh::count_selected_entities(
        fixture.meta_data().locally_owned_part() & krino_mesh->get_child_part(),
        mesh.buckets(stk::topology::ELEMENT_RANK)));

    krino_mesh = std::make_shared<CDMesh>(mesh, krino_mesh);

    auto & node_buckets = mesh.buckets(stk::topology::NODE_RANK);
    for(auto && b_ptr : node_buckets)
    {
      for(auto && node : *b_ptr)
      {
        const double * coords = field_data<double>(coord_field, node);
        double * ls_data = field_data<double>(ls_policy.ls_isovar, node);
        if(ls_data) *ls_data = coords[2]-0.5;
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
  std::shared_ptr<CDMesh> krino_mesh;
  Block_Surface_Connectivity block_surface_info;
  LS_FIELD_POLICY ls_policy;
  LogRedirecter log;
};

namespace {
template <class MESH_FIXTURE, class LS_FIELD_POLICY>
void build_two_tri3_mesh_on_one_or_two_procs(CompleteDecompositionFixture<MESH_FIXTURE, LS_FIELD_POLICY> & fixture,
    stk::mesh::Part & block_part, const int parallel_rank)
{
  stk::mesh::BulkData & mesh = fixture.fixture.bulk_data();
  stk::mesh::MetaData & meta = fixture.fixture.meta_data();
  const int parallel_size = mesh.parallel_size();
  mesh.modification_begin();
  {
    /*
     *   4---3
     *   |\ 2|
     *   | \ |
     *   |1 \|
     *   1---2
     */
    stk::mesh::PartVector elem_parts;
    elem_parts.push_back(&block_part);
    elem_parts.push_back(&fixture.cdfemSupport.get_active_part());

    std::vector<stk::mesh::EntityId> elem1_nodes = {1, 2, 4};
    std::vector<stk::mesh::EntityId> elem2_nodes = {2, 3, 4};

    if(parallel_rank == 0) fixture.create_element(elem_parts, 1, elem1_nodes);
    if(parallel_rank == 1 || parallel_size == 1) fixture.create_element(elem_parts, 2, elem2_nodes);

    if(parallel_size > 1)
    {
      const int opp_rank = parallel_rank == 0 ? 1 : 0;
      const auto node2 = mesh.get_entity(stk::topology::NODE_RANK, 2);
      const auto node4 = mesh.get_entity(stk::topology::NODE_RANK, 4);
      mesh.add_node_sharing(node2, opp_rank);
      mesh.add_node_sharing(node4, opp_rank);
    }
  }
  mesh.modification_end();

  EXPECT_EQ(2u, stk::mesh::count_selected_entities(meta.universal_part(), mesh.buckets(stk::topology::ELEMENT_RANK)));

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

  node1_coords[0] = 0.;
  node1_coords[1] = 0.;

  node2_coords[0] = 1.;
  node2_coords[1] = 0.;

  node3_coords[0] = 1.;
  node3_coords[1] = 1.;

  node4_coords[0] = 0.;
  node4_coords[1] = 1.;
}

template <class MESH_FIXTURE, class LS_FIELD_POLICY>
void build_two_tri3_mesh_on_one_or_two_procs_with_sides(CompleteDecompositionFixture<MESH_FIXTURE, LS_FIELD_POLICY> & fixture,
    stk::mesh::Part & blockPart, const stk::mesh::PartVector & sideParts, const int parallel_rank)
{
  build_two_tri3_mesh_on_one_or_two_procs(fixture, blockPart, parallel_rank);

  ThrowRequire(sideParts.size() == 4);
  stk::mesh::BulkData & mesh = fixture.fixture.bulk_data();
  const int parallel_size = mesh.parallel_size();

  mesh.modification_begin();
  if(parallel_rank == 0)
  {
    const auto element1 = mesh.get_entity(stk::topology::ELEMENT_RANK,1);
    mesh.declare_element_side(element1, 0, stk::mesh::PartVector{sideParts[0]});
    mesh.declare_element_side(element1, 2, stk::mesh::PartVector{sideParts[3]});
  }
  if(parallel_rank == 1 || parallel_size == 1)
  {
    const auto element2 = mesh.get_entity(stk::topology::ELEMENT_RANK,2);
    mesh.declare_element_side(element2, 0, stk::mesh::PartVector{sideParts[1]});
    mesh.declare_element_side(element2, 1, stk::mesh::PartVector{sideParts[2]});
  }
  mesh.modification_end();
}

void check_two_tri3_snapped_mesh_np2(CompleteDecompositionFixture<SimpleStkFixture2d, SingleLSPolicy> & fixture, const int parallel_rank)
{
  stk::mesh::BulkData & mesh = fixture.fixture.bulk_data();
  stk::mesh::MetaData & meta = fixture.fixture.meta_data();
  AuxMetaData & aux_meta = AuxMetaData::get(meta);
  std::vector<stk::mesh::Entity> entities;

  // Should be no new nodes because of snapped interface, 3 nodes owned by P0, 1 by P1
  mesh.get_entities(stk::topology::NODE_RANK, meta.universal_part(), entities);
  EXPECT_EQ(4u, entities.size());
  mesh.get_entities(stk::topology::NODE_RANK, aux_meta.active_part(), entities);
  EXPECT_EQ(4u, entities.size());
  mesh.get_entities(stk::topology::NODE_RANK, aux_meta.active_part() & meta.locally_owned_part(), entities);
  if(parallel_rank == 0)
  {
    EXPECT_EQ(3u, entities.size());
  }
  else
  {
    EXPECT_EQ(1u, entities.size());
  }

  // Should be 1 interface edge
  mesh.get_entities(stk::topology::EDGE_RANK, aux_meta.active_part() &
      aux_meta.get_part("surface_block_1_A_B") & aux_meta.get_part("surface_block_1_B_A"), entities);
  EXPECT_EQ(1u, entities.size());

  // Should be 2 coincident subelements, no parents
  mesh.get_entities(stk::topology::ELEMENT_RANK, meta.universal_part(), entities);
  EXPECT_EQ(2u, entities.size());
  mesh.get_entities(stk::topology::ELEMENT_RANK, aux_meta.active_part(), entities);
  EXPECT_EQ(2u, entities.size());
  mesh.get_entities(stk::topology::ELEMENT_RANK, !aux_meta.active_part(), entities);
  EXPECT_EQ(0u, entities.size());
  mesh.get_entities(stk::topology::ELEMENT_RANK, aux_meta.get_part("block_1_A"), entities);
  EXPECT_EQ(1u, entities.size());
  mesh.get_entities(stk::topology::ELEMENT_RANK, aux_meta.get_part("block_1_B"), entities);
  EXPECT_EQ(1u, entities.size());
}

}

typedef CompleteDecompositionFixture<SimpleStkFixture2d, SingleLSPolicy> CDMeshTests2D;
TEST_F(CDMeshTests2D, IsovariableNotDefinedOnDecomposedBlock)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_size = stk::parallel_machine_size(pm);

  if (parallel_size > 1) return;

  setup_ls_field();

  auto & block_part = declare_input_block("block_1", stk::topology::TRIANGLE_3_2D);
  register_ls_on_blocks({&block_part}, false);

  commit();

  EXPECT_ANY_THROW(krino_mesh->check_isovariable_field_existence_on_decomposed_blocks(true));
}

TEST_F(CDMeshTests2D, IsovariableNotDefinedOnBlock1)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_size = stk::parallel_machine_size(pm);

  if (parallel_size > 1) return;

  setup_ls_field();

  auto & block_part = declare_input_block("block_1", stk::topology::TRIANGLE_3_2D);
  register_ls_on_blocks({&block_part}, false);

  auto & meta = fixture.meta_data();
  auto & aux_meta = AuxMetaData::get(meta);

  auto * A_part = meta.get_part("block_1_A");
  ASSERT_TRUE(A_part != nullptr);
  auto * B_part = meta.get_part("block_1_B");
  ASSERT_TRUE(B_part != nullptr);

  // Catch the case where the field exists on both conformal parts, but not
  // the initial un-decomposed part so we can't do the initial decomposition.
  aux_meta.register_field("LS", FieldType::REAL, stk::topology::NODE_RANK,
      1u, 1u, *A_part);
  aux_meta.register_field("LS", FieldType::REAL, stk::topology::NODE_RANK,
      1u, 1u, *B_part);

  commit();

  EXPECT_ANY_THROW(krino_mesh->check_isovariable_field_existence_on_decomposed_blocks(true));
}

TEST_F(CDMeshTests2D, IsovariableOnlyOnBlock1SteadyState)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_size = stk::parallel_machine_size(pm);

  if (parallel_size > 1) return;

  setup_ls_field();

  auto & block_part = declare_input_block("block_1", stk::topology::TRIANGLE_3_2D);
  register_ls_on_blocks({&block_part}, false);

  auto & meta = fixture.meta_data();
  auto & aux_meta = AuxMetaData::get(meta);

  aux_meta.register_field("LS", FieldType::REAL, stk::topology::NODE_RANK,
      1u, 1u, block_part);

  commit();

  EXPECT_NO_THROW(krino_mesh->check_isovariable_field_existence_on_decomposed_blocks(false));
}

TEST_F(CDMeshTests2D, DeathIsovariableNotDefinedOnDecomposedBlock)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_size = stk::parallel_machine_size(pm);

  if (parallel_size > 1) return;

  const bool is_death = true;
  setup_ls_field(is_death);

  auto & block_part = declare_input_block("block_1", stk::topology::TRIANGLE_3_2D);
  register_ls_on_blocks({&block_part}, false);

  commit();

  EXPECT_ANY_THROW(krino_mesh->check_isovariable_field_existence_on_decomposed_blocks(true));
}

TEST_F(CDMeshTests2D, DeathIsovariableNotDefinedOnDeadBlock)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_size = stk::parallel_machine_size(pm);

  if (parallel_size > 1) return;

  const bool is_death = true;
  setup_ls_field(is_death);

  auto & block_part = declare_input_block("block_1", stk::topology::TRIANGLE_3_2D);
  register_ls_on_blocks({&block_part}, false);

  auto & meta = fixture.meta_data();
  auto & aux_meta = AuxMetaData::get(meta);
  auto * alive_part = meta.get_part("block_1_A");
  ASSERT_TRUE(alive_part != nullptr);

  aux_meta.register_field("LS", FieldType::REAL, stk::topology::NODE_RANK,
      1u, 1u, block_part);
  aux_meta.register_field("LS", FieldType::REAL, stk::topology::NODE_RANK,
      1u, 1u, *alive_part);

  commit();

  EXPECT_NO_THROW(krino_mesh->check_isovariable_field_existence_on_decomposed_blocks(true));
}

typedef CompleteDecompositionFixture<SimpleStkFixture2d, SingleLSPolicy> CDMeshTests2D;
TEST_F(CDMeshTests2D, Single_Tri3_Decomposition)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_size = stk::parallel_machine_size(pm);

  if (parallel_size > 1) return;

  setup_ls_field();

  auto & block_part = declare_input_block("block_1", stk::topology::TRIANGLE_3_2D);
  register_ls_on_blocks({&block_part});

  commit();

  stk::mesh::BulkData & mesh = fixture.bulk_data();
  stk::mesh::MetaData & meta = fixture.meta_data();
  AuxMetaData & aux_meta = AuxMetaData::get(meta);

  mesh.modification_begin();
  {
    stk::mesh::PartVector elem_parts;
    elem_parts.push_back(&block_part);
    elem_parts.push_back(&aux_meta.active_part());
    std::vector<stk::mesh::EntityId> node_ids = {1, 2, 3};

    create_element(elem_parts, 1, node_ids);
  }
  mesh.modification_end();

  stk::mesh::create_exposed_block_boundary_sides(mesh, meta.universal_part(), {&aux_meta.exposed_boundary_part()});

  const auto node1 = mesh.get_entity(stk::topology::NODE_RANK, 1);
  const auto node2 = mesh.get_entity(stk::topology::NODE_RANK, 2);
  const auto node3 = mesh.get_entity(stk::topology::NODE_RANK, 3);

  double * node1_coords = field_data<double>(coord_field, node1);
  double * node2_coords = field_data<double>(coord_field, node2);
  double * node3_coords = field_data<double>(coord_field, node3);

  node1_coords[0] = 0.;
  node1_coords[1] = 0.;

  node2_coords[0] = 1.;
  node2_coords[1] = 0.;

  node3_coords[0] = 0.;
  node3_coords[1] = 1.;

  *field_data<double>(ls_policy.ls_isovar, node1) = -1;
  *field_data<double>(ls_policy.ls_isovar, node2) = 1;
  *field_data<double>(ls_policy.ls_isovar, node3) = -1;

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

  std::vector<stk::mesh::Entity> entities;

  // Should have added 2 nodes at the cutting locations
  mesh.get_entities(stk::topology::NODE_RANK, meta.universal_part(), entities);
  EXPECT_EQ(5u, entities.size());
  mesh.get_entities(stk::topology::NODE_RANK, aux_meta.active_part(), entities);
  EXPECT_EQ(5u, entities.size());

  // Should be 1 interface edge
  mesh.get_entities(stk::topology::EDGE_RANK, aux_meta.active_part() &
      aux_meta.get_part("surface_block_1_A_B") & aux_meta.get_part("surface_block_1_B_A"), entities);
  EXPECT_EQ(1u, entities.size());

  // Should be 3 conformal elements plus the parent element
  mesh.get_entities(stk::topology::ELEMENT_RANK, meta.universal_part(), entities);
  EXPECT_EQ(4u, entities.size());
  mesh.get_entities(stk::topology::ELEMENT_RANK, aux_meta.active_part(), entities);
  EXPECT_EQ(3u, entities.size());
  mesh.get_entities(stk::topology::ELEMENT_RANK, !aux_meta.active_part(), entities);
  EXPECT_EQ(1u, entities.size());
  mesh.get_entities(stk::topology::ELEMENT_RANK, aux_meta.get_part("block_1_A"), entities);
  EXPECT_EQ(1u, entities.size());
  mesh.get_entities(stk::topology::ELEMENT_RANK, aux_meta.get_part("block_1_B"), entities);
  EXPECT_EQ(2u, entities.size());
}

TEST_F(CDMeshTests2D, Two_Tri3_Snapped_Interface_Decomposition_NP2)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_size = stk::parallel_machine_size(pm);
  const int parallel_rank = stk::parallel_machine_rank(pm);

  if (parallel_size != 2) return;

  stk::mesh::BulkData & mesh = fixture.bulk_data();
  stk::mesh::MetaData & meta = fixture.meta_data();
  AuxMetaData & aux_meta = AuxMetaData::get(meta);

  setup_ls_field();

  auto & block_part = declare_input_block("block_1", stk::topology::TRIANGLE_3_2D);
  register_ls_on_blocks({&block_part});

  commit();

  build_two_tri3_mesh_on_one_or_two_procs(*this, block_part, parallel_rank);

  stk::mesh::create_exposed_block_boundary_sides(mesh, meta.universal_part(), {&aux_meta.exposed_boundary_part()});

  {
    const auto node1 = mesh.get_entity(stk::topology::NODE_RANK, 1);
    const auto node2 = mesh.get_entity(stk::topology::NODE_RANK, 2);
    const auto node3 = mesh.get_entity(stk::topology::NODE_RANK, 3);
    const auto node4 = mesh.get_entity(stk::topology::NODE_RANK, 4);

    *field_data<double>(ls_policy.ls_isovar, node1) = -1;
    *field_data<double>(ls_policy.ls_isovar, node2) = 0;
    *field_data<double>(ls_policy.ls_isovar, node3) = 1;
    *field_data<double>(ls_policy.ls_isovar, node4) = 0;

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

    check_two_tri3_snapped_mesh_np2(*this, parallel_rank);

    std::vector<stk::mesh::Entity> entities;
    mesh.get_entities(stk::topology::ELEMENT_RANK, aux_meta.get_part("block_1_A") & meta.locally_owned_part(), entities);
    if(parallel_rank == 0)
    {
      EXPECT_TRUE(entities.empty());
    }
    else
    {
      EXPECT_EQ(1u, entities.size());
    }
    mesh.get_entities(stk::topology::ELEMENT_RANK, aux_meta.get_part("block_1_B") & meta.locally_owned_part(), entities);
    if(parallel_rank == 0)
    {
      EXPECT_EQ(1u, entities.size());
    }
    else
    {
      EXPECT_TRUE(entities.empty());
    }

    krino_mesh = std::make_shared<CDMesh>(mesh, krino_mesh);
  }

  // Swap the A and B elements
  {
    const auto node1 = mesh.get_entity(stk::topology::NODE_RANK, 1);
    const auto node2 = mesh.get_entity(stk::topology::NODE_RANK, 2);
    const auto node3 = mesh.get_entity(stk::topology::NODE_RANK, 3);
    const auto node4 = mesh.get_entity(stk::topology::NODE_RANK, 4);

    *field_data<double>(ls_policy.ls_isovar, node1) = 1;
    *field_data<double>(ls_policy.ls_isovar, node2) = 0;
    *field_data<double>(ls_policy.ls_isovar, node3) = -1;
    *field_data<double>(ls_policy.ls_isovar, node4) = 0;

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

    check_two_tri3_snapped_mesh_np2(*this, parallel_rank);

    std::vector<stk::mesh::Entity> entities;
    mesh.get_entities(stk::topology::ELEMENT_RANK, aux_meta.get_part("block_1_A") & meta.locally_owned_part(), entities);
    if(parallel_rank == 0)
    {
      EXPECT_EQ(1u, entities.size());
    }
    else
    {
      EXPECT_TRUE(entities.empty());
    }
    mesh.get_entities(stk::topology::ELEMENT_RANK, aux_meta.get_part("block_1_B") & meta.locally_owned_part(), entities);
    if(parallel_rank == 0)
    {
      EXPECT_TRUE(entities.empty());
    }
    else
    {
      EXPECT_EQ(1u, entities.size());
    }
  }
}

TEST_F(CDMeshTests2D, Two_Tri3_Death_Snapped_Interface_Decomposition_NP2)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_size = stk::parallel_machine_size(pm);
  const int parallel_rank = stk::parallel_machine_rank(pm);

  if (parallel_size != 2) return;

  stk::mesh::BulkData & mesh = fixture.bulk_data();
  stk::mesh::MetaData & meta = fixture.meta_data();
  AuxMetaData & aux_meta = AuxMetaData::get(meta);

  setup_ls_field(true);

  auto & block_part = declare_input_block("block_1", stk::topology::TRIANGLE_3_2D);
  register_ls_on_blocks({&block_part});

  commit();

  build_two_tri3_mesh_on_one_or_two_procs(*this, block_part, parallel_rank);

  stk::mesh::create_exposed_block_boundary_sides(mesh, meta.universal_part(), {&aux_meta.exposed_boundary_part()});

  {
    const auto node1 = mesh.get_entity(stk::topology::NODE_RANK, 1);
    const auto node2 = mesh.get_entity(stk::topology::NODE_RANK, 2);
    const auto node3 = mesh.get_entity(stk::topology::NODE_RANK, 3);
    const auto node4 = mesh.get_entity(stk::topology::NODE_RANK, 4);

    *field_data<double>(ls_policy.ls_isovar, node1) = -1;
    *field_data<double>(ls_policy.ls_isovar, node2) = 0;
    *field_data<double>(ls_policy.ls_isovar, node3) = 1;
    *field_data<double>(ls_policy.ls_isovar, node4) = 0;

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

    check_two_tri3_snapped_mesh_np2(*this, parallel_rank);

    std::vector<stk::mesh::Entity> entities;
    mesh.get_entities(stk::topology::ELEMENT_RANK, aux_meta.get_part("block_1_A") & meta.locally_owned_part(), entities);
    if(parallel_rank == 0)
    {
      EXPECT_TRUE(entities.empty());
    }
    else
    {
      EXPECT_EQ(1u, entities.size());
    }
    mesh.get_entities(stk::topology::ELEMENT_RANK, aux_meta.get_part("block_1_B") & meta.locally_owned_part(), entities);
    if(parallel_rank == 0)
    {
      EXPECT_EQ(1u, entities.size());
    }
    else
    {
      EXPECT_TRUE(entities.empty());
    }
  }
}

TEST_F(CDMeshTests2D, Two_Tri3_Death_Snapped_Interface_Decomposition_NP2_Opposite_Ownership)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_size = stk::parallel_machine_size(pm);
  const int parallel_rank = stk::parallel_machine_rank(pm);

  if (parallel_size != 2) return;

  stk::mesh::BulkData & mesh = fixture.bulk_data();
  stk::mesh::MetaData & meta = fixture.meta_data();
  AuxMetaData & aux_meta = AuxMetaData::get(meta);

  setup_ls_field(true);

  auto & block_part = declare_input_block("block_1", stk::topology::TRIANGLE_3_2D);
  register_ls_on_blocks({&block_part});

  commit();

  build_two_tri3_mesh_on_one_or_two_procs(*this, block_part, parallel_rank);

  stk::mesh::create_exposed_block_boundary_sides(mesh, meta.universal_part(), {&aux_meta.exposed_boundary_part()});

  {
    const auto node1 = mesh.get_entity(stk::topology::NODE_RANK, 1);
    const auto node2 = mesh.get_entity(stk::topology::NODE_RANK, 2);
    const auto node3 = mesh.get_entity(stk::topology::NODE_RANK, 3);
    const auto node4 = mesh.get_entity(stk::topology::NODE_RANK, 4);

    *field_data<double>(ls_policy.ls_isovar, node1) = 1;
    *field_data<double>(ls_policy.ls_isovar, node2) = 0;
    *field_data<double>(ls_policy.ls_isovar, node3) = -1;
    *field_data<double>(ls_policy.ls_isovar, node4) = 0;

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

    check_two_tri3_snapped_mesh_np2(*this, parallel_rank);

    std::vector<stk::mesh::Entity> entities;
    mesh.get_entities(stk::topology::ELEMENT_RANK, aux_meta.get_part("block_1_A") & meta.locally_owned_part(), entities);
    if(parallel_rank == 0)
    {
      EXPECT_EQ(1u, entities.size());
    }
    else
    {
      EXPECT_TRUE(entities.empty());
    }
    mesh.get_entities(stk::topology::ELEMENT_RANK, aux_meta.get_part("block_1_B") & meta.locally_owned_part(), entities);
    if(parallel_rank == 0)
    {
      EXPECT_TRUE(entities.empty());
    }
    else
    {
      EXPECT_EQ(1u, entities.size());
    }
  }
}

TEST_F(CDMeshTests2D, Two_Tri3_Periodic)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_size = stk::parallel_machine_size(pm);
  const int parallel_rank = stk::parallel_machine_rank(pm);

  if (parallel_size > 2) return;

  cdfemSupport.set_cdfem_edge_tol(0.1);

  stk::mesh::BulkData & mesh = fixture.bulk_data();
  AuxMetaData & aux_meta = AuxMetaData::get(fixture.meta_data());

  setup_ls_field();

  auto & block_part = declare_input_block("block_1", stk::topology::TRIANGLE_3_2D);
  register_ls_on_blocks({&block_part});

  commit();

  build_two_tri3_mesh_on_one_or_two_procs(*this, block_part, parallel_rank);

  stk::mesh::create_exposed_block_boundary_sides(mesh, fixture.meta_data().universal_part(), {&aux_meta.exposed_boundary_part()});

  {
    const auto node1 = mesh.get_entity(stk::topology::NODE_RANK, 1);
    const auto node2 = mesh.get_entity(stk::topology::NODE_RANK, 2);
    const auto node3 = mesh.get_entity(stk::topology::NODE_RANK, 3);
    const auto node4 = mesh.get_entity(stk::topology::NODE_RANK, 4);

    *field_data<double>(ls_policy.ls_isovar, node1) = -0.99;
    *field_data<double>(ls_policy.ls_isovar, node2) = -0.8;
    *field_data<double>(ls_policy.ls_isovar, node3) = 0.2;
    *field_data<double>(ls_policy.ls_isovar, node4) = 0.01;

    krino_mesh->add_periodic_node_pair(node3, node4);

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

    // Periodic constraint should cause interface to snap to both 3 and 4 so no elements
    // get cut.
    std::vector<stk::mesh::Entity> entities;
    mesh.get_entities(stk::topology::ELEMENT_RANK, aux_meta.get_part("block_1_A"), entities);
    EXPECT_EQ(0u, entities.size());
    mesh.get_entities(stk::topology::ELEMENT_RANK, aux_meta.get_part("block_1_B"), entities);
    EXPECT_EQ(2u, entities.size());
  }
}

TEST_F(CDMeshTests2D, Two_Tri3_PeriodicParallelNonSharedNode)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_size = stk::parallel_machine_size(pm);
  const int parallel_rank = stk::parallel_machine_rank(pm);

  if (parallel_size == 1 || (parallel_size % 2) != 0) return;

  cdfemSupport.set_cdfem_edge_tol(0.1);

  stk::mesh::BulkData & mesh = fixture.bulk_data();
  stk::mesh::MetaData & meta = fixture.meta_data();
  AuxMetaData & aux_meta = AuxMetaData::get(meta);

  setup_ls_field();

  auto & block_part = declare_input_block("block_1", stk::topology::TRIANGLE_3_2D);
  register_ls_on_blocks({&block_part});

  commit();

  const stk::mesh::EntityId starting_id = 1 + 3*parallel_rank;
  mesh.modification_begin();
  {
    stk::mesh::PartVector elem_parts;
    elem_parts.push_back(&block_part);
    elem_parts.push_back(&cdfemSupport.get_active_part());

    std::vector<stk::mesh::EntityId> elem_nodes = {starting_id, starting_id + 1, starting_id + 2};

    create_element(elem_parts, parallel_rank + 1, elem_nodes);
  }
  mesh.modification_end();

  stk::mesh::create_exposed_block_boundary_sides(mesh, meta.universal_part(), {&aux_meta.exposed_boundary_part()});

  EXPECT_EQ(1u, stk::mesh::count_selected_entities(meta.universal_part(), mesh.buckets(stk::topology::ELEMENT_RANK)));

  const auto node1 = mesh.get_entity(stk::topology::NODE_RANK, starting_id);
  const auto node2 = mesh.get_entity(stk::topology::NODE_RANK, starting_id + 1);
  const auto node3 = mesh.get_entity(stk::topology::NODE_RANK, starting_id + 2);

  ASSERT_TRUE(mesh.is_valid(node1));
  ASSERT_TRUE(mesh.is_valid(node2));
  ASSERT_TRUE(mesh.is_valid(node3));

  double * node1_coords = field_data<double>(coord_field, node1);
  double * node2_coords = field_data<double>(coord_field, node2);
  double * node3_coords = field_data<double>(coord_field, node3);

  const double x0 = 1.1 * parallel_rank;
  node1_coords[0] = x0;
  node1_coords[1] = 0.;

  node2_coords[0] = x0 + 1.;
  node2_coords[1] = 0.;

  node3_coords[0] = x0 + 1.;
  node3_coords[1] = 1.;

  // Going to add "periodic" constraint between node1 on P0 and P1, so ghost those nodes appropriately
  mesh.modification_begin();
  auto & ghosting = mesh.create_ghosting("test_ghosting");
  mesh.modification_end();
  stk::mesh::EntityProcVec add_send;
  const int mod = parallel_rank % 2;
  if(mod == 0)
  {
    add_send.push_back(std::make_pair(node1, parallel_rank + 1));
  }
  else
  {
    add_send.push_back(std::make_pair(node1, parallel_rank - 1));
  }
  mesh.modification_begin();
  mesh.change_ghosting(ghosting, add_send);
  mesh.modification_end();

  ASSERT_EQ(4u,
      stk::mesh::count_selected_entities(meta.universal_part(), mesh.buckets(stk::topology::NODE_RANK)));
  {
    // Procs with mod == 0 will setup isovars so that the 1-2 edge has a crossing within the snap tolerance
    // of node 1. Procs with mod == 1 will have the crossing outside the snap tolerance
    *field_data<double>(ls_policy.ls_isovar, node1) = (mod == 0) ? -0.01 : -1.;
    *field_data<double>(ls_policy.ls_isovar, node2) = 1.;
    *field_data<double>(ls_policy.ls_isovar, node3) = 1.;

    auto other_node1 = mesh.get_entity(stk::topology::NODE_RANK,
        (mod == 0) ? starting_id + 3 : starting_id - 3);
    krino_mesh->add_periodic_node_pair(node1, other_node1);

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

    // Periodic constraint should cause interface to snap to node1 so the element is entirely in the
    // A phase.
    std::vector<stk::mesh::Entity> entities;
    mesh.get_entities(stk::topology::ELEMENT_RANK, aux_meta.get_part("block_1_A"), entities);
    EXPECT_EQ(1u, entities.size());
    mesh.get_entities(stk::topology::ELEMENT_RANK, aux_meta.get_part("block_1_B"), entities);
    EXPECT_EQ(0u, entities.size());
  }
}

void set_level_set(const stk::mesh::BulkData & mesh, FieldRef lsField, const std::vector<stk::mesh::EntityId> & nodeIds, const std::vector<double> & nodeLS)
{
  ThrowRequire(nodeIds.size() == nodeLS.size());
  for (size_t i=0; i<nodeIds.size(); ++i)
  {
    const auto node = mesh.get_entity(stk::topology::NODE_RANK, nodeIds[i]);
    ThrowRequire(mesh.is_valid(node));
    *field_data<double>(lsField, node) = nodeLS[i];
  }
}

TEST_F(CDMeshTests2D, Two_Tri3_Check_Compatibility_When_Snapping)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_size = stk::parallel_machine_size(pm);
  const int parallel_rank = stk::parallel_machine_rank(pm);

  if (parallel_size > 2) return;

  stk::mesh::BulkData & mesh = fixture.bulk_data();
  AuxMetaData & aux_meta = AuxMetaData::get(fixture.meta_data());

  setup_ls_field();

  auto & block_part = declare_input_block("block_1", stk::topology::TRIANGLE_3_2D);
  register_ls_on_blocks({&block_part});

  commit();

  build_two_tri3_mesh_on_one_or_two_procs(*this, block_part, parallel_rank);

  stk::mesh::create_exposed_block_boundary_sides(mesh, fixture.meta_data().universal_part(), {&aux_meta.exposed_boundary_part()});

  set_level_set(mesh, ls_policy.ls_isovar, {1,2,3,4} , {-1., -1., 1.e20, 1.});

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

  std::vector<stk::mesh::Entity> entities;
  mesh.get_entities(stk::topology::ELEMENT_RANK, aux_meta.get_part("block_1_A"), entities);
  EXPECT_EQ(2u, entities.size());
  mesh.get_entities(stk::topology::ELEMENT_RANK, aux_meta.get_part("block_1_B"), entities);
  EXPECT_EQ(1u, entities.size());

  //fixture.write_results("Two_Tri3_Check_Compatibility_When_Snapping.e");
}

void set_level_sets(const stk::mesh::BulkData & mesh, std::vector<FieldRef> & lsFields, const std::vector<stk::mesh::EntityId> & nodeIds, const std::vector<std::vector<double>> & nodeLS)
{
  ThrowRequire(nodeIds.size() == nodeLS.size());
  const size_t numFields = lsFields.size();
  for (size_t i=0; i<nodeIds.size(); ++i)
  {
    const auto node = mesh.get_entity(stk::topology::NODE_RANK, nodeIds[i]);
    ThrowRequire(mesh.is_valid(node));
    ThrowRequire(numFields == nodeLS[i].size());
    for (size_t j=0; j<numFields; ++j)
      *field_data<double>(lsFields[j], node) = nodeLS[i][j];
  }
}

typedef CompleteDecompositionFixture<SimpleStkFixture2d, LSPerPhasePolicy<3> > CDMeshTests2DLSPerPhase;
TEST_F(CDMeshTests2DLSPerPhase, Two_Tri3_Check_Compatibility_When_Snapping)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_size = stk::parallel_machine_size(pm);
  const int parallel_rank = stk::parallel_machine_rank(pm);

  if (parallel_size > 2) return;

  cdfemSupport.set_cdfem_edge_degeneracy_handling(SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE);

  stk::mesh::BulkData & mesh = fixture.bulk_data();
  AuxMetaData & aux_meta = AuxMetaData::get(fixture.meta_data());

  aux_meta.set_assert_32bit_flag();
  aux_meta.clear_force_64bit_flag();

  setup_ls_field();

  auto & blockPart = declare_input_block("block_1", stk::topology::TRIANGLE_3_2D);
  register_ls_on_blocks({&blockPart});
  const auto & sideParts = declare_input_surfaces_touching_block(4, blockPart);

  commit();

  build_two_tri3_mesh_on_one_or_two_procs_with_sides(*this, blockPart, sideParts, parallel_rank);

  stk::mesh::create_exposed_block_boundary_sides(mesh, fixture.meta_data().universal_part(), {&aux_meta.exposed_boundary_part()});

  const double eps = 1.e-13;
  set_level_sets(mesh, ls_policy.ls_isovars, {1,2,3,4} , {{-1.,1.,1.+eps}, {-1.,1.,1.+eps}, {1.e2,1.,-1.}, {1.,-1.,-1.+eps}});

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

  //std::cout << log.get_log() << std::endl;

  std::vector<stk::mesh::Entity> entities;
  mesh.get_entities(stk::topology::NODE_RANK, fixture.meta_data().universal_part(), entities);
  EXPECT_EQ(7u, entities.size());
  mesh.get_entities(stk::topology::ELEMENT_RANK, aux_meta.get_part("block_1_P0"), entities);
  EXPECT_EQ(3u, entities.size());
  mesh.get_entities(stk::topology::ELEMENT_RANK, aux_meta.get_part("block_1_P1"), entities);
  EXPECT_EQ(1u, entities.size());
  mesh.get_entities(stk::topology::ELEMENT_RANK, aux_meta.get_part("block_1_P2"), entities);
  EXPECT_EQ(2u, entities.size());

  //fixture.write_results("Two_Tri3_Check_Compatibility_When_Snapping_LSPerPhase.e");
}

typedef CompleteDecompositionFixture<SimpleStkFixture3d, SingleLSPolicy> CDMeshTests3D;
TEST_F(CDMeshTests3D, Write_Results_No_Side)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_size = stk::parallel_machine_size(pm);
  const int parallel_rank = stk::parallel_machine_rank(pm);

  if (parallel_size > 2) return;

  const stk::topology tet4 = stk::topology::TETRAHEDRON_4;
  auto & block1_part = declare_input_block("block_1", tet4);
  auto & block2_part = declare_input_block("block_2", tet4);
  auto & surface_part = declare_input_surface("surface_1", tet4.side_topology(), {block1_part.mesh_meta_data_ordinal(), block2_part.mesh_meta_data_ordinal()});

  commit();

  build_two_tet4_mesh_np2(*this, block1_part, block2_part, surface_part, parallel_rank, parallel_size, false);

  fixture.write_results("Write_Results_No_Side.e");
}

TEST_F(CDMeshTests3D, Write_Results_With_Side)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_size = stk::parallel_machine_size(pm);
  const int parallel_rank = stk::parallel_machine_rank(pm);

  if (parallel_size > 2) return;

  const stk::topology tet4 = stk::topology::TETRAHEDRON_4;
  auto & block1_part = declare_input_block("block_1", tet4);
  auto & block2_part = declare_input_block("block_2", tet4);
  auto & surface_part = declare_input_surface("surface_1", tet4.side_topology(), {block1_part.mesh_meta_data_ordinal(), block2_part.mesh_meta_data_ordinal()});

  commit();

  build_two_tet4_mesh_np2(*this, block1_part, block2_part, surface_part, parallel_rank, parallel_size);

  fixture.write_results("Write_Results_With_Side.e");
}

namespace
{
void randomize_ls_field(const stk::mesh::BulkData & mesh,
    const FieldRef & field,
    std::mt19937 & mt,
    std::uniform_real_distribution<double> & dist)
{
  auto & node_buckets = mesh.buckets(stk::topology::NODE_RANK);
  for(auto && b_ptr : node_buckets)
  {
    for(auto && node : *b_ptr)
    {
      double * ls_data = field_data<double>(field, node);
      if(ls_data) *ls_data = dist(mt);
    }
  }
}
void set_ls_field_on_part(const stk::mesh::BulkData & mesh,
    const stk::mesh::Part & part,
    const FieldRef & ls_field,
    const double ls_value)
{
  auto & node_buckets = mesh.get_buckets(stk::topology::NODE_RANK, part);
  for(auto && b_ptr : node_buckets)
  {
    for(auto && node : *b_ptr)
    {
      double * ls_data = field_data<double>(ls_field, node);
      if(ls_data) *ls_data = ls_value;
    }
  }
}

}

TEST_F(CDMeshTests3D, Random_TwoTet4_InternalSideset)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_size = stk::parallel_machine_size(pm);
  const int parallel_rank = stk::parallel_machine_rank(pm);

  stk::mesh::BulkData & mesh = fixture.bulk_data();
  stk::mesh::MetaData & meta = fixture.meta_data();
  AuxMetaData & aux_meta = AuxMetaData::get(meta);


  if (parallel_size > 2) return;

  // Use a large snap tolerance to make snapping more common since it is a frequent source
  // of parallel bugs
  cdfemSupport.set_cdfem_edge_tol(0.1);
  cdfemSupport.set_simplex_generation_method(CUT_QUADS_BY_GLOBAL_IDENTIFIER);

  setup_ls_field();

  const stk::topology tet4 = stk::topology::TETRAHEDRON_4;
  auto & block1_part = declare_input_block("block_1", tet4);
  auto & block2_part = declare_input_block("block_2", tet4);
  auto & surface_part = declare_input_surface("surface_1", tet4.side_topology(), {block1_part.mesh_meta_data_ordinal(), block2_part.mesh_meta_data_ordinal()});

  register_ls_on_blocks({&block1_part, &block2_part});

  commit();

  build_two_tet4_mesh_np2(*this, block1_part, block2_part, surface_part, parallel_rank, parallel_size);

  stk::mesh::create_exposed_block_boundary_sides(mesh, meta.universal_part(), {&aux_meta.exposed_boundary_part()});

  std::mt19937 mt(std::mt19937::default_seed + parallel_rank);
  std::uniform_real_distribution<double> dist(-1., 1.);
#ifdef NDEBUG
  const int num_cases = 5000;
#else
  const int num_cases = 1000;
#endif
  for(int i=0; i < num_cases; ++i)
  {
    if (i%1000 == 0) std::cout << "Testing random configuration " << i << std::endl;

    randomize_ls_field(mesh, ls_policy.ls_isovar, mt, dist);
    stk::mesh::communicate_field_data(mesh, {&ls_policy.ls_isovar.field(), &coord_field.field()});

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

    if(HasNonfatalFailure())
    {
      std::cout << "Failure on iteration " << i << std::endl;
      krino_mesh->debug_output();
      std::cout << log.get_log() << std::endl;
      std::ostringstream fname;
      fname << "Random_TwoTet4_InternalSideset_iter_" << i << ".e";
      fixture.write_results(fname.str());
      ASSERT_TRUE(false);
    }

    krino_mesh = std::make_shared<CDMesh>(mesh, krino_mesh);
  }
}

typedef CompleteDecompositionFixture<BoundingBoxMeshTri3, SingleLSPolicy> CDMeshTestsBboxMesh2D;
TEST_F(CDMeshTestsBboxMesh2D, Random_SnapMesh)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_size = stk::parallel_machine_size(pm);
  const int parallel_rank = stk::parallel_machine_rank(pm);

  stk::mesh::MetaData & meta = fixture.meta_data();
  AuxMetaData & aux_meta = AuxMetaData::get(meta);

  if (parallel_size > 2) return;

  cdfemSupport.set_cdfem_edge_degeneracy_handling(SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE);
  const double approxMinRelativeSize = 0.25;

  setup_ls_field();

  auto & block1_part = aux_meta.get_part("block_1");

  register_ls_on_blocks({&block1_part});

  typename BoundingBoxMesh::BoundingBoxType domain(Vector3d::ZERO, Vector3d(1.,1.,0.));
  const double mesh_size = 1./3.;
  fixture.set_domain(domain, mesh_size);
  fixture.populate_mesh();
  stk::mesh::BulkData & mesh = fixture.bulk_data();
  stk::mesh::create_exposed_block_boundary_sides(mesh, meta.universal_part(), {&aux_meta.exposed_boundary_part(),&aux_meta.active_part()});

  double overallMinEdgeLength = std::numeric_limits<double>::max();
  double overallMaxEdgeLength = -std::numeric_limits<double>::max();
  double overallMinVolume = std::numeric_limits<double>::max();
  double overallMaxVolume = -std::numeric_limits<double>::max();
  const double expectedMinLength = mesh_size*approxMinRelativeSize;
  const double expectedMinVol = std::pow(mesh_size*approxMinRelativeSize, 2.) / 2.;

  std::mt19937 mt(std::mt19937::default_seed + parallel_rank);
  std::uniform_real_distribution<double> dist(-1., 1.);
#ifdef NDEBUG
  const int num_cases = 5000;
#else
  const int num_cases = 1000;
#endif
  for(int i=0; i < num_cases; ++i)
  {
    if (i%1000 == 0) std::cout << "Testing random configuration " << i << std::endl;

    MeshClone::stash_or_restore_mesh(mesh, 0); // restore original uncut mesh
    commit(); // new krino_mesh each time

    randomize_ls_field(mesh, ls_policy.ls_isovar, mt, dist);
    set_ls_field_on_part(mesh, aux_meta.exposed_boundary_part(), ls_policy.ls_isovar, 1.);
    stk::mesh::communicate_field_data(mesh, {&ls_policy.ls_isovar.field(), &coord_field.field()});

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
    if (minVolume < 0.5*expectedMinVol)
    {
      failedQuality = true;
      std::cout << "Failed quality requirements: minEdgeLength=" << minEdgeLength
      << ", maxEdgeLength=" << maxEdgeLength
      << ", minVolume=" << minVolume
      << ", maxVolume=" << maxVolume << std::endl;
    }

    if(HasNonfatalFailure() || failedQuality)
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
  std::cout << "Expected quality: minEdgeLength~=" << expectedMinLength << ", minVolume~=" << expectedMinVol << std::endl;
  std::cout << "Quality results: minEdgeLength=" << overallMinEdgeLength
      << ", maxEdgeLength=" << overallMaxEdgeLength
      << ", minVolume=" << overallMinVolume
      << ", maxVolume=" << overallMaxVolume << std::endl;
}

typedef CompleteDecompositionFixture<BoundingBoxMeshTri3, LSPerPhasePolicy<3>> CDMeshTestsBboxMesh2DLSPerPhase;
TEST_F(CDMeshTestsBboxMesh2DLSPerPhase, Random_SnapMesh)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_size = stk::parallel_machine_size(pm);
  const int parallel_rank = stk::parallel_machine_rank(pm);

  stk::mesh::MetaData & meta = fixture.meta_data();
  AuxMetaData & aux_meta = AuxMetaData::get(meta);

  if (parallel_size > 2) return;

  cdfemSupport.set_cdfem_edge_degeneracy_handling(SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE);
  const double approxMinRelativeSize = 0.25;

  setup_ls_field();

  auto & block1_part = aux_meta.get_part("block_1");

  register_ls_on_blocks({&block1_part});

  typename BoundingBoxMesh::BoundingBoxType domain(Vector3d::ZERO, Vector3d(1.,1.,0.));
  const double mesh_size = 1./3.;
  fixture.set_domain(domain, mesh_size);
  fixture.populate_mesh();
  fixture.create_domain_sides();
  stk::mesh::BulkData & mesh = fixture.bulk_data();
  stk::mesh::create_exposed_block_boundary_sides(mesh, meta.universal_part(), {&aux_meta.exposed_boundary_part(),&aux_meta.active_part()});

  const auto & ls_isovars = ls_policy.ls_isovars;
  std::vector<const stk::mesh::FieldBase *> sync_fields = {&coord_field.field()};
  for(auto && isovar : ls_isovars) sync_fields.push_back(&isovar.field());

  double overallMinEdgeLength = std::numeric_limits<double>::max();
  double overallMaxEdgeLength = -std::numeric_limits<double>::max();
  double overallMinVolume = std::numeric_limits<double>::max();
  double overallMaxVolume = -std::numeric_limits<double>::max();
  const double expectedMinLength = mesh_size*approxMinRelativeSize;
  const double expectedMinVol = std::pow(mesh_size*approxMinRelativeSize, 2.) / 2.;

  std::mt19937 mt(std::mt19937::default_seed + parallel_rank);
  std::uniform_real_distribution<double> dist(-1., 1.);
#ifdef NDEBUG
  const int num_cases = 5000;
#else
  const int num_cases = 1000;
#endif
  for(int i=0; i < num_cases; ++i)
  {
    if (i%1000 == 0) std::cout << "Testing random configuration " << i << std::endl;

    MeshClone::stash_or_restore_mesh(mesh, 0); // restore original uncut mesh
    commit(); // new krino_mesh each time

    for(auto && isovar : ls_isovars)
    {
      randomize_ls_field(mesh, isovar, mt, dist);
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

    if(HasNonfatalFailure())
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
  std::cout << "Expected quality: minEdgeLength~=" << expectedMinLength << ", minVolume~=" << expectedMinVol << std::endl;
  std::cout << "Quality results: minEdgeLength=" << overallMinEdgeLength
      << ", maxEdgeLength=" << overallMaxEdgeLength
      << ", minVolume=" << overallMinVolume
      << ", maxVolume=" << overallMaxVolume << std::endl;
}

typedef CompleteDecompositionFixture<BoundingBoxMeshTet4, SingleLSPolicy> CDMeshTestsBboxMesh3D;
TEST_F(CDMeshTestsBboxMesh3D, Random_SnapMesh)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_size = stk::parallel_machine_size(pm);
  const int parallel_rank = stk::parallel_machine_rank(pm);

  stk::mesh::MetaData & meta = fixture.meta_data();
  AuxMetaData & aux_meta = AuxMetaData::get(meta);

  if (parallel_size > 2) return;

  cdfemSupport.set_cdfem_edge_degeneracy_handling(SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE);
  const double approxMinRelativeSize = 0.25;

  setup_ls_field();

  auto & block1_part = aux_meta.get_part("block_1");

  register_ls_on_blocks({&block1_part});

  typename BoundingBoxMesh::BoundingBoxType domain(Vector3d::ZERO, Vector3d(1.,1.,1.));
  const double mesh_size = 1./3.;
  fixture.set_domain(domain, mesh_size);
  fixture.populate_mesh();
  stk::mesh::BulkData & mesh = fixture.bulk_data();
  stk::mesh::create_exposed_block_boundary_sides(mesh, meta.universal_part(), {&aux_meta.exposed_boundary_part(),&aux_meta.active_part()});

  double overallMinEdgeLength = std::numeric_limits<double>::max();
  double overallMaxEdgeLength = -std::numeric_limits<double>::max();
  double overallMinVolume = std::numeric_limits<double>::max();
  double overallMaxVolume = -std::numeric_limits<double>::max();
  const double expectedMinLength = mesh_size*approxMinRelativeSize;
  const double expectedMinVol = std::pow(mesh_size*approxMinRelativeSize, 3.) / 6.;

  std::mt19937 mt(std::mt19937::default_seed + parallel_rank);
  std::uniform_real_distribution<double> dist(-1., 1.);
#ifdef NDEBUG
  const int num_cases = 1000;
#else
  const int num_cases = 250;
#endif
  for(int i=0; i < num_cases; ++i)
  {
    if (i%250 == 0) std::cout << "Testing random configuration " << i << std::endl;

    MeshClone::stash_or_restore_mesh(mesh, 0); // restore original uncut mesh
    commit(); // new krino_mesh each time

    randomize_ls_field(mesh, ls_policy.ls_isovar, mt, dist);
    set_ls_field_on_part(mesh, aux_meta.exposed_boundary_part(), ls_policy.ls_isovar, 1.);
    stk::mesh::communicate_field_data(mesh, {&ls_policy.ls_isovar.field(), &coord_field.field()});

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
    if (minVolume < 0.5*expectedMinVol)
    {
      failedQuality = true;
      std::cout << "Failed quality requirements: minEdgeLength=" << minEdgeLength
      << ", maxEdgeLength=" << maxEdgeLength
      << ", minVolume=" << minVolume
      << ", maxVolume=" << maxVolume << std::endl;
    }

    if(HasNonfatalFailure() || failedQuality)
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
  std::cout << "Expected quality: minEdgeLength~=" << expectedMinLength << ", minVolume~=" << expectedMinVol << std::endl;
  std::cout << "Actual quality results: minEdgeLength=" << overallMinEdgeLength
      << ", maxEdgeLength=" << overallMaxEdgeLength
      << ", minVolume=" << overallMinVolume
      << ", maxVolume=" << overallMaxVolume << std::endl;
}

typedef CompleteDecompositionFixture<BoundingBoxMeshTet4, SingleLSPolicy> CDMeshTestsBboxMesh3DBCC;
TEST_F(CDMeshTestsBboxMesh3DBCC, Random_SnapMesh)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_size = stk::parallel_machine_size(pm);
  const int parallel_rank = stk::parallel_machine_rank(pm);

  stk::mesh::MetaData & meta = fixture.meta_data();
  AuxMetaData & aux_meta = AuxMetaData::get(meta);

  if (parallel_size > 2) return;

  cdfemSupport.set_cdfem_edge_degeneracy_handling(SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE);
  const double approxMinRelativeSize = 0.25;

  setup_ls_field();

  auto & block1_part = aux_meta.get_part("block_1");

  register_ls_on_blocks({&block1_part});

  typename BoundingBoxMesh::BoundingBoxType domain(Vector3d::ZERO, Vector3d(1.,1.,1.));
  const double mesh_size = 1./3.;
  fixture.set_domain(domain, mesh_size);
  fixture.set_mesh_structure_type(BCC_BOUNDING_BOX_MESH);
  fixture.populate_mesh();
  stk::mesh::BulkData & mesh = fixture.bulk_data();
  stk::mesh::create_exposed_block_boundary_sides(mesh, meta.universal_part(), {&aux_meta.exposed_boundary_part(),&aux_meta.active_part()});

  const double BCCsize = std::sqrt(3.)/2.*mesh_size;
  double overallMinEdgeLength = std::numeric_limits<double>::max();
  double overallMaxEdgeLength = -std::numeric_limits<double>::max();
  double overallMinVolume = std::numeric_limits<double>::max();
  double overallMaxVolume = -std::numeric_limits<double>::max();
  const double expectedMinLength = BCCsize*approxMinRelativeSize;
  const double expectedMinVol = std::pow(BCCsize*approxMinRelativeSize, 3.) / (6.*std::sqrt(2.));

  std::mt19937 mt(std::mt19937::default_seed + parallel_rank);
  std::uniform_real_distribution<double> dist(-1., 1.);
#ifdef NDEBUG
  const int num_cases = 1000;
#else
  const int num_cases = 250;
#endif
  for(int i=0; i < num_cases; ++i)
  {
    if (i%250 == 0) std::cout << "Testing random configuration " << i << std::endl;

    MeshClone::stash_or_restore_mesh(mesh, 0); // restore original uncut mesh
    commit(); // new krino_mesh each time

    randomize_ls_field(mesh, ls_policy.ls_isovar, mt, dist);
    set_ls_field_on_part(mesh, aux_meta.exposed_boundary_part(), ls_policy.ls_isovar, 1.);
    stk::mesh::communicate_field_data(mesh, {&ls_policy.ls_isovar.field(), &coord_field.field()});

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
    if (minVolume < 0.5*expectedMinVol)
    {
      failedQuality = true;
      std::cout << "Failed quality requirements: minEdgeLength=" << minEdgeLength
      << ", maxEdgeLength=" << maxEdgeLength
      << ", minVolume=" << minVolume
      << ", maxVolume=" << maxVolume << std::endl;
    }

    if(HasNonfatalFailure() || failedQuality)
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
  std::cout << "Expected quality: minEdgeLength~=" << expectedMinLength << ", minVolume~=" << expectedMinVol << std::endl;
  std::cout << "Actual quality results: minEdgeLength=" << overallMinEdgeLength
      << ", maxEdgeLength=" << overallMaxEdgeLength
      << ", minVolume=" << overallMinVolume
      << ", maxVolume=" << overallMaxVolume << std::endl;
}

TEST_F(CDMeshTests2DLSPerPhase, Tri3_3LS_SnappedTriplePoint)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_size = stk::parallel_machine_size(pm);

  if (parallel_size > 1) return;

  stk::mesh::BulkData & mesh = fixture.bulk_data();
  stk::mesh::MetaData & meta = fixture.meta_data();
  AuxMetaData & aux_meta = AuxMetaData::get(meta);

  cdfemSupport.set_cdfem_edge_tol(0.15);

  setup_ls_field();

  auto & block_part = declare_input_block("block_1", stk::topology::TRIANGLE_3_2D);
  register_ls_on_blocks({&block_part});

  commit();

  mesh.modification_begin();
  {
    stk::mesh::PartVector elem_parts;
    elem_parts.push_back(&block_part);
    elem_parts.push_back(&aux_meta.active_part());
    std::vector<stk::mesh::EntityId> node_ids = {1, 2, 3};

    create_element(elem_parts, 1, node_ids);
  }
  mesh.modification_end();

  stk::mesh::create_exposed_block_boundary_sides(mesh, meta.universal_part(), {&aux_meta.exposed_boundary_part()});

  const auto node1 = mesh.get_entity(stk::topology::NODE_RANK, 1);
  const auto node2 = mesh.get_entity(stk::topology::NODE_RANK, 2);
  const auto node3 = mesh.get_entity(stk::topology::NODE_RANK, 3);

  double * node1_coords = field_data<double>(coord_field, node1);
  double * node2_coords = field_data<double>(coord_field, node2);
  double * node3_coords = field_data<double>(coord_field, node3);

  node1_coords[0] = 0.;
  node1_coords[1] = 0.;

  node2_coords[0] = 1.;
  node2_coords[1] = 0.;

  node3_coords[0] = 0.;
  node3_coords[1] = 1.;

  const auto & ls_isovars = ls_policy.ls_isovars;
  *field_data<double>(ls_isovars[0], node1) = 0.;
  *field_data<double>(ls_isovars[0], node2) = 0.;
  *field_data<double>(ls_isovars[0], node3) = 0.;

  *field_data<double>(ls_isovars[1], node1) = 0.1;
  *field_data<double>(ls_isovars[1], node2) = -0.2;
  *field_data<double>(ls_isovars[1], node3) = -0.01;

  *field_data<double>(ls_isovars[2], node1) = 0.2;
  *field_data<double>(ls_isovars[2], node2) = -0.25;
  *field_data<double>(ls_isovars[2], node3) = -0.005;

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

  std::vector<stk::mesh::Entity> entities;

  // This assumes that the element is cut first by the (0,1) interface,
  // then the (0, 2) virtual interface, and then the (1, 2) interface.
  // THere are other valid decompositions if the element is cut in a different order.
  // Should have added 3 nodes at the cutting locations
  mesh.get_entities(stk::topology::NODE_RANK, meta.universal_part(), entities);
  EXPECT_EQ(5u, entities.size());
  mesh.get_entities(stk::topology::NODE_RANK, aux_meta.active_part(), entities);
  EXPECT_EQ(5u, entities.size());

  mesh.get_entities(stk::topology::EDGE_RANK, aux_meta.active_part() &
      aux_meta.get_part("surface_block_1_P0_P1") & aux_meta.get_part("surface_block_1_P1_P0"), entities);
  EXPECT_EQ(1u, entities.size());
  mesh.get_entities(stk::topology::EDGE_RANK, aux_meta.active_part() &
      aux_meta.get_part("surface_block_1_P1_P2") & aux_meta.get_part("surface_block_1_P2_P1"), entities);
  EXPECT_EQ(1u, entities.size());
  mesh.get_entities(stk::topology::EDGE_RANK, aux_meta.active_part() &
      aux_meta.get_part("surface_block_1_P0_P2") & aux_meta.get_part("surface_block_1_P2_P0"), entities);
  EXPECT_EQ(0u, entities.size());

  // Should be 3 conformal elements plus the parent element
  mesh.get_entities(stk::topology::ELEMENT_RANK, meta.universal_part(), entities);
  EXPECT_EQ(4u, entities.size());
  mesh.get_entities(stk::topology::ELEMENT_RANK, aux_meta.active_part(), entities);
  EXPECT_EQ(3u, entities.size());
  mesh.get_entities(stk::topology::ELEMENT_RANK, !aux_meta.active_part(), entities);
  EXPECT_EQ(1u, entities.size());
  mesh.get_entities(stk::topology::ELEMENT_RANK, aux_meta.get_part("block_1_P0"), entities);
  EXPECT_EQ(1u, entities.size());
  mesh.get_entities(stk::topology::ELEMENT_RANK, aux_meta.get_part("block_1_P1"), entities);
  EXPECT_EQ(1u, entities.size());
  mesh.get_entities(stk::topology::ELEMENT_RANK, aux_meta.get_part("block_1_P2"), entities);
  EXPECT_EQ(1u, entities.size());
}

TEST_F(CDMeshTests2DLSPerPhase, Tri3_3LS_TriplePointDebug)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_size = stk::parallel_machine_size(pm);

  if (parallel_size > 1) return;

  stk::mesh::BulkData & mesh = fixture.bulk_data();
  stk::mesh::MetaData & meta = fixture.meta_data();
  AuxMetaData & aux_meta = AuxMetaData::get(meta);

  cdfemSupport.set_cdfem_edge_tol(0.01);

  setup_ls_field();

  auto & block_part = declare_input_block("block_1", stk::topology::TRIANGLE_3_2D);
  register_ls_on_blocks({&block_part});

  commit();

  mesh.modification_begin();
  {
    stk::mesh::PartVector elem_parts;
    elem_parts.push_back(&block_part);
    elem_parts.push_back(&aux_meta.active_part());
    std::vector<stk::mesh::EntityId> node_ids = {1, 2, 3};

    create_element(elem_parts, 1, node_ids);
  }
  mesh.modification_end();

  stk::mesh::create_exposed_block_boundary_sides(mesh, meta.universal_part(), {&aux_meta.exposed_boundary_part()});

  const auto node1 = mesh.get_entity(stk::topology::NODE_RANK, 1);
  const auto node2 = mesh.get_entity(stk::topology::NODE_RANK, 2);
  const auto node3 = mesh.get_entity(stk::topology::NODE_RANK, 3);

  double * node1_coords = field_data<double>(coord_field, node1);
  double * node2_coords = field_data<double>(coord_field, node2);
  double * node3_coords = field_data<double>(coord_field, node3);

  node1_coords[0] = 0.;
  node1_coords[1] = 0.;

  node2_coords[0] = 1.;
  node2_coords[1] = 0.;

  node3_coords[0] = 0.;
  node3_coords[1] = 1.;

  const auto & ls_isovars = ls_policy.ls_isovars;
  *field_data<double>(ls_isovars[0], node1) = -0.2;
  *field_data<double>(ls_isovars[1], node1) =  0.1;
  *field_data<double>(ls_isovars[2], node1) =  0.6;

  *field_data<double>(ls_isovars[0], node2) =  0.7;
  *field_data<double>(ls_isovars[1], node2) = -0.5;
  *field_data<double>(ls_isovars[2], node2) =  0.4;

  *field_data<double>(ls_isovars[0], node3) =  0.1;
  *field_data<double>(ls_isovars[1], node3) =  0.3;
  *field_data<double>(ls_isovars[2], node3) = -0.5;

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

  // Test output, but remove output unless actually debugging.
  fixture.write_results("debug_2d.e");
  std::remove("debug_2d.e");
}

typedef CompleteDecompositionFixture<SimpleStkFixture2d, LSPerPhasePolicy<3> > CDMeshTests2DLSPerPhase;
TEST_F(CDMeshTests2DLSPerPhase, Random_TwoTri3_InternalSideset_Snap)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_size = stk::parallel_machine_size(pm);
  const int parallel_rank = stk::parallel_machine_rank(pm);

  if (parallel_size > 2) return;

  stk::mesh::BulkData & mesh = fixture.bulk_data();
  stk::mesh::MetaData & meta = fixture.meta_data();
  AuxMetaData & aux_meta = AuxMetaData::get(meta);

  cdfemSupport.set_cdfem_edge_degeneracy_handling(SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE);

  setup_ls_field();

  const stk::topology tri3 = stk::topology::TRIANGLE_3_2D;
  auto & block1_part = declare_input_block("block_1", tri3);
  auto & block2_part = declare_input_block("block_2", tri3);
  auto & surface_part = declare_input_surface("surface_1", tri3.side_topology(), {block1_part.mesh_meta_data_ordinal(), block2_part.mesh_meta_data_ordinal()});

  register_ls_on_blocks({&block1_part, &block2_part});

  commit();

  build_two_tri3_mesh_np2(*this, block1_part, block2_part, surface_part, parallel_rank, parallel_size);

  stk::mesh::create_exposed_block_boundary_sides(mesh, meta.universal_part(), {&aux_meta.exposed_boundary_part()});

  const auto & ls_isovars = ls_policy.ls_isovars;
  std::vector<const stk::mesh::FieldBase *> sync_fields = {&coord_field.field()};
  for(auto && isovar : ls_isovars) sync_fields.push_back(&isovar.field());

  std::mt19937 mt(std::mt19937::default_seed + parallel_rank);
  std::uniform_real_distribution<double> dist(-1., 1.);
#ifdef NDEBUG
  const int num_cases = 5000;
#else
  const int num_cases = 1000;
#endif
  for(int i=0; i < num_cases; ++i)
  {
    if (i%1000 == 0) std::cout << "Testing random configuration " << i << std::endl;

    MeshClone::stash_or_restore_mesh(mesh, 0); // restore original uncut mesh
    commit(); // new krino_mesh each time

    for(auto && isovar : ls_isovars)
    {
      randomize_ls_field(mesh, isovar, mt, dist);
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

    if (false)
    {
      std::ostringstream fname;
      fname << "Random_TwoTri3_InternalSideset_Snap_iter_" << i << ".e";
      fixture.write_results(fname.str());
    }

    if(HasNonfatalFailure())
    {
      std::cout << "Failure on iteration " << i << std::endl;
      std::ostringstream fname;
      fname << "Random_TwoTri3_InternalSideset_Snap_iter_" << i << ".e";
      fixture.write_results(fname.str());
      ASSERT_TRUE(false);
    }

    krino_mesh = std::make_shared<CDMesh>(mesh, krino_mesh);
  }
}

typedef CompleteDecompositionFixture<SimpleStkFixture3d, LSPerInterfacePolicy<3> > CDMeshTests3DLSPerInterface;
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

  setup_ls_field();

  const stk::topology tet4 = stk::topology::TETRAHEDRON_4;
  auto & block1_part = declare_input_block("block_1", tet4);

  register_ls_on_blocks({&block1_part});

  commit();

  build_one_tet4_mesh(*this, block1_part, parallel_rank, parallel_size);

  stk::mesh::create_exposed_block_boundary_sides(mesh, meta.universal_part(), {&aux_meta.exposed_boundary_part()});

  const auto & ls_isovars = ls_policy.ls_isovars;
  std::vector<const stk::mesh::FieldBase *> sync_fields = {&coord_field.field()};
  for(auto && isovar : ls_isovars) sync_fields.push_back(&isovar.field());

  std::mt19937 mt(std::mt19937::default_seed + parallel_rank);
  std::uniform_real_distribution<double> dist(-1., 1.);
#ifdef NDEBUG
  const int num_cases = 5000;
#else
  const int num_cases = 1000;
#endif
  for(int i=0; i < num_cases; ++i)
  {
    if (i%1000 == 0) std::cout << "Testing random configuration " << i << std::endl;

    MeshClone::stash_or_restore_mesh(mesh, 0); // restore original uncut mesh
    commit(); // new krino_mesh each time

    for(auto && isovar : ls_isovars)
    {
      randomize_ls_field(mesh, isovar, mt, dist);
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

    if(HasNonfatalFailure())
    {
      std::cout << log.get_log() << std::endl;
      std::cout << "Failure on iteration " << i << std::endl;
      std::ostringstream fname;
      fname << "Random_TwoTet4_InternalSideset_iter_" << i << ".e";
      fixture.write_results(fname.str());
      ASSERT_TRUE(false);
    }

    krino_mesh = std::make_shared<CDMesh>(mesh, krino_mesh);
  }
}

typedef CompleteDecompositionFixture<SimpleStkFixture3d, LSPerInterfacePolicy<3> > CDMeshTests3DLSPerInterface;
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

  setup_ls_field();

  const stk::topology tet4 = stk::topology::TETRAHEDRON_4;
  auto & block1_part = declare_input_block("block_1", tet4);

  register_ls_on_blocks({&block1_part});

  commit();

  build_one_tet4_mesh(*this, block1_part, parallel_rank, parallel_size);

  stk::mesh::create_exposed_block_boundary_sides(mesh, meta.universal_part(), {&aux_meta.exposed_boundary_part()});

  const std::array<stk::mesh::Entity,4> nodes = {{ mesh.get_entity(stk::topology::NODE_RANK, 1), mesh.get_entity(stk::topology::NODE_RANK, 2), mesh.get_entity(stk::topology::NODE_RANK, 3), mesh.get_entity(stk::topology::NODE_RANK, 4) }};
  const auto & ls_isovars = ls_policy.ls_isovars;
  *field_data<double>(ls_isovars[0], nodes[0]) = -1;
  *field_data<double>(ls_isovars[0], nodes[1]) = 1.;
  *field_data<double>(ls_isovars[0], nodes[2]) = 2.;
  *field_data<double>(ls_isovars[0], nodes[3]) = 3.;

  commit();
  decompose_mesh();


  EXPECT_TRUE(check_induced_parts(mesh));
  EXPECT_TRUE(check_face_and_edge_ownership(mesh));
  EXPECT_TRUE(check_face_and_edge_relations(mesh));
  EXPECT_TRUE(check_shared_entity_nodes(mesh));
  EXPECT_TRUE(krino_mesh->check_element_side_parts());

  EXPECT_EQ(1u+1u, mesh.num_elements(nodes[0]));
  EXPECT_EQ(1u+1u, mesh.num_elements(nodes[1]));
  EXPECT_EQ(2u+1u, mesh.num_elements(nodes[2]));
  EXPECT_EQ(3u+1u, mesh.num_elements(nodes[3]));

  // regression test
  const ScaledJacobianQualityMetric qualityMetric;
  const double quality = determine_quality(mesh, krino_mesh->get_active_part(), qualityMetric);
  const double goldQuality = 0.21;
  EXPECT_GT(quality, goldQuality);
}

typedef CompleteDecompositionFixture<SimpleStkFixture3d, LSPerPhasePolicy<3> > CDMeshTests3DLSPerPhase;
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

  setup_ls_field();

  const stk::topology tet4 = stk::topology::TETRAHEDRON_4;
  auto & block1_part = declare_input_block("block_1", tet4);
  auto & block2_part = declare_input_block("block_2", tet4);
  auto & surface_part = declare_input_surface("surface_1", tet4.side_topology(), {block1_part.mesh_meta_data_ordinal(), block2_part.mesh_meta_data_ordinal()});

  register_ls_on_blocks({&block1_part, &block2_part});

  commit();

  build_two_tet4_mesh_np2(*this, block1_part, block2_part, surface_part, parallel_rank, parallel_size);

  stk::mesh::create_exposed_block_boundary_sides(mesh, meta.universal_part(), {&aux_meta.exposed_boundary_part()});

  const auto & ls_isovars = ls_policy.ls_isovars;
  std::vector<const stk::mesh::FieldBase *> sync_fields = {&coord_field.field()};
  for(auto && isovar : ls_isovars) sync_fields.push_back(&isovar.field());

  std::mt19937 mt(std::mt19937::default_seed + parallel_rank);
  std::uniform_real_distribution<double> dist(-1., 1.);
#ifdef NDEBUG
  const int num_cases = 5000;
#else
  const int num_cases = 1000;
#endif
  for(int i=0; i < num_cases; ++i)
  {
    if (i%1000 == 0) std::cout << "Testing random configuration " << i << std::endl;

    for(auto && isovar : ls_isovars)
    {
      randomize_ls_field(mesh, isovar, mt, dist);
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

    if(HasNonfatalFailure())
    {
      std::cout << "Failure on iteration " << i << std::endl;
      std::ostringstream fname;
      fname << "Random_TwoTet4_InternalSideset_iter_" << i << ".e";
      fixture.write_results(fname.str());
      ASSERT_TRUE(false);
    }

    krino_mesh = std::make_shared<CDMesh>(mesh, krino_mesh);
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

  cdfemSupport.set_cdfem_edge_degeneracy_handling(SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE);

  setup_ls_field();

  const stk::topology tet4 = stk::topology::TETRAHEDRON_4;
  auto & block1_part = declare_input_block("block_1", tet4);
  auto & block2_part = declare_input_block("block_2", tet4);
  auto & surface_part = declare_input_surface("surface_1", tet4.side_topology(), {block1_part.mesh_meta_data_ordinal(), block2_part.mesh_meta_data_ordinal()});

  register_ls_on_blocks({&block1_part, &block2_part});

  commit();

  build_two_tet4_mesh_np2(*this, block1_part, block2_part, surface_part, parallel_rank, parallel_size);

  stk::mesh::create_exposed_block_boundary_sides(mesh, meta.universal_part(), {&aux_meta.exposed_boundary_part()});

  const auto & ls_isovars = ls_policy.ls_isovars;
  std::vector<const stk::mesh::FieldBase *> sync_fields = {&coord_field.field()};
  for(auto && isovar : ls_isovars) sync_fields.push_back(&isovar.field());

  std::mt19937 mt(std::mt19937::default_seed + parallel_rank);
  std::uniform_real_distribution<double> dist(-1., 1.);
#ifdef NDEBUG
  const int num_cases = 5000;
#else
  const int num_cases = 1000;
#endif
  for(int i=0; i < num_cases; ++i)
  {
    if (i%1000 == 0) std::cout << "Testing random configuration " << i << std::endl;

    MeshClone::stash_or_restore_mesh(mesh, 0); // restore original uncut mesh
    commit(); // new krino_mesh each time

    for(auto && isovar : ls_isovars)
    {
      randomize_ls_field(mesh, isovar, mt, dist);
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

    if(HasNonfatalFailure())
    {
      std::cout << "Failure on iteration " << i << std::endl;
      std::ostringstream fname;
      fname << "Random_TwoTet4_InternalSideset_iter_" << i << ".e";
      fixture.write_results(fname.str());
      ASSERT_TRUE(false);
    }

    krino_mesh = std::make_shared<CDMesh>(mesh, krino_mesh);
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

  setup_ls_field();

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

  const auto & ls_isovars = ls_policy.ls_isovars;
  std::vector<const stk::mesh::FieldBase *> sync_fields = {&coord_field.field()};
  for (auto && isovar : ls_isovars)
    sync_fields.push_back(&isovar.field());

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
      for (auto && isovar : ls_isovars)
      {
        randomize_ls_field(mesh, isovar, mt, dist);
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
    krino_mesh = std::make_shared<CDMesh>(mesh, krino_mesh);

    // Step 2.
    randomize_all_fields();
    run_decomp();
    krino::MeshClone::stash_or_restore_mesh(mesh, i);

    // Step 3.
    ASSERT_NO_THROW(krino_mesh->rebuild_after_rebalance());
    krino_mesh = std::make_shared<CDMesh>(mesh, krino_mesh);
  }
}

TEST_F(CDMeshTests3D, Rebalance_with_rcb)
{
  run_rebalance_with("rcb");
}

TEST_F(CDMeshTests3D, Rebalance_with_parmetis)
{
  run_rebalance_with("parmetis");
}

typedef CompleteDecompositionFixture<SimpleStkFixture3d, SingleLSPolicy> NonconformalAdaptivityTest;
TEST_F(NonconformalAdaptivityTest, InternalSidePositivePermutationNonOwnedElement)
{
  /* This tests that percept can correctly handle an internal side where the element
   * with the same owning proc as the side has a negative permutation and the second
   * connected element (with a different owning processor) has a positive permutation.
   * In serial it just tests use of the adaptivity interface.
   */

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_size = stk::parallel_machine_size(pm);
  const int parallel_rank = stk::parallel_machine_rank(pm);
  stk::mesh::BulkData & mesh = fixture.bulk_data();

  if (parallel_size > 2) return;

  cdfemSupport.set_cdfem_edge_tol(0.1);
  cdfemSupport.set_simplex_generation_method(CUT_QUADS_BY_GLOBAL_IDENTIFIER);

  setup_ls_field();

  const stk::topology tet4 = stk::topology::TETRAHEDRON_4;
  auto & block1_part = declare_input_block("block_1", tet4);
  auto & block2_part = declare_input_block("block_2", tet4);
  auto & surface_part = declare_input_surface("surface_1", tet4.side_topology(), {block1_part.mesh_meta_data_ordinal(), block2_part.mesh_meta_data_ordinal()});

  register_ls_on_blocks({&block1_part, &block2_part});

  auto & aux_meta = AuxMetaData::get(fixture.meta_data());
  FieldRef marker_field =
      aux_meta.register_field("refine_marker", FieldType::INTEGER, stk::topology::ELEMENT_RANK,
          1u, 1, fixture.meta_data().universal_part());

  auto & meta = fixture.meta_data();
  auto & active_part = aux_meta.active_part();
  stk::diag::TimerSet enabledTimerSet(0);
  stk::diag::Timer root_timer = createRootTimer("test", enabledTimerSet);
  HAdapt::setup(meta, active_part, root_timer);

  commit();

  const bool build_all_on_P0 = false;
  build_two_tet4_mesh_np2(*this, block1_part, block2_part, surface_part, parallel_rank,
      parallel_size, true, build_all_on_P0);

  auto elem_1 = mesh.get_entity(stk::topology::ELEMENT_RANK, 1u);
  int & elem_1_marker = *field_data<int>(marker_field, elem_1);
  elem_1_marker = 1;

  if(parallel_size == 2)
  {
    stk::mesh::EntityProcVec changes;
    if(parallel_rank == 0)
    {
      changes.push_back(std::make_pair(mesh.get_entity(stk::topology::FACE_RANK, 7), 1));
    }
    mesh.change_entity_owner(changes);
  }

  EXPECT_NO_THROW(HAdapt::do_adaptive_refinement(meta, marker_field.name()));
}

typedef CompleteDecompositionFixture<SimpleStkFixture3d, SingleLSPolicy> MeshCloneTest;
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
  auto & surface_1_part = declare_input_surface("surface_1", tet4.side_topology(), {block1_part.mesh_meta_data_ordinal(), block2_part.mesh_meta_data_ordinal()});
  auto & surface_2_part = declare_input_surface("surface_2", tet4.side_topology(), {block1_part.mesh_meta_data_ordinal(), block2_part.mesh_meta_data_ordinal()});

  auto & aux_meta = AuxMetaData::get(fixture.meta_data());
  aux_meta.register_field("side_field", FieldType::REAL, stk::topology::FACE_RANK,
      1u, 1u, surface_1_part);

  commit();

  build_two_tet4_mesh_np2(*this, block1_part, block2_part, surface_1_part, parallel_rank, parallel_size);

  // Stash the initial mesh.
  ASSERT_NO_THROW(MeshClone::stash_or_restore_mesh(mesh, 0));

  // Change the parallel owner of the shared face, then change its surface part from surface_1 to
  // surface_2.
  const auto side_1 = mesh.get_entity(stk::topology::FACE_RANK, 7);
  const auto current_owner = mesh.parallel_owner_rank(side_1);
  stk::mesh::EntityProcVec owner_changes;
  if(parallel_rank == current_owner)
  {
    owner_changes.emplace_back(side_1, (current_owner + 1) % 2);
  }
  mesh.change_entity_owner(owner_changes);
  ASSERT_TRUE(mesh.is_valid(side_1));
  const auto new_owner = mesh.parallel_owner_rank(side_1);
  ASSERT_NE(new_owner, current_owner);

  mesh.modification_begin();
  if(new_owner == parallel_rank)
  {
    mesh.change_entity_parts(side_1, stk::mesh::ConstPartVector{&surface_2_part}, stk::mesh::ConstPartVector{&surface_1_part});
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
  if(new_owner == parallel_rank)
  {
    mesh.change_entity_parts(side_1, stk::mesh::ConstPartVector{&surface_1_part}, stk::mesh::ConstPartVector{&surface_2_part});
  }
  mesh.modification_end();
  ASSERT_TRUE(mesh.bucket(side_1).member(surface_1_part));
  MeshClone::stash_or_restore_mesh(mesh, 1);
  const auto side_1_new = mesh.get_entity(stk::topology::FACE_RANK, 7);
  EXPECT_TRUE(mesh.bucket(side_1_new).member(surface_2_part));
  EXPECT_FALSE(mesh.bucket(side_1_new).member(surface_1_part));
}

}
