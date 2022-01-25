// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_Unit_Part_Decomposition_Fixture.hpp>
#include <Akri_AuxMetaData.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_Interface_Name_Generator.hpp>

namespace krino
{

Part_Decomposition_Fixture::Part_Decomposition_Fixture()
: fixture(3)
{
  AuxMetaData & aux_meta = AuxMetaData::get(get_meta_data());
  aux_meta.declare_io_part_with_topology("block_1", stk::topology::TETRAHEDRON_4);
  aux_meta.declare_io_part_with_topology("block_2", stk::topology::TETRAHEDRON_4);
  aux_meta.declare_io_part_with_topology("surface_1", stk::topology::TRIANGLE_3);
}

Part_Decomposition_Fixture::~Part_Decomposition_Fixture()
{
}

Block_Surface_Connectivity Part_Decomposition_Fixture::addOneSidedSideset()
{
  Block_Surface_Connectivity block_surface_info;
  stk::mesh::PartOrdinal surface_1_ordinal = get_meta_data().get_part("surface_1")->mesh_meta_data_ordinal();
  stk::mesh::PartOrdinal block_1_ordinal = get_meta_data().get_part("block_1")->mesh_meta_data_ordinal();
  block_surface_info.add_surface(surface_1_ordinal, {block_1_ordinal});
  return block_surface_info;
}

Block_Surface_Connectivity Part_Decomposition_Fixture::addTwoSidedSideset()
{
  Block_Surface_Connectivity block_surface_info;
  stk::mesh::PartOrdinal block_1_ordinal = get_meta_data().get_part("block_1")->mesh_meta_data_ordinal();
  stk::mesh::PartOrdinal block_2_ordinal = get_meta_data().get_part("block_2")->mesh_meta_data_ordinal();

  stk::mesh::Part & surface_1 = *get_meta_data().get_part("surface_1");
  stk::mesh::Part & surface_1_block_1 = get_meta_data().declare_part_with_topology("surface_block_1_tri3_1", stk::topology::TRIANGLE_3);
  stk::mesh::Part & surface_1_block_2 = get_meta_data().declare_part_with_topology("surface_block_2_tri3_1", stk::topology::TRIANGLE_3);
  get_meta_data().declare_part_subset(surface_1, surface_1_block_1);
  get_meta_data().declare_part_subset(surface_1, surface_1_block_2);
  block_surface_info.add_surface(surface_1.mesh_meta_data_ordinal(), {block_1_ordinal, block_2_ordinal});
  block_surface_info.add_surface(surface_1_block_1.mesh_meta_data_ordinal(), {block_1_ordinal});
  block_surface_info.add_surface(surface_1_block_2.mesh_meta_data_ordinal(), {block_2_ordinal});
  return block_surface_info;
}

stk::mesh::MetaData & Part_Decomposition_Fixture::get_meta_data()
{
  return fixture.meta_data();
}

stk::mesh::Part * Part_Decomposition_Fixture::findPart(const std::string & part_name)
{
  stk::mesh::Part * result = nullptr;
  const stk::mesh::PartVector & mesh_parts = fixture.meta_data().get_parts();
  stk::mesh::PartVector::const_iterator found;
  found = std::find_if(mesh_parts.begin(), mesh_parts.end(), PartNameIs(part_name));
  if( found != mesh_parts.end() )
  {
    result = *found;
  }
  return result;
}

stk::mesh::Part * Part_Decomposition_Fixture::findSuperset(const std::string & superset_name, const stk::mesh::Part * const part)
{
  stk::mesh::Part * result = nullptr;
  stk::mesh::PartVector::const_iterator found;
  found = std::find_if(part->supersets().begin(), part->supersets().end(), PartNameIs(superset_name));
  if( found != part->supersets().end() )
  {
    result = *found;
  }
  return result;
}

PhaseVec Part_Decomposition_Fixture::ls_phases(int num_ls, bool one_phase_per_ls)
{
  static bool init = false;
  static std::vector<PhaseTag> phase_tags;
  static std::vector<NamedPhase> named_phases;
  if(!init)
  {
    const LevelSet_Identifier id0(0);
    const LevelSet_Identifier id1(1);
    const LevelSet_Identifier id2(2);
    const LevelSet_Identifier id3(3);
    PhaseTag pp, nn, pn, np;
    pp.add(id0,1); pp.add(id1,1);
    nn.add(id0,-1); nn.add(id1,-1);
    pn.add(id0,1); pn.add(id1,-1);
    np.add(id0,-1); np.add(id1,1);

    phase_tags.push_back(pp);
    phase_tags.push_back(nn);
    phase_tags.push_back(pn);
    phase_tags.push_back(np);

    named_phases.push_back(NamedPhase("A", pp));
    named_phases.push_back(NamedPhase("B", nn));
    named_phases.push_back(NamedPhase("C", pn));
    named_phases.push_back(NamedPhase("D", np));

    PhaseTag ls1, ls2, ls3, ls4;
    ls1.add(id0,-1);
    ls2.add(id1,-1);
    ls3.add(id2,-1);
    ls4.add(id3,-1);
    named_phases.push_back(NamedPhase("LS1", ls1));
    named_phases.push_back(NamedPhase("LS2", ls2));
    named_phases.push_back(NamedPhase("LS3", ls3));
    named_phases.push_back(NamedPhase("LS4", ls4));

    init = true;
  }

  PhaseVec result;
  if(!one_phase_per_ls)
  {
    result.push_back(named_phases[0]);
    result.push_back(named_phases[1]);
    if(num_ls == 2)
    {
      result.push_back(named_phases[2]);
      result.push_back(named_phases[3]);
    }
  }
  else
  {
    for(int i=0; i < num_ls; ++i)
    {
      result.push_back(named_phases[4+i]);
    }
  }
  return result;
}

PhaseVec Part_Decomposition_Fixture::death_phases()
{
  static bool init = false;
  const LevelSet_Identifier id0(0);
  static PhaseTag pos, neg;
  pos.add(id0,1); neg.add(id0,-1);
  static NamedPhase dead("dead", pos), alive("", neg);
  static PhaseVec named_phases;
  if(!init)
  {
    named_phases.push_back(alive);
    named_phases.push_back(dead);
    init = true;
  }

  return named_phases;
}

const Interface_Name_Generator & Part_Decomposition_Fixture::ls_name_generator()
{
  static LS_Name_Generator name_gen;
  return name_gen;
}

const Interface_Name_Generator & Part_Decomposition_Fixture::death_name_generator()
{
  static Death_Name_Generator name_gen("test");
  return name_gen;
}

void Part_Decomposition_Fixture::performDecomposition(const stk::mesh::PartVector & used_blocks,
    const Block_Surface_Connectivity & input_block_surface_info,
    bool cdfem_death, int num_ls, bool one_ls_per_phase)
{
  Phase_Support & phase_support = Phase_Support::get(get_meta_data());
  phase_support.set_input_block_surface_connectivity(input_block_surface_info);
  std::vector<std::tuple<stk::mesh::PartVector, 
      std::shared_ptr<Interface_Name_Generator>, PhaseVec>> ls_sets;
  if(cdfem_death)
  {

    auto interface_name_gen = std::shared_ptr<Interface_Name_Generator>(new Death_Name_Generator("test"));
    ls_sets.push_back(std::make_tuple(used_blocks,interface_name_gen,Part_Decomposition_Fixture::death_phases()));
    phase_support.decompose_blocks(ls_sets);
    phase_support.build_decomposed_block_surface_connectivity();
  }
  else
  {
    auto interface_name_gen = std::shared_ptr<Interface_Name_Generator>(new LS_Name_Generator());
    ls_sets.push_back(std::make_tuple(used_blocks,interface_name_gen,Part_Decomposition_Fixture::ls_phases(num_ls, one_ls_per_phase)));
    phase_support.decompose_blocks(ls_sets);
    phase_support.build_decomposed_block_surface_connectivity();
  }
}

} // namespace krino
