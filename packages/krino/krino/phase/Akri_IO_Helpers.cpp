// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_IO_Helpers.hpp>

#include <stk_io/IossBridge.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <stk_util/util/tokenize.hpp>
#include <Ioss_ElementTopology.h>

#include <Akri_DiagWriter.hpp>

namespace krino{

Block_Surface_Connectivity::Block_Surface_Connectivity(const stk::mesh::MetaData & meta)
{
  /* %TRACE[ON]% */
  Trace trace__("Block_Surface_Connectivity::Block_Surface_Connectivity(const stk::mesh::MetaData & meta)"); /* %TRACE% */

  std::vector<const stk::mesh::Part *> surfacesInMap = meta.get_surfaces_in_surface_to_block_map();
  for(auto && surfaceInMap : surfacesInMap)
  {
    std::set<stk::mesh::PartOrdinal> touching_block_ordinals;
    for (auto && touching_part : meta.get_blocks_touching_surface(surfaceInMap))
    {
      touching_block_ordinals.insert(touching_part->mesh_meta_data_ordinal());
    }
    add_surface(surfaceInMap->mesh_meta_data_ordinal(), touching_block_ordinals);
  }
}

void Block_Surface_Connectivity::dump_surface_connectivity(const stk::mesh::MetaData & meta)
{
  const std::vector<const stk::mesh::Part *> surfacesInMap = meta.get_surfaces_in_surface_to_block_map();
  for(auto && surface : surfacesInMap)
  {
    krinolog << "Surface " << surface->name() << " touches blocks ";
    for (auto && touchingBlock : meta.get_blocks_touching_surface(surface))
      krinolog << touchingBlock->name() << " ";
    krinolog << stk::diag::dendl;
  }
}

std::set<stk::mesh::PartOrdinal> Block_Surface_Connectivity::get_surfaces_touching_block(const stk::mesh::PartOrdinal & blockOrdinal) const
{
  auto it = block_to_surface_map.find(blockOrdinal);
  if(it != block_to_surface_map.end())
    return it->second;

  std::set<stk::mesh::PartOrdinal> emptySurfaces;
  return emptySurfaces;
}

std::set<stk::mesh::PartOrdinal> Block_Surface_Connectivity::get_blocks_touching_surface(const stk::mesh::PartOrdinal & surfaceOrdinal) const
{
  auto it = surface_to_block_map.find(surfaceOrdinal);
  if(it != surface_to_block_map.end())
    return it->second;

  std::set<stk::mesh::PartOrdinal> emptyBlocks;
  return emptyBlocks;
}

void Block_Surface_Connectivity::add_surface(const stk::mesh::PartOrdinal & surf_ordinal, const std::set<stk::mesh::PartOrdinal> touching_blocks)
{
  surface_to_block_map[surf_ordinal].insert(touching_blocks.begin(), touching_blocks.end());
  for(auto && block : touching_blocks)
  {
    block_to_surface_map[block].insert(surf_ordinal);
  }
}

stk::mesh::Part * get_parent_sideset(const stk::mesh::Part & part)
{
  stk::mesh::Part * result = nullptr;
  if(part.primary_entity_rank() == part.mesh_meta_data().side_rank())
  {
    for (auto superset : part.supersets())
    {
      if (stk::io::is_part_io_part(superset) && !stk::io::is_part_assembly_io_part(*superset))
      {
        STK_ThrowRequireMsg(result == nullptr, "krino::get_parent_sideset: Part " << part.name() << " has more than 1 parent IO part");
        result = superset;
      }
    }
  }
  return result;
}

std::string get_ioss_topology_name(const stk::topology topology)
{
  const Ioss::ElementTopology * iossTopology = Ioss::ElementTopology::factory(stk::io::map_stk_topology_to_ioss(topology));
  STK_ThrowRequireMsg(iossTopology != nullptr, "krino::get_ioss_topology_name: IOSS topology factory must return a topology type for " << topology.name());
  return iossTopology->name();
}

bool is_sideblock_split_by_topology(const std::string & sideblockName, const stk::topology blockTopology)
{
  // Ideally this would not need string comparisons.
  return sideblockName.find("_" + get_ioss_topology_name(blockTopology) + "_") != std::string::npos;
}

std::string generate_ioss_sideblock_name(const std::string & baseName,
    const std::string & blockOrTopologyName, const std::string & sideTopologyName)
{
  std::vector<std::string> tokens;
  stk::util::tokenize(baseName, "_", tokens);
  const bool isCanonicallyNamed = tokens.size() == 2 && tokens.back().find_first_not_of("0123456789") == std::string::npos;
  if (isCanonicallyNamed)
    return tokens[0] + "_" + blockOrTopologyName + "_" + sideTopologyName + "_" + tokens[1];
  return baseName + "_" + blockOrTopologyName + "_" + sideTopologyName;
}

std::set<unsigned>
get_blocks_and_touching_surfaces(const stk::mesh::MetaData & mesh_meta,
  const stk::mesh::PartVector& input_blocks,
  const Block_Surface_Connectivity & input_block_surface_info)
{
  // For each block that is input_blocks, insert into the set that is returned.
  // Also add surfaces that touch input_blocks to the set.
  std::set<unsigned> blocks_and_touching_sides;
  for (auto && block_ptr : input_blocks)
  {
    blocks_and_touching_sides.insert(block_ptr->mesh_meta_data_ordinal());
    const std::set<stk::mesh::PartOrdinal> touching_surface_ordinals = input_block_surface_info.get_surfaces_touching_block(block_ptr->mesh_meta_data_ordinal());
    for (auto && surf_ordinal : touching_surface_ordinals)
    {
      stk::mesh::Part & surf_part = mesh_meta.get_part(surf_ordinal);
      blocks_and_touching_sides.insert(surf_part.mesh_meta_data_ordinal());

      // Add all subset parts to the set of parts
      if (!(surf_part.subsets().empty()))
      {
        for (stk::mesh::PartVector::const_iterator subset = surf_part.subsets().begin(); subset != surf_part.subsets().end(); ++subset)
        {
          blocks_and_touching_sides.insert((*subset)->mesh_meta_data_ordinal());
        }
      }
    }
  }
  return blocks_and_touching_sides;
}

 // namespace krino
}
