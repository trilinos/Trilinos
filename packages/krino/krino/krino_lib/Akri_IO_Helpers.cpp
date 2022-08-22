// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_IO_Helpers.hpp>

#include <Ioss_Region.h>
#include <Ioss_ElementBlock.h>
#include <stk_io/IossBridge.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>

#include <Akri_AuxMetaData.hpp>
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

std::vector<std::string>
get_input_mesh_block_names( const Ioss::Region & io_region )
{
  std::vector<std::string> mesh_elem_blocks;
  for(auto && elem_block : io_region.get_element_blocks())
  {
    if (stk::io::include_entity(elem_block))
    {
      mesh_elem_blocks.push_back(elem_block->name());
    }
  }
  return mesh_elem_blocks;
}
 // namespace krino
}
