// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_IO_Helpers.hpp>

#include <Ioss_Utils.h>
#include <Ioss_Region.h>
#include <Ioss_ElementBlock.h>
#include <Ioss_SideSet.h>
#include <Ioss_SideBlock.h>
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

Block_Surface_Connectivity::Block_Surface_Connectivity(const stk::mesh::MetaData & meta, const Ioss::Region & io_region)
{
  /* %TRACE[ON]% */
  Trace trace__("Block_Surface_Connectivity::Block_Surface_Connectivity(const Ioss::Region & reg)"); /* %TRACE% */

  std::vector<std::string> side_block_names;
  std::vector<stk::mesh::PartOrdinal> side_block_ordinals;

  for(auto sideset : io_region.get_sidesets())
  {
    side_block_names.clear();
    sideset->block_membership(side_block_names);
    side_block_ordinals.clear();
    for (auto && block_name : side_block_names)
    {
      const stk::mesh::Part * side_block_part = meta.get_part(block_name);
      ThrowRequire(nullptr != side_block_part);
      side_block_ordinals.push_back(side_block_part->mesh_meta_data_ordinal());
    }
    const stk::mesh::Part * side_part = meta.get_part(sideset->name());
    ThrowRequire(nullptr != side_part);
    add_surface(side_part->mesh_meta_data_ordinal(), std::set<stk::mesh::PartOrdinal>(side_block_ordinals.begin(), side_block_ordinals.end()));

    if (!sideset->get_side_blocks().empty())
    {
      for (auto&& side_subset : sideset->get_side_blocks())
      {
        // Fmwk only creates subset if more than 1 sideblock, but stk always creates them, so just check.
        const stk::mesh::Part * side_subset_part = meta.get_part(side_subset->name());
        if (nullptr == side_subset_part) continue;
        side_block_names.clear();
        side_subset->block_membership(side_block_names);
        side_block_ordinals.clear();
        for (auto && block_name : side_block_names)
        {
          const stk::mesh::Part * side_block_part = meta.get_part(block_name);
          ThrowRequire(nullptr != side_block_part);
          side_block_ordinals.push_back(side_block_part->mesh_meta_data_ordinal());
        }
        add_surface(side_subset_part->mesh_meta_data_ordinal(), std::set<stk::mesh::PartOrdinal>(side_block_ordinals.begin(), side_block_ordinals.end()));
      }
    }
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
