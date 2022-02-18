// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_IO_Helpers_h
#define Akri_IO_Helpers_h
//
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Selector.hpp>

#include <Akri_PhaseTag.hpp>

#include <map>

namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class MetaData; } }
namespace stk { namespace diag { class Timer; } }
namespace Ioss { class Region; }

namespace krino {

class AuxMetaData;

class Block_Surface_Connectivity
{
public:
  Block_Surface_Connectivity() {}
  Block_Surface_Connectivity(const stk::mesh::MetaData & meta);
  Block_Surface_Connectivity(const stk::mesh::MetaData & meta, const Ioss::Region & io_region);
  void get_surfaces_touching_block(const stk::mesh::PartOrdinal & block_ordinal,
      std::set<stk::mesh::PartOrdinal> & surface_ordinals) const
  {
    auto it = block_to_surface_map.find(block_ordinal);
    if(it == block_to_surface_map.end()) return;
    surface_ordinals.insert(it->second.begin(), it->second.end());
  }
  void get_blocks_touching_surface(const stk::mesh::PartOrdinal & surface_ordinal,
      std::set<stk::mesh::PartOrdinal> & block_ordinals) const
  {
    block_ordinals.clear();
    auto it = surface_to_block_map.find(surface_ordinal);
    if(it == surface_to_block_map.end()) return;
    block_ordinals.insert(it->second.begin(), it->second.end());
  }

  void add_surface(const stk::mesh::PartOrdinal & surf_ordinal, const std::set<stk::mesh::PartOrdinal> touching_blocks)
  {
    surface_to_block_map[surf_ordinal].insert(touching_blocks.begin(), touching_blocks.end());
    for(auto && block : touching_blocks)
    {
      block_to_surface_map[block].insert(surf_ordinal);
    }
  }
private:
  std::map< stk::mesh::PartOrdinal, std::set<stk::mesh::PartOrdinal> > block_to_surface_map;
  std::map< stk::mesh::PartOrdinal, std::set<stk::mesh::PartOrdinal> > surface_to_block_map;
};

std::vector<std::string> get_input_mesh_block_names( const Ioss::Region & io_region );
} // namespace krino

#endif // Akri_IO_Helpers_h
