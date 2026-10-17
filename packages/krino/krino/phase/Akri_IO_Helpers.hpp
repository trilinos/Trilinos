// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_IO_Helpers_h
#define Akri_IO_Helpers_h

#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_topology/topology.hpp>

#include <map>
#include <set>
#include <string>
#include <vector>

namespace stk { namespace mesh { class MetaData; } }

namespace krino {

class Block_Surface_Connectivity
{
public:
  Block_Surface_Connectivity() {}
  Block_Surface_Connectivity(const stk::mesh::MetaData & meta);

  std::set<stk::mesh::PartOrdinal> get_surfaces_touching_block(const stk::mesh::PartOrdinal & blockOrdinal) const;
  std::set<stk::mesh::PartOrdinal> get_blocks_touching_surface(const stk::mesh::PartOrdinal & surfaceOrdinal) const;
  void add_surface(const stk::mesh::PartOrdinal & surf_ordinal, const std::set<stk::mesh::PartOrdinal> touching_blocks);

  static void dump_surface_connectivity(const stk::mesh::MetaData & meta);

private:
  std::map< stk::mesh::PartOrdinal, std::set<stk::mesh::PartOrdinal> > block_to_surface_map;
  std::map< stk::mesh::PartOrdinal, std::set<stk::mesh::PartOrdinal> > surface_to_block_map;
};

// The non-assembly IO sideset that owns a side part (e.g. a sideblock subset), or nullptr if there is none.
stk::mesh::Part * get_parent_sideset(const stk::mesh::Part & part);

// Ioss name of an stk topology (e.g. "tri3" for TRIANGLE_3).  Throws if Ioss does not support it.
std::string get_ioss_topology_name(const stk::topology topology);

// True if the sideblock name shows that it was split by element topology rather than by element block.
bool is_sideblock_split_by_topology(const std::string & sideblockName, const stk::topology blockTopology);

// Ioss sideblock name for a sideset base name, the block or topology name it is split by, and its side
// topology name.  For example, ("surface_1", "block_1", "tri3") gives "surface_block_1_tri3_1".
std::string generate_ioss_sideblock_name(const std::string & baseName,
    const std::string & blockOrTopologyName, const std::string & sideTopologyName);

std::set<unsigned> get_blocks_and_touching_surfaces(const stk::mesh::MetaData & mesh_meta,
  const stk::mesh::PartVector& input_blocks,
  const Block_Surface_Connectivity & input_block_surface_info);

} // namespace krino

#endif // Akri_IO_Helpers_h
