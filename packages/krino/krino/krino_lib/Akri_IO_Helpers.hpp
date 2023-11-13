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
#include <map>

namespace stk { namespace mesh { class MetaData; } }
namespace Ioss { class Region; }

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

std::vector<std::string> get_input_mesh_block_names( const Ioss::Region & io_region );
} // namespace krino

#endif // Akri_IO_Helpers_h
