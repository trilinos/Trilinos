// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_LevelSet_Parser_h
#define Akri_LevelSet_Parser_h

namespace krino { namespace Parser { class Node; } }
namespace stk { namespace diag { class Timer; } }
namespace stk { namespace mesh { class MetaData; } }

#include <string>

namespace krino {
// Name of a level set or bounding surface from the name key of its input: either given directly
//   name: LS01
// or composed of a prefix and an index zero-padded to a number of digits, given directly or as the
// number of digits of the largest index
//   name: {prefix: LS, index: 1, digits: 2}       -> LS01
//   name: {prefix: LS, index: 1, max_index: 39}   -> LS01
// Returns an empty name if name is absent or invalid (dooming if invalid).
std::string parse_level_set_or_bounding_surface_name(const Parser::Node & surfaceNode);

namespace LevelSet_Parser {
  void parse(const Parser::Node & region_node, stk::mesh::MetaData & meta, const stk::diag::Timer & parentTimer);
}

namespace BoundingSurface_Parser {
  void parse(const Parser::Node & region_node, stk::mesh::MetaData & meta, const stk::diag::Timer & parentTimer);
}
}

#endif // Akri_LevelSet_Parser_h
