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

namespace krino {
namespace LevelSet_Parser {
  void parse(const Parser::Node & region_node, stk::mesh::MetaData & meta, const stk::diag::Timer & parentTimer);
}

namespace BoundingSurface_Parser {
  void parse(const Parser::Node & region_node, stk::mesh::MetaData & meta, const stk::diag::Timer & parentTimer);
}
}

#endif // Akri_LevelSet_Parser_h
