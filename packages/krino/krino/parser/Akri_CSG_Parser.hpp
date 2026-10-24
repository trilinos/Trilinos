// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_PARSER_AKRI_CSG_PARSER_HPP_
#define KRINO_KRINO_PARSER_AKRI_CSG_PARSER_HPP_

#include <memory>
#include <string>

#include <stk_mesh/base/Types.hpp>

namespace stk { namespace diag { class Timer; } }

namespace krino { namespace Parser { class Node; } }
namespace krino { class CSGTree; }

namespace krino {
namespace CSG_Parser {
  CSGTree parse_CSG_tree(const Parser::Node & node, const stk::mesh::MetaData & meta, const stk::diag::Timer &lsTimer);
}
}

#endif /* KRINO_KRINO_PARSER_AKRI_CSG_PARSER_HPP_ */
