// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_Surface_Parser_h
#define Akri_Surface_Parser_h
#include <functional>

namespace stk { namespace diag { class Timer; } }
namespace stk { namespace mesh { class MetaData; } }
namespace krino { namespace Parser { class Node; } }
namespace krino { class Surface; }

namespace krino {
namespace Surface_Parser {
  Surface * parse(const Parser::Node & parserNode, const stk::mesh::MetaData & meta, const stk::diag::Timer &parentTimer);
}
}

#endif // Akri_Surface_Parser_h
