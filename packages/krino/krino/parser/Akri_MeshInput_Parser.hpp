// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_MeshInput_Parser_h
#define Akri_MeshInput_Parser_h

namespace krino { class MeshInputOptions; }
namespace krino { namespace Parser { class Node; } }

namespace krino {
namespace MeshInput_Parser {
 void parse(const Parser::Node & node);
 bool parse_generated_mesh(const Parser::Node & fem_node, MeshInputOptions & options);
}
}

#endif // Akri_MeshInput_Parser_h
