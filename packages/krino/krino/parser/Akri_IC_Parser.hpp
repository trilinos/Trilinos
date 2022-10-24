// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_IC_Parser_h
#define Akri_IC_Parser_h

namespace krino { class LevelSet; }
namespace krino { namespace Parser { class Node; } }

namespace krino {
namespace IC_Parser {
  void parse(const Parser::Node & node, LevelSet & ls);
}
}

#endif // Akri_IC_Parser_h
