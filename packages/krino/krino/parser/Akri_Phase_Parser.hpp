// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_Phase_Parser_h
#define Akri_Phase_Parser_h

#include <string>

namespace krino { namespace Parser { class Node; } }

namespace krino {
namespace Phase_Parser {
  void parse(const Parser::Node & fem_node, const std::string & fem_model_name);
}
}

#endif // Akri_Phase_Parser_h
