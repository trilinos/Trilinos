// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_Phase_Parser_h
#define Akri_Phase_Parser_h

#include <Akri_PhaseInput.hpp>

#include <string>
#include <vector>

namespace krino { namespace Parser { class Node; } }

namespace krino {
namespace Phase_Parser {
  // Parses the subdomains and parts of a finite element model.  Level set and bounding surface names
  // are given surface identifiers in the returned input's surfaceNames.
  PhaseInput parse(const Parser::Node & fem_node, const std::string & fem_model_name);

  // Generates the subdomains of each between_consecutive_level_sets and per_smallest_level_set entry of
  // phaseInput from the given level set and bounding surface names, inserting them where the entry was,
  // then sets the elsewhere phase.  The names matching an entry's level_set_prefix (case-insensitive) are
  // ordered case-insensitively, so numbered names must be padded (e.g. LS00 ... LS10).  Subdomains are
  // named subdomain_prefix followed by a level set's name without level_set_prefix: for
  // between_consecutive_level_sets, the upper level set of each consecutive pair; for
  // per_smallest_level_set, each level set.  The prefixes are given either nested under the keyword or
  // beside it in the same subdomain entry.
  void expand_generated_subdomains(PhaseInput & phaseInput, const std::vector<std::string> & surfaceNames);
}
}

#endif // Akri_Phase_Parser_h
