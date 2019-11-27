// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef percept_SmoothingHelperFunctions_hpp
#define percept_SmoothingHelperFunctions_hpp

#include <string>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/BoundaryAnalysis.hpp>
#include <stk_mesh/base/BulkModification.hpp>

#include <percept/PerceptMesh.hpp> // FIXME don't include all this!!!

namespace percept
{
  /// return the set of nodes that are on the outer skin or on shared boundaries between
  /// "blocks" which are defined as element_rank() parts
  stk::mesh::Part* get_skin_part(stk::mesh::BulkData * bulkData, const std::string& part_name, bool remove_previous_part_nodes=true);
}

#endif
