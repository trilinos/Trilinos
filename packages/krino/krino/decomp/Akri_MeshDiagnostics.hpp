// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_MeshDiagnostics_h
#define Akri_MeshDiagnostics_h

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Selector.hpp>

namespace krino {

void print_volume_or_surface_area(const stk::mesh::BulkData& mesh, const stk::mesh::EntityRank entity_rank, const stk::mesh::Selector & active_selector, const stk::mesh::PartVector & parts);
double compute_volume_or_surface_area(const stk::mesh::BulkData& mesh, const stk::mesh::EntityRank entity_rank, const stk::mesh::Selector & selector);

} // namespace krino

#endif // Akri_MeshDiagnostics_h
