// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_INCLUDE_AKRI_CDMESH_UTILS_H_
#define KRINO_INCLUDE_AKRI_CDMESH_UTILS_H_
#include <stk_mesh/base/BulkData.hpp>

namespace krino {

class CDFEM_Support;
class PhaseTag;
class InterfaceID;
class AuxMetaData;
class Phase_Support;

bool parts_are_compatible_for_snapping(const stk::mesh::BulkData & mesh, stk::mesh::Entity possible_snap_node, stk::mesh::Entity fixed_node);
bool parts_are_compatible_for_snapping_when_ignoring_phase(const stk::mesh::BulkData & mesh, const AuxMetaData & auxMeta, const Phase_Support & phaseSupport, stk::mesh::Entity possible_snap_node, stk::mesh::Entity fixed_node);
bool phase_matches_interface(const CDFEM_Support & cdfemSupport, const PhaseTag & phase, const InterfaceID interface);
bool determine_phase_from_parts(PhaseTag & phase, const stk::mesh::PartVector & parts, const Phase_Support & phaseSupport);
PhaseTag determine_phase_for_entity(const stk::mesh::BulkData & mesh, stk::mesh::Entity entity, const Phase_Support & phaseSupport);
bool nodes_are_on_any_interface(const stk::mesh::BulkData & mesh, const Phase_Support & phaseSupport, const stk::mesh::Bucket & nodeBucket);
bool node_is_on_any_interface(const stk::mesh::BulkData & mesh, const Phase_Support & phaseSupport, const stk::mesh::Entity node);
bool node_is_on_interface(const stk::mesh::BulkData & mesh, const CDFEM_Support & cdfemSupport, const Phase_Support & phaseSupport, stk::mesh::Entity node, const InterfaceID & interface);

}



#endif /* KRINO_INCLUDE_AKRI_CDMESH_UTILS_H_ */
