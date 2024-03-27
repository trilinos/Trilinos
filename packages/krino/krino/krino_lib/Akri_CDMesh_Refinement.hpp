// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_INCLUDE_AKRI_CDMESH_REFINEMENT_H_
#define KRINO_INCLUDE_AKRI_CDMESH_REFINEMENT_H_

#include <stk_mesh/base/BulkData.hpp>

namespace krino {

class InterfaceGeometry;
class CDFEM_Support;
class CDFEM_Snapper;
class RefinementSupport;
class Phase_Support;

void
mark_possible_cut_elements_for_adaptivity(const stk::mesh::BulkData& mesh,
      const RefinementInterface & refinement,
      const InterfaceGeometry & interfaceGeometry,
      const RefinementSupport & refinementSupport,
      const int numRefinements);

void
mark_elements_that_intersect_interval(const stk::mesh::BulkData& mesh,
      const RefinementInterface & refinement,
      const InterfaceGeometry & interfaceGeometry,
      const std::array<double,2> refinementInterval,
      const int numRefineLevels,
      const bool isDefaultCoarsen);

void
mark_interface_elements_for_adaptivity(const stk::mesh::BulkData& mesh,
      const RefinementInterface & refinement,
      const InterfaceGeometry & interfaceGeometry,
      const RefinementSupport & refinementSupport,
      const FieldRef coords_field,
      const int num_refinements);

std::vector<std::pair<stk::mesh::Entity,unsigned>> get_owned_adaptivity_parents_and_their_element_part(const stk::mesh::BulkData& mesh,
    const RefinementInterface & refinement,
    const Phase_Support & phaseSupport);

}

#endif /* KRINO_INCLUDE_AKRI_CDMESH_REFINEMENT_H_ */
