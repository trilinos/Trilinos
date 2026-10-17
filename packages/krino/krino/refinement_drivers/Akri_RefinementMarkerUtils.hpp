// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_REFINEMENT_DRIVERS_AKRI_REFINEMENTMARKERUTILS_HPP_
#define KRINO_KRINO_REFINEMENT_DRIVERS_AKRI_REFINEMENTMARKERUTILS_HPP_

#include <Akri_FieldRef.hpp>
#include <Akri_RefinementDrivers.hpp>
#include <stk_mesh/base/Types.hpp>

#include <vector>

namespace krino {

class RefinementManager;

// Element marker helpers for RefinementManager's marker field, used with perform_multilevel_adaptivity.

// Fills the marker field with COARSEN if isDefaultCoarsen, else NOTHING.
void initialize_marker(const stk::mesh::BulkData& mesh,
    const RefinementManager & refinement,
    const bool isDefaultCoarsen);

// Marker of an element: REFINE if it is selected for refinement and below minRefineLevel, NOTHING if it is
// selected and at or above it, else COARSEN if isDefaultCoarsen, else NOTHING.
Refinement_Marker determine_refinement_marker(const bool isElementSelectedForRefinement,
    const int minRefineLevel,
    const int elementRefineLevel,
    const bool isDefaultCoarsen);

// Sets the marker of each of the given elements by determine_refinement_marker, as selected for refinement.
void mark_elements_to_reach_refinement_level(const stk::mesh::BulkData& mesh,
    const RefinementManager & refinement,
    const std::vector<stk::mesh::Entity> elementsToMark,
    const int minRefineLevel,
    const bool isDefaultCoarsen);

// Whether to keep coarsening by default: only for the first targetRefineLevel refinement iterations, to avoid
// cycling when the elements selected for refinement change with refinement.
bool should_continue_to_coarsen(const int refinementIterCount, const int targetRefineLevel);

// Logs the minimum and maximum edge length of the given (non-transition) elements at each refinement level.
void write_refinement_level_sizes(const stk::mesh::BulkData& mesh,
    const RefinementManager & refinement,
    const FieldRef coordsField,
    const std::vector<stk::mesh::Entity> & elements,
    const int maxRefineLevel);

} // namespace krino

#endif /* KRINO_KRINO_REFINEMENT_DRIVERS_AKRI_REFINEMENTMARKERUTILS_HPP_ */
