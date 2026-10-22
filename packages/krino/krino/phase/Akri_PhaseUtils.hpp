// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_PHASE_AKRI_PHASEUTILS_HPP_
#define KRINO_KRINO_PHASE_AKRI_PHASEUTILS_HPP_
#include <Akri_PhaseTag.hpp>

#include <map>
#include <vector>

namespace stk { namespace mesh { class MetaData; } }

namespace krino {

class PhasePartInfo;
class Surface_Identifier;

std::map<int,int> build_phase_conversion_map_for_element_and_side_parts(
    const stk::mesh::MetaData & meta,
    const PhasePartInfo & phasePartInfo,
    const PhaseTag & fromPhase,
    const PhaseTag & toPhase);

std::vector<unsigned> get_indices_of_phases_that_depend_on_surface(
    const bool oneLSPerPhase,
    const PhaseVec & meshPhases,
    const Surface_Identifier & surfaceID);

}

#endif /* KRINO_KRINO_PHASE_AKRI_PHASEUTILS_HPP_ */
