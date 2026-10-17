// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.
#include <stk_mesh/base/MetaData.hpp>

#include <Akri_PhasePartInfo.hpp>
#include <Akri_PhaseUtils.hpp>

#include <Akri_PhaseTag.hpp>
#include <Akri_Surface_Identifier.hpp>

namespace krino {

static int get_converted_interface_part_ordinal(const PhasePartInfo & phasePartInfo, const unsigned convertedTouchingPartOrd, const unsigned convertedOppositePartOrd, const int valueToIndicateRemoval)
{
  if (convertedTouchingPartOrd == convertedOppositePartOrd)
    return valueToIndicateRemoval;

  return phasePartInfo.get_interface_part(convertedTouchingPartOrd, convertedOppositePartOrd);
}

std::map<int,int> build_phase_conversion_map_for_element_and_side_parts(
    const stk::mesh::MetaData & meta,
    const PhasePartInfo & phasePartInfo,
    const PhaseTag & fromPhase,
    const PhaseTag & toPhase)
{
  const int valueToIndicateRemoval = -1; // negative value used to indicate that we will remove this interface part
  std::map<int,int> partOrdinalMapping;
  for (auto * part : meta.get_mesh_parts())
  {
    const unsigned partOrd = part->mesh_meta_data_ordinal();
    if ((part->primary_entity_rank() == stk::topology::ELEMENT_RANK || part->primary_entity_rank() == meta.side_rank()))
    {
      if (phasePartInfo.is_conforming_part(partOrd))
      {
        if (fromPhase == phasePartInfo.get_conforming_part_phase(partOrd))
          partOrdinalMapping[partOrd] = phasePartInfo.get_conforming_part(partOrd, toPhase);
      }
      else if (phasePartInfo.is_interface_part(partOrd))
      {
        const unsigned touchingPartOrd = phasePartInfo.get_interface_part_touching_part(partOrd);
        const unsigned oppositePartOrd = phasePartInfo.get_interface_part_opposite_part(partOrd);
        const unsigned toPhasePartOrd = phasePartInfo.get_conforming_part(phasePartInfo.get_nonconforming_part(partOrd), toPhase);

        int convertedInterfacePartOrd = valueToIndicateRemoval;
        bool interfacePartConverted = false;
        if (fromPhase == phasePartInfo.get_conforming_part_phase(touchingPartOrd))
        {
          convertedInterfacePartOrd = get_converted_interface_part_ordinal(phasePartInfo, toPhasePartOrd, oppositePartOrd, valueToIndicateRemoval);
          interfacePartConverted = true;
        }
        else if (fromPhase == phasePartInfo.get_conforming_part_phase(oppositePartOrd))
        {
          convertedInterfacePartOrd = get_converted_interface_part_ordinal(phasePartInfo, touchingPartOrd, toPhasePartOrd, valueToIndicateRemoval);
          interfacePartConverted = true;
        }

        if (interfacePartConverted)
        {
          partOrdinalMapping[partOrd] = convertedInterfacePartOrd;
          if (valueToIndicateRemoval == convertedInterfacePartOrd)
          {
            // Also need to remove interface superset part
            for(const auto * superset : part->supersets())
              if (phasePartInfo.is_interface_superset_part(superset->mesh_meta_data_ordinal()))
                partOrdinalMapping[superset->mesh_meta_data_ordinal()] = valueToIndicateRemoval;
          }
        }
      }
    }
  }
  return partOrdinalMapping;
}

static std::vector<unsigned>
get_indices_of_phases_that_explicity_depend_on_surface(const PhaseVec & meshPhases, const Surface_Identifier & surfaceID)
{
  // gather up the set of phases that depend on this surface
  PhaseTag negPhase;
  PhaseTag posPhase;
  negPhase.add(surfaceID,-1);
  posPhase.add(surfaceID,+1);

  std::vector<unsigned> phases;
  for(unsigned phaseIndex=0; phaseIndex<meshPhases.size(); ++phaseIndex)
  {
    const NamedPhase & phase = meshPhases[phaseIndex];
    if (phase.tag().contain(negPhase) || phase.tag().contain(posPhase))
      phases.push_back(phaseIndex);
  }

  return phases;
}

static std::vector<unsigned> build_unsigned_vector_from_0_to_N(const size_t N)
{
  std::vector<unsigned> v;
  v.reserve(N);
  for (unsigned i=0; i<N; ++i)
    v.push_back(i);
  return v;
}

std::vector<unsigned>
get_indices_of_phases_that_depend_on_surface(const bool oneLSPerPhase, const PhaseVec & meshPhases, const Surface_Identifier & surfaceID)
{
  if (oneLSPerPhase)
    return build_unsigned_vector_from_0_to_N(meshPhases.size()); // Every level set impacts every phase
  return get_indices_of_phases_that_explicity_depend_on_surface(meshPhases, surfaceID);
}

}

