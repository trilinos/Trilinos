// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_PHASE_AKRI_PHASEINPUT_HPP_
#define KRINO_KRINO_PHASE_AKRI_PHASEINPUT_HPP_

#include <Akri_PhaseTag.hpp>
#include <Akri_Surface_Identifier.hpp>

#include <map>
#include <string>
#include <vector>

namespace krino {

typedef std::map<std::string,std::vector<std::string>> PartnamePhasenameMap;

// Surface names in identifier order, so that Surface_Identifier(i) names surfaceNames[i].
class SurfaceNames
{
public:
  // Identifier of the named surface (case-insensitive), adding it if it is new.
  Surface_Identifier get_identifier(const std::string & name);
  const std::vector<std::string> & get_names() const { return myNames; }
private:
  std::vector<std::string> myNames;
};

// The phases parsed from a finite element model's input, before its mesh exists.  Parsers build this,
// and it is given to each mesh created from that model: Surface_Manager::set_surface_names(surfaceNames)
// first, so the identifiers in meshPhases are valid there, then Phase_Support::set_phase_input.
// A subdomain entry that generates subdomains from the level sets (or bounding surfaces) whose names
// start with levelSetPrefix.  It is expanded once the level set names are known
// (Phase_Parser::expand_generated_subdomains), with the generated subdomains inserted at phaseIndex in
// meshPhases.
struct GeneratedSubdomains
{
  enum Type
  {
    BETWEEN_CONSECUTIVE_LEVEL_SETS, // one subdomain between each pair of consecutive level sets
    PER_SMALLEST_LEVEL_SET          // one subdomain where each level set is the smallest
  };
  Type type;
  std::string subdomainPrefix;
  std::string levelSetPrefix;
  size_t phaseIndex;
};

struct PhaseInput
{
  SurfaceNames surfaceNames;
  PhaseVec meshPhases;
  PartnamePhasenameMap blockPhasesByName;
  bool oneLevelSetPerPhase = false;
  std::vector<GeneratedSubdomains> generatedSubdomains; // not yet expanded
};

} // namespace krino

#endif /* KRINO_KRINO_PHASE_AKRI_PHASEINPUT_HPP_ */
