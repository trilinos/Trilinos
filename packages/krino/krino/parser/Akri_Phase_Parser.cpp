// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_Phase_Parser.hpp>

#include <stk_util/environment/RuntimeDoomed.hpp>
#include <Akri_Parser.hpp>

#include <algorithm>
#include <cctype>
#include <set>
#include <string>
#include <vector>

namespace krino {

void set_elsewhere_phase(PhaseVec & meshPhases)
{
  NamedPhase * elsewherePhase = nullptr;
  for (auto && meshPhase : meshPhases)
  {
    if (meshPhase.tag().empty())
    {
      if (nullptr != elsewherePhase)
        stk::RuntimeDoomedAdHoc() << "More than one phase with no levelset or bounding surface specification.";
      elsewherePhase = &meshPhase;
    }
  }

  if (nullptr != elsewherePhase)
  {
    for (auto && meshPhase : meshPhases)
      if (&meshPhase != elsewherePhase)
        elsewherePhase->tag().add_opposite_sides(meshPhase.tag());
  }
}

static const char * generated_subdomains_keyword(const GeneratedSubdomains::Type type)
{
  return (type == GeneratedSubdomains::BETWEEN_CONSECUTIVE_LEVEL_SETS) ? "between_consecutive_level_sets" : "per_smallest_level_set";
}

static void parse_generated_subdomains(const Parser::Node & subdomNode, const Parser::Node & generatedNode, const GeneratedSubdomains::Type type, PhaseInput & phaseInput)
{
  const auto keyword = generated_subdomains_keyword(type);
  for (const auto otherKey : {"name", "level_set_regions", "smallest_level_set", "bounding_surface"})
    if (subdomNode.get_if_present(otherKey))
      stk::RuntimeDoomedAdHoc() << "Cannot combine " << keyword << " with " << otherKey << " in a subdomain.\n";

  GeneratedSubdomains generated;
  generated.type = type;
  generatedNode.get_if_present("subdomain_prefix", generated.subdomainPrefix);
  generatedNode.get_if_present("level_set_prefix", generated.levelSetPrefix);
  if (generated.subdomainPrefix.empty())
    stk::RuntimeDoomedAdHoc() << "Missing subdomain_prefix for " << keyword << ".\n";
  if (generated.levelSetPrefix.empty())
    stk::RuntimeDoomedAdHoc() << "Missing level_set_prefix for " << keyword << ".\n";

  // Like smallest_level_set, which these subdomains are, this makes the model one level set per phase.
  if (type == GeneratedSubdomains::PER_SMALLEST_LEVEL_SET)
    phaseInput.oneLevelSetPerPhase = true;

  generated.phaseIndex = phaseInput.meshPhases.size();
  phaseInput.generatedSubdomains.push_back(generated);
}

static std::string to_upper(std::string str)
{
  std::transform(str.begin(), str.end(), str.begin(), [](const unsigned char c) { return std::toupper(c); });
  return str;
}

static bool starts_with_case_insensitive(const std::string & name, const std::string & prefix)
{
  return name.size() >= prefix.size() && to_upper(name.substr(0, prefix.size())) == to_upper(prefix);
}

static std::vector<std::string> get_sorted_matching_names(const std::vector<std::string> & surfaceNames, const std::string & prefix)
{
  std::vector<std::string> matchingNames;
  for (const std::string & name : surfaceNames)
    if (starts_with_case_insensitive(name, prefix))
      matchingNames.push_back(name);
  // Surface names are case-insensitive, so order them that way.  Numbered names must be padded
  // (LS00 ... LS10) to order correctly.
  std::sort(matchingNames.begin(), matchingNames.end(),
      [](const std::string & a, const std::string & b) { return to_upper(a) < to_upper(b); });
  return matchingNames;
}

static std::string generated_subdomain_name(const GeneratedSubdomains & generated, const std::string & levelSetName)
{
  return generated.subdomainPrefix + levelSetName.substr(generated.levelSetPrefix.size());
}

static PhaseVec generate_subdomains_between_consecutive_level_sets(const GeneratedSubdomains & generated, const std::vector<std::string> & matchingNames, SurfaceNames & phaseSurfaceNames)
{
  PhaseVec subdomains;
  for (size_t k=1; k<matchingNames.size(); ++k)
  {
    const std::string & lowerName = matchingNames[k-1];
    const std::string & upperName = matchingNames[k];
    NamedPhase subdomain(generated_subdomain_name(generated, upperName));
    subdomain.tag().add(phaseSurfaceNames.get_identifier(lowerName), +1);
    subdomain.tag().add(phaseSurfaceNames.get_identifier(upperName), -1);
    subdomains.push_back(subdomain);
  }
  return subdomains;
}

static PhaseVec generate_subdomains_per_smallest_level_set(const GeneratedSubdomains & generated, const std::vector<std::string> & matchingNames, SurfaceNames & phaseSurfaceNames)
{
  PhaseVec subdomains;
  for (const std::string & levelSetName : matchingNames)
  {
    NamedPhase subdomain(generated_subdomain_name(generated, levelSetName));
    subdomain.tag().add(phaseSurfaceNames.get_identifier(levelSetName), -1);
    subdomains.push_back(subdomain);
  }
  return subdomains;
}

static PhaseVec generate_subdomains(const GeneratedSubdomains & generated, const std::vector<std::string> & surfaceNames, SurfaceNames & phaseSurfaceNames)
{
  const std::vector<std::string> matchingNames = get_sorted_matching_names(surfaceNames, generated.levelSetPrefix);
  const size_t minNumMatches = (generated.type == GeneratedSubdomains::BETWEEN_CONSECUTIVE_LEVEL_SETS) ? 2 : 1;
  if (matchingNames.size() < minNumMatches)
  {
    stk::RuntimeDoomedAdHoc() << generated_subdomains_keyword(generated.type) << " with level_set_prefix " << generated.levelSetPrefix
        << " matched " << matchingNames.size() << " level sets or bounding surfaces; at least " << minNumMatches << " needed.\n";
    return {};
  }

  if (generated.type == GeneratedSubdomains::BETWEEN_CONSECUTIVE_LEVEL_SETS)
    return generate_subdomains_between_consecutive_level_sets(generated, matchingNames, phaseSurfaceNames);
  return generate_subdomains_per_smallest_level_set(generated, matchingNames, phaseSurfaceNames);
}

static void check_for_duplicate_subdomain_names(const PhaseVec & meshPhases)
{
  std::set<std::string> names;
  for (const NamedPhase & phase : meshPhases)
    if (!names.insert(to_upper(phase.name())).second)
      stk::RuntimeDoomedAdHoc() << "Subdomain " << phase.name() << " is defined more than once.\n";
}

void
Phase_Parser::expand_generated_subdomains(PhaseInput & phaseInput, const std::vector<std::string> & surfaceNames)
{
  if (phaseInput.generatedSubdomains.empty())
    return;

  // Insert from the last entry, so the phaseIndex of earlier entries stays valid.
  for (auto iter = phaseInput.generatedSubdomains.rbegin(); iter != phaseInput.generatedSubdomains.rend(); ++iter)
  {
    const PhaseVec subdomains = generate_subdomains(*iter, surfaceNames, phaseInput.surfaceNames);
    phaseInput.meshPhases.insert(phaseInput.meshPhases.begin() + iter->phaseIndex, subdomains.begin(), subdomains.end());
  }
  phaseInput.generatedSubdomains.clear();

  check_for_duplicate_subdomain_names(phaseInput.meshPhases);
  set_elsewhere_phase(phaseInput.meshPhases);
}

PhaseInput
Phase_Parser::parse(const Parser::Node & fem_node, const std::string & fem_model_name)
{
  PhaseInput phaseInput;
  PhaseVec & meshPhases = phaseInput.meshPhases;
  SurfaceNames & surfaceNames = phaseInput.surfaceNames;

  const Parser::Node subdom_nodes = fem_node.get_sequence_if_present("subdomains");

  if (subdom_nodes)
  {
    for ( auto && subdom_node : subdom_nodes )
    {
      if (const Parser::Node betweenNode = subdom_node.get_options_if_present("between_consecutive_level_sets"))
      {
        parse_generated_subdomains(subdom_node, betweenNode, GeneratedSubdomains::BETWEEN_CONSECUTIVE_LEVEL_SETS, phaseInput);
        continue;
      }
      if (const Parser::Node perSmallestNode = subdom_node.get_options_if_present("per_smallest_level_set"))
      {
        parse_generated_subdomains(subdom_node, perSmallestNode, GeneratedSubdomains::PER_SMALLEST_LEVEL_SET, phaseInput);
        continue;
      }

      std::string subdom_name;
      subdom_node.get_if_present("name", subdom_name);
      if (subdom_name.empty())
      {
        stk::RuntimeDoomedAdHoc() << "Missing subdomain name in finite element model " << fem_model_name;
      }
      meshPhases.push_back(NamedPhase(subdom_name));

      const Parser::Node where_nodes = subdom_node.get_sequence_if_present("level_set_regions");

      if (where_nodes)
      {
        for ( auto && where_node : where_nodes )
        {
          std::string negative_ls_name;
          if (where_node.get_if_present("negative_level_set", negative_ls_name))
          {
            meshPhases.back().tag().add(surfaceNames.get_identifier(negative_ls_name),-1);
          }

          std::string positive_ls_name;
          if (where_node.get_if_present("positive_level_set", positive_ls_name))
          {
            meshPhases.back().tag().add(surfaceNames.get_identifier(positive_ls_name),+1);
          }
        }
      }

      std::string smallest_ls_name;
      if (subdom_node.get_if_present("smallest_level_set", smallest_ls_name))
      {
        if ( where_nodes )
        {
          stk::RuntimeDoomedAdHoc() << "Cannot combine \"negative_level_set|positive_level_set\" syntax with \"smallest_level_set\" syntax.";
        }
        phaseInput.oneLevelSetPerPhase = true;
        meshPhases.back().tag().add(surfaceNames.get_identifier(smallest_ls_name),-1);
      }

      std::string boundingSurfaceName;
      if (subdom_node.get_if_present("bounding_surface", boundingSurfaceName))
      {
        if ( !meshPhases.back().tag().empty() )
        {
          stk::RuntimeDoomedAdHoc() << "Cannot combine \"negative_level_set|positive_level_set|smallest_level_set\" syntax with bounding_surface syntax.";
        }
        meshPhases.back().tag().add(surfaceNames.get_identifier(boundingSurfaceName),-1);
      }
    }
  }

  // With between_consecutive_level_sets or per_smallest_level_set entries the elsewhere phase also
  // excludes their generated subdomains, so it is set when they are expanded.
  if (phaseInput.generatedSubdomains.empty())
    set_elsewhere_phase(meshPhases);

  const Parser::Node part_nodes = fem_node.get_sequence_if_present("parts");

  if (part_nodes)
  {
    auto & modelBlockPhaseNames = phaseInput.blockPhasesByName;

    for ( auto && part_node : part_nodes )
    {
      std::string part_name;
      part_node.get_if_present("name", part_name);
      if (part_name.empty())
      {
        stk::RuntimeDoomedAdHoc() << "Missing part name in finite element model " << fem_model_name;
      }

      const Parser::Node part_subdom_nodes = part_node.get_sequence_if_present("subdomains");
      if (part_subdom_nodes)
      {
        std::vector<std::string> & blockPhaseNames = modelBlockPhaseNames[part_name];
        for ( auto && subdom_node : part_subdom_nodes )
        {
          const std::string subdom_name = subdom_node.as<std::string>();
          blockPhaseNames.push_back(subdom_name);
        }
      }
    }
  }
  return phaseInput;
}

} // namespace krino
