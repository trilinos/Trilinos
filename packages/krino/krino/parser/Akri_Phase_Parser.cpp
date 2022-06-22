// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_Phase_Parser.hpp>

#include <Akri_DiagWriter.hpp>
#include <Akri_Phase_Support.hpp>

#include <stk_util/environment/RuntimeDoomed.hpp>
#include <Akri_Surface_Manager.hpp>
#include <Akri_Parser.hpp>

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

void
Phase_Parser::parse(const Parser::Node & fem_node, const std::string & fem_model_name)
{
  Phase_Support & phaseSupport = Phase_Support::get_or_create(fem_model_name);
  PhaseVec & meshPhases = phaseSupport.get_mesh_phases();
  Surface_Manager & surfaceManager = Surface_Manager::get_or_create(fem_model_name);

  const Parser::Node subdom_nodes = fem_node.get_sequence_if_present("subdomains");

  if (subdom_nodes)
  {
    for ( auto && subdom_node : subdom_nodes )
    {
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
            meshPhases.back().tag().add(surfaceManager.get_identifier(negative_ls_name),-1);
          }

          std::string positive_ls_name;
          if (where_node.get_if_present("positive_level_set", positive_ls_name))
          {
            meshPhases.back().tag().add(surfaceManager.get_identifier(positive_ls_name),+1);
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
        phaseSupport.set_one_levelset_per_phase(true);
        meshPhases.back().tag().add(surfaceManager.get_identifier(smallest_ls_name),-1);
      }

      std::string boundingSurfaceName;
      if (subdom_node.get_if_present("bounding_surface", boundingSurfaceName))
      {
        if ( !meshPhases.back().tag().empty() )
        {
          stk::RuntimeDoomedAdHoc() << "Cannot combine \"negative_level_set|positive_level_set|smallest_level_set\" syntax with bounding_surface syntax.";
        }
        meshPhases.back().tag().add(surfaceManager.get_identifier(boundingSurfaceName),-1);
      }
    }
  }

  set_elsewhere_phase(meshPhases);

  const Parser::Node part_nodes = fem_node.get_sequence_if_present("parts");

  if (part_nodes)
  {
    auto & modelBlockPhaseNames = phaseSupport.get_block_phases_by_name();

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
}

} // namespace krino
