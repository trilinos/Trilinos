// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_Phase_Parser.hpp>

#include <Akri_DiagWriter.hpp>
#include <Akri_LevelSet.hpp>
#include <Akri_Phase_Support.hpp>
#include <Akri_YAML_Parser.hpp>

#include <stk_util/environment/RuntimeDoomed.hpp>

namespace krino {

void
Phase_Parser::parse(const YAML::Node & fem_node, const std::string & fem_model_name)
{
  PhaseVec & mesh_phases = Phase_Support::get_phases(fem_model_name);

  const YAML::Node subdom_nodes = YAML_Parser::get_sequence_if_present(fem_node, "subdomains");

  if (subdom_nodes)
  {
    for ( auto && subdom_node : subdom_nodes )
    {
      std::string subdom_name;
      YAML_Parser::get_if_present(subdom_node, "name", subdom_name);
      if (subdom_name.empty())
      {
        stk::RuntimeDoomedAdHoc() << "Missing subdomain name in finite element model " << fem_model_name;
      }
      mesh_phases.push_back(NamedPhase(subdom_name));

      const YAML::Node where_nodes = YAML_Parser::get_sequence_if_present(subdom_node, "level_set_regions");

      if (where_nodes)
      {
        for ( auto && where_node : where_nodes )
        {
          std::string negative_ls_name;
          if (YAML_Parser::get_if_present(where_node, "negative_level_set", negative_ls_name))
          {
            mesh_phases.back().tag().add(LevelSet::get_identifier(negative_ls_name),-1);
          }

          std::string positive_ls_name;
          if (YAML_Parser::get_if_present(where_node, "positive_level_set", positive_ls_name))
          {
            mesh_phases.back().tag().add(LevelSet::get_identifier(positive_ls_name),+1);
          }
        }
      }

      std::string smallest_ls_name;
      if (YAML_Parser::get_if_present(subdom_node, "smallest_level_set", smallest_ls_name))
      {
        if ( where_nodes )
        {
          stk::RuntimeDoomedAdHoc() << "Cannot combine \"negative_level_set|positive_level_set\" syntax with \"smallest_level_set\" syntax.";
        }
        Phase_Support::set_one_levelset_per_phase(true);
        mesh_phases.back().tag().add(LevelSet::get_identifier(smallest_ls_name),-1);
      }
    }
  }

  const YAML::Node part_nodes = YAML_Parser::get_sequence_if_present(fem_node, "parts");

  if (part_nodes)
  {
    std::map<std::string,std::vector<std::string>> & FEmodel_block_phase_names = Phase_Support::get_block_phases_by_name(fem_model_name);

    for ( auto && part_node : part_nodes )
    {
      std::string part_name;
      YAML_Parser::get_if_present(part_node, "name", part_name);
      if (part_name.empty())
      {
        stk::RuntimeDoomedAdHoc() << "Missing part name in finite element model " << fem_model_name;
      }

      const YAML::Node part_subdom_nodes = YAML_Parser::get_sequence_if_present(part_node, "subdomains");
      if (part_subdom_nodes)
      {
        std::vector<std::string> & block_phase_names = FEmodel_block_phase_names[part_name];
        for ( auto && subdom_node : part_subdom_nodes )
        {
          const std::string subdom_name = subdom_node.as<std::string>();
          block_phase_names.push_back(subdom_name);
        }
      }
    }
  }
}

} // namespace krino
