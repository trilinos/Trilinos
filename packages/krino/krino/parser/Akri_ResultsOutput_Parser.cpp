// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_ResultsOutput_Parser.hpp>

#include <Akri_DiagWriter.hpp>
#include <Akri_Region.hpp>
#include <Akri_ResultsOutputOptions.hpp>
#include <Akri_YAML_Parser.hpp>

#include <stk_util/environment/RuntimeDoomed.hpp>

namespace krino {

void
ResultsOutput_Parser::parse(const YAML::Node & region_node, Region & region)
{
  const YAML::Node results_node = YAML_Parser::get_map_if_present(region_node, "output");
  if ( results_node )
  {
    ResultsOutputOptions * options = region.get_results_options();

    std::string results_database_name;
    if (YAML_Parser::get_if_present(results_node, "database_name", results_database_name))
    {
      options->set_filename(results_database_name);
    }
    else
    {
      stk::RuntimeDoomedAdHoc() << "Missing results database_name for region " << region.name() << ".\n";
    }

    std::string results_title;
    if (YAML_Parser::get_if_present(results_node, "title", results_title))
    {
      options->set_title(results_title);
    }

    bool compose_results;
    if (YAML_Parser::get_if_present(results_node, "compose_results", compose_results))
    {
      if (compose_results) options->add_property(Ioss::Property("COMPOSE_RESULTS", 1));
    }

    int output_frequency = 1;
    if (YAML_Parser::get_if_present(results_node, "output_frequency", output_frequency))
    {
      options->add_step_increment(0, output_frequency);
    }

    const YAML::Node nodal_variable_nodes = YAML_Parser::get_sequence_if_present(results_node, "nodal_variables");
    for ( auto && variable_node : nodal_variable_nodes )
    {
      const std::string field_name = variable_node.as<std::string>();
      options->add_nodal_field(field_name, field_name);
    }

    const YAML::Node element_variable_nodes = YAML_Parser::get_sequence_if_present(results_node, "element_variables");
    for ( auto && variable_node : element_variable_nodes )
    {
      const std::string field_name = variable_node.as<std::string>();
      options->add_element_field(field_name, field_name);
    }
  }
}

} // namespace krino
