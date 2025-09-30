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

#include <stk_util/environment/RuntimeDoomed.hpp>
#include "Akri_Parser.hpp"

namespace krino {

void
ResultsOutput_Parser::parse(const Parser::Node & region_node, Region & region)
{
  const Parser::Node results_node = region_node.get_map_if_present("output");
  if ( results_node )
  {
    ResultsOutputOptions * options = region.get_results_options();

    std::string results_database_name;
    if (results_node.get_if_present("database_name", results_database_name))
    {
      options->set_filename(results_database_name);
    }
    else
    {
      stk::RuntimeDoomedAdHoc() << "Missing results database_name for region " << region.name() << ".\n";
    }

    std::string results_title;
    if (results_node.get_if_present("title", results_title))
    {
      options->set_title(results_title);
    }

    bool compose_results{false};
    if (results_node.get_if_present("compose_results", compose_results))
    {
      if (compose_results) options->add_property(Ioss::Property("COMPOSE_RESULTS", 1));
    }

    int output_frequency = 1;
    if (results_node.get_if_present("output_frequency", output_frequency))
    {
      options->add_step_increment(0, output_frequency);
    }

    const Parser::Node nodal_variable_nodes = results_node.get_sequence_if_present("nodal_variables");
    for ( auto && variable_node : nodal_variable_nodes )
    {
      const std::string field_name = variable_node.as<std::string>();
      options->add_nodal_field(field_name, field_name);
    }

    const Parser::Node element_variable_nodes = results_node.get_sequence_if_present("element_variables");
    for ( auto && variable_node : element_variable_nodes )
    {
      const std::string field_name = variable_node.as<std::string>();
      options->add_element_field(field_name, field_name);
    }
  }
}


} // namespace krino
