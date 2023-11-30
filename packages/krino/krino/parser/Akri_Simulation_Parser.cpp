// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_Simulation_Parser.hpp>

#include <Akri_DiagWriter.hpp>
#include <Akri_Simulation.hpp>
#include <Akri_Region_Parser.hpp>
#include <Akri_Parser.hpp>

#include <stk_util/environment/RuntimeDoomed.hpp>

namespace krino {

void
Simulation_Parser::parse(Simulation & simulation, const Parser::Node & node)
{
  const Parser::Node sim_node = node.get_map_if_present("simulation");
  if ( sim_node )
  {
    double start_time = 0.0;
    if (sim_node.get_if_present("start_time", start_time))
    {
      simulation.set_current_time(start_time);
    }

    double stop_time = 0.0;
    if (sim_node.get_if_present("stop_time", stop_time))
    {
      simulation.set_stop_time(stop_time);
    }

    double time_step = 0.0;
    if (sim_node.get_if_present("time_step", time_step))
    {
      simulation.set_time_step(time_step);
    }

    Region_Parser::parse(sim_node, simulation);
  }
  else
  {
    stk::RuntimeDoomedAdHoc() << "Missing simulation.\n";
  }
}

} // namespace krino
