// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_Region_Parser.hpp>

#include <Akri_CDFEM_Options_Parser.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_IO_Helpers.hpp>
#include <Akri_LevelSet_Parser.hpp>
#include <Akri_Phase_Support.hpp>
#include <Akri_Region.hpp>
#include <Akri_RegionInterface.hpp>
#include <Akri_ResultsOutput_Parser.hpp>
#include "Akri_Parser.hpp"

#include <stk_util/environment/RuntimeDoomed.hpp>

namespace krino {

void
Region_Parser::parse(const Parser::Node & simulation_node, Simulation & simulation)
{
  const Parser::Node region_nodes = simulation_node.get_sequence_if_present("regions");
  if ( !region_nodes ) return;

  for ( auto && region_node : region_nodes )
  {
    std::string region_name;
    region_node.get_if_present("name", region_name);
    if (region_name.empty())
    {
      stk::RuntimeDoomedAdHoc() << "Blank or missing region name.\n";
    }
    Region * region = new Region(simulation, region_name );

    int initial_refinement_levels = 0;
    if (region_node.get_if_present("initial_uniform_refinement_levels", initial_refinement_levels))
    {
      region->set_initial_refinement_levels(initial_refinement_levels);
    }

    bool use_32bit_ids = false;
    region_node.get_if_present("use_32bit_ids", use_32bit_ids);

    bool force_64bit_ids = !use_32bit_ids;
    region_node.get_if_present("force_64bit_ids", force_64bit_ids);

    if (use_32bit_ids && force_64bit_ids)
    {
      stk::RuntimeDoomedAdHoc() << "Can't specify options use_32bit_ids=true and force_64bit_ids=true together in region " << region_name << "\n";
    }

    std::string fem_model_name;
    region_node.get_if_present("finite_element_model", fem_model_name);
    if (fem_model_name.empty())
    {
      stk::RuntimeDoomedAdHoc() << "Blank or missing finite element model for region " << region_name << "\n";
    }
    region->associate_input_mesh(fem_model_name, use_32bit_ids, force_64bit_ids);

    RegionInterface::set_currently_parsed_region(region->get_stk_mesh_meta_data(), region->name(), region->getRegionTimer(), region->name_of_input_mesh(), region->get_input_io_region());

    Phase_Support & phase_support = krino::Phase_Support::get(region->get_stk_mesh_meta_data());
    const PhaseVec & mesh_phases = Phase_Support::get_phases(region->name_of_input_mesh());
    const PartnamePhasenameMap & mesh_block_phase_names = Phase_Support::get_block_phases_by_name(region->name_of_input_mesh());
    phase_support.determine_block_phases(mesh_phases, mesh_block_phase_names);


    const Block_Surface_Connectivity block_surf_info(region->get_stk_mesh_meta_data());
    phase_support.set_input_block_surface_connectivity(block_surf_info);
    phase_support.setup_phases(mesh_phases);

    LevelSet_Parser::parse(region_node, RegionInterface::get_currently_parsed_region());
    CDFEM_Options_Parser::parse(region_node, RegionInterface::get_currently_parsed_region());
    ResultsOutput_Parser::parse(region_node, *region);
  }
}

} // namespace krino
