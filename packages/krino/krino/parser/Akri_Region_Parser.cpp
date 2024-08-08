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
#include <Akri_RefinementSupport.hpp>
#include <Akri_Region.hpp>
#include <Akri_ResultsOutput_Parser.hpp>
#include <Akri_Surface_Manager.hpp>
#include "Akri_Parser.hpp"

#include <stk_util/environment/RuntimeDoomed.hpp>

namespace krino {

static void parse_postprocessors(const Parser::Node & regionNode, PostProcessors & postprocessors)
{
  const Parser::Node ppNodes = regionNode.get_sequence_if_present("postprocessors");
  if ( ppNodes )
  {
    for ( auto && ppNode : ppNodes )
    {
      std::string fieldName;
      ppNode.get_if_present("field", fieldName);
      if (fieldName.empty())
      {
        stk::RuntimeDoomedAdHoc() << "Blank or missing field name in postprocessor.\n";
      }
      std::string analyticalExpr;
      ppNode.get_if_present("expression", analyticalExpr);
      if (analyticalExpr.empty())
      {
        stk::RuntimeDoomedAdHoc() << "Blank or missing expression in postprocessor.\n";
      }
      postprocessors.add_scalar_postprocesor(fieldName, analyticalExpr);
    }
  }
}

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

    stk::mesh::MetaData & meta = region->mesh_meta_data();
    RefinementSupport & refinementSupport = RefinementSupport::get(meta);

    int initial_refinement_levels = 0;
    if (region_node.get_if_present("initial_uniform_refinement_levels", initial_refinement_levels))
    {
      refinementSupport.set_initial_refinement_levels(initial_refinement_levels);
    }

    int levelsetRefinementLevels = 0;
    region_node.get_if_present("level_set_refinement_levels", levelsetRefinementLevels);
    if (levelsetRefinementLevels < 0)
    {
      stk::RuntimeDoomedAdHoc() << "Invalid negative level_set_refinement_levels " << levelsetRefinementLevels << "\n";
    }
    refinementSupport.activate_nonconformal_adaptivity(levelsetRefinementLevels);

    std::vector<double> refinementIntervalData;
    if (region_node.get_if_present("level_set_refinement_interval", refinementIntervalData))
    {
      if (refinementIntervalData.size() != 2)
      {
        stk::RuntimeDoomedAdHoc() << "level_set_refinement_interval must be a vector of length 2 (lo,hi).\n";
      }
      const std::array<double,2> refinementInterval{refinementIntervalData[0], refinementIntervalData[1]};
      if (!refinementSupport.is_valid_interval(refinementInterval))
      {
        stk::RuntimeDoomedAdHoc() << "level_set_refinement_interval must be a vector of length 2 (lo,hi) with lo <= hi.\n";
      }
      refinementSupport.set_refinement_interval(refinementInterval);
    }

    parse_postprocessors(region_node, region->get_postprocessors());

    const stk::diag::Timer & regionTimer = region->getRegionTimer();
    Phase_Support::associate_FEModel_and_metadata(fem_model_name, meta);
    Surface_Manager::associate_FEModel_and_metadata(fem_model_name, meta);
    Phase_Support & phase_support = krino::Phase_Support::get(meta);
    phase_support.determine_block_phases();

    const Block_Surface_Connectivity block_surf_info(meta);
    phase_support.set_input_block_surface_connectivity(block_surf_info);
    phase_support.setup_phases();

    LevelSet_Parser::parse(region_node, meta, regionTimer);
    BoundingSurface_Parser::parse(region_node, meta, regionTimer);
    CDFEM_Options_Parser::parse(region_node, meta);
    ResultsOutput_Parser::parse(region_node, *region);
  }
}

} // namespace krino
