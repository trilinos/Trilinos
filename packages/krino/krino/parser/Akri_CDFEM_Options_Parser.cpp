// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_CDFEM_Options_Parser.hpp>

#include <Akri_DiagWriter.hpp>
#include <Akri_CDFEM_Support.hpp>
#include <Akri_Parser.hpp>

#include <stk_util/environment/RuntimeDoomed.hpp>
#include <Akri_RefinementSupport.hpp>

namespace krino {

void
CDFEM_Options_Parser::parse(const Parser::Node & region_node, stk::mesh::MetaData & meta)
{
  const Parser::Node cdfem_node = region_node.get_map_if_present("cdfem_options");
  if ( cdfem_node )
  {
    CDFEM_Support & cdfem_support = CDFEM_Support::get(meta);
    RefinementSupport & refinementSupport = RefinementSupport::get(meta);

    std::string cdfem_edge_degeneracy_handling_string;
    if (cdfem_node.get_if_present("cdfem_edge_degeneracy_handling", cdfem_edge_degeneracy_handling_string))
    {
      std::transform(cdfem_edge_degeneracy_handling_string.begin(), cdfem_edge_degeneracy_handling_string.end(), cdfem_edge_degeneracy_handling_string.begin(), ::toupper);
      static std::map<std::string, Edge_Degeneracy_Handling> valid_entries =
        { {"SNAP_TO_NODE", SNAP_TO_NODE}, {"SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS", SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE} };
      auto it = valid_entries.find(cdfem_edge_degeneracy_handling_string);
      if (it == valid_entries.end())
      {
        stk::RuntimeDoomedAdHoc() << "Invalid cdfem_edge_degeneracy_handling type: " << cdfem_node.info();
      }
      else
      {
        cdfem_support.set_cdfem_edge_degeneracy_handling( it->second );
        if (it->second == SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE)
        {
          const double defaultEdgeTolWhenSnapping = 0.05;
          cdfem_support.set_cdfem_edge_tol( defaultEdgeTolWhenSnapping );
        }
      }
    }

    bool myFlagUseFacetsInsteadOfLevelSets = false;
    cdfem_node.get_if_present("use_facets_to_perform_decomposition", myFlagUseFacetsInsteadOfLevelSets);
    cdfem_support.set_use_facets_instead_of_levelset_fields( myFlagUseFacetsInsteadOfLevelSets );

    double cdfem_edge_tol = 0.0;
    if (cdfem_node.get_if_present("cdfem_edge_tolerance", cdfem_edge_tol))
    {
      const double minimum_edge_tol = 1.e-12;
      if (cdfem_edge_tol < minimum_edge_tol)
      {
        krinolog << "Using minimum edge tolerance of  " << minimum_edge_tol << " instead of specified tolerance of " << cdfem_edge_tol << stk::diag::dendl;
        cdfem_edge_tol = minimum_edge_tol;
      }
      cdfem_support.set_cdfem_edge_tol( cdfem_edge_tol );
    }

    double maxEdgeSnap = 1.0;
    if (cdfem_node.get_if_present("max_edge_snap", maxEdgeSnap))
    {
      if (maxEdgeSnap < 0. || maxEdgeSnap > 1.)
      {
        stk::RuntimeDoomedAdHoc() << "max_edge_snap must be in the range of 0 to 1: " << cdfem_node.info();
      }
      cdfem_support.set_max_edge_snap( maxEdgeSnap );
    }

    double snappingSharpFeatureAngleInDegrees = 135.;
    if (cdfem_node.get_if_present("sharp_feature_angle_when_snapping", snappingSharpFeatureAngleInDegrees))
    {
      cdfem_support.set_snapping_sharp_feature_angle_in_degrees( snappingSharpFeatureAngleInDegrees );
    }

    std::string cdfem_simplex_generation_method_string;
    if (cdfem_node.get_if_present("cdfem_simplex_generation_method", cdfem_simplex_generation_method_string))
    {
      std::transform(cdfem_simplex_generation_method_string.begin(), cdfem_simplex_generation_method_string.end(), cdfem_simplex_generation_method_string.begin(), ::toupper);
      static std::map<std::string, Simplex_Generation_Method> valid_entries = {
        {"CUT_QUADS_BY_GLOBAL_IDENTIFIER", CUT_QUADS_BY_GLOBAL_IDENTIFIER},
        {"CUT_QUADS_BY_LARGEST_ANGLE", CUT_QUADS_BY_LARGEST_ANGLE},
        {"CUT_QUADS_BY_NEAREST_EDGE_CUT", CUT_QUADS_BY_NEAREST_EDGE_CUT}
      };
      auto it = valid_entries.find(cdfem_simplex_generation_method_string);
      if (it == valid_entries.end())
      {
        stk::RuntimeDoomedAdHoc() << "Invalid cdfem_simplex_generation_method type: " << cdfem_node.info();
      }
      else
      {
        cdfem_support.set_simplex_generation_method( it->second );
      }
    }

    int num_init_decomp_cycles = 0;
    if (cdfem_node.get_if_present("number_of_initial_decomposition_cycles", num_init_decomp_cycles))
    {
      cdfem_support.set_num_initial_decomposition_cycles( num_init_decomp_cycles );
    }

    bool interface_refinement_specified = false;
    int interface_minimum_refinement_level = 0;
    if(cdfem_node.get_if_present("cdfem_interface_minimum_refinement_level", interface_minimum_refinement_level))
    {
      interface_refinement_specified = true;
    }

    int interface_maximum_refinement_level = 0;
    if(cdfem_node.get_if_present("cdfem_interface_maximum_refinement_level", interface_maximum_refinement_level))
    {
      interface_refinement_specified = true;
    }

    if (interface_refinement_specified)
    {
      refinementSupport.activate_interface_refinement( interface_minimum_refinement_level, interface_maximum_refinement_level );
    }

    int nonconformal_adapt_levels = 0;
    if (cdfem_node.get_if_present("cdfem_nonconformal_adaptivity_levels", nonconformal_adapt_levels))
    {
      refinementSupport.activate_nonconformal_adaptivity( nonconformal_adapt_levels );
    }

    int post_adapt_refine_levels = 0;
    if (cdfem_node.get_if_present("post_adaptivity_uniform_refinement_levels", post_adapt_refine_levels))
    {
      refinementSupport.set_post_adapt_refinement_levels( post_adapt_refine_levels );
    }

    uint64_t nonconformal_adapt_target_element_count = 0;
    if (cdfem_node.get_if_present("nonconformal_adaptivity_target_element_count", nonconformal_adapt_target_element_count))
    {
      refinementSupport.activate_nonconformal_adapt_target_count( nonconformal_adapt_target_element_count );
    }

    int post_cdfem_refinement_levels = 0;
    if (cdfem_node.get_if_present("post_cdfem_refinement_levels", post_cdfem_refinement_levels))
    {
      cdfem_support.set_post_cdfem_refinement_levels( post_cdfem_refinement_levels );
    }

    const Parser::Node post_cdfem_refinement_blocks_seq = cdfem_node.get_sequence_if_present("post_cdfem_refinement_blocks");
    if (post_cdfem_refinement_blocks_seq)
    {
      std::vector<std::string> post_cdfem_refinement_blocks;
      for ( auto && post_cdfem_refinement_block_node : post_cdfem_refinement_blocks_seq )
      {
        const std::string post_cdfem_refinement_block = post_cdfem_refinement_block_node.as<std::string>();
        post_cdfem_refinement_blocks.push_back(post_cdfem_refinement_block);
      }
      cdfem_support.set_post_cdfem_refinement_blocks( post_cdfem_refinement_blocks );
    }

    bool myFlagDoNearbyRefinementBeforeInterfaceRefinement = false;
    cdfem_node.get_if_present("perform_initial_refinement_near_interfaces", myFlagDoNearbyRefinementBeforeInterfaceRefinement);
    refinementSupport.do_nearby_refinement_before_interface_refinement( myFlagDoNearbyRefinementBeforeInterfaceRefinement );
  }
}

} // namespace krino
