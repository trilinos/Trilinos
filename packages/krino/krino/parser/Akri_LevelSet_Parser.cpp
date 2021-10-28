// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_LevelSet_Parser.hpp>

#include <Akri_AuxMetaData.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_IC_Parser.hpp>
#include <Akri_IO_Helpers.hpp>
#include <Akri_LevelSet.hpp>
#include <Akri_Phase_Support.hpp>
#include <Akri_RegionInterface.hpp>
#include <Akri_Vec.hpp>
#include <Akri_Parser.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_util/environment/RuntimeDoomed.hpp>
#include <stk_util/environment/RuntimeWarning.hpp>

namespace krino {

namespace {

void
register_blocks_for_level_set(RegionInterface & reg, LevelSet & ls)
{
  const std::string composite_name = ls.get_composite_name();
  if (!composite_name.empty())
  {
    LS_SideTag::declare_composite(ls.get_identifier(), LevelSet::get_identifier(composite_name));
  }

  Phase_Support & phase_support = Phase_Support::get(reg.get_stk_mesh_meta_data());
  const PhaseVec & mesh_phases = Phase_Support::get_phases(reg.name_of_input_mesh());
  const std::vector<unsigned> ls_phases = Phase_Support::get_level_set_phases(mesh_phases, ls);
  const std::vector<stk::mesh::Part*> decomposed_blocks = phase_support.get_blocks_decomposed_by_levelset(ls_phases);
  phase_support.register_blocks_for_level_set(&ls, decomposed_blocks);
}

}

void
LevelSet_Parser::parse(const Parser::Node & region_node, RegionInterface & region)
{
  const Parser::Node ls_nodes = region_node.get_sequence_if_present("level_set_interfaces");
  if ( ls_nodes )
  {
    for ( auto && ls_node : ls_nodes )
    {
      std::string ls_name;
      ls_node.get_if_present("name", ls_name);
      if (ls_name.empty())
      {
        stk::RuntimeDoomedAdHoc() << "Blank or missing levelset name.\n";
      }
      LevelSet & ls = LevelSet::build(region.get_stk_mesh_meta_data(), ls_name, region.getRegionTimer());

      std::string distance_name;
      if (ls_node.get_if_present("distance_variable", distance_name))
      {
        ls.set_distance_name(distance_name);
      }

      std::vector<double> extension_velocity;
      if (ls_node.get_if_present("extension_velocity", extension_velocity))
      {
        if (extension_velocity.size() != ls.spatial_dimension)
        {
          stk::RuntimeDoomedAdHoc() << "Expecting " << ls.spatial_dimension << " real values for extension_velocity for level set " << ls_name << ".\n";
        }
        ls.set_extension_velocity(Vector3d(extension_velocity.data(), extension_velocity.size()));
      }

      double narrow_band_multiplier = 0.0;
      if (ls_node.get_if_present("narrow_band_element_size_multiplier", narrow_band_multiplier))
      {
        if( narrow_band_multiplier < 0. )
        {
          stk::RuntimeDoomedAdHoc() << "Error: Narrow band element size multiplier must be >= 0.\n";
        }
        else if ( narrow_band_multiplier < 1. )
        {
          stk::RuntimeWarningAdHoc() << "Narrow band element size multiplier is less than 1.  "
              << "Except in certain cases of adaptive refinement around the interface, this will produce errors in the distance field."
              << std::endl;
        }
        ls.narrow_band_multiplier(narrow_band_multiplier);
      }

      std::string redistance_method_name;
      if (ls_node.get_if_present("redistance_method", redistance_method_name))
      {
        std::transform(redistance_method_name.begin(), redistance_method_name.end(), redistance_method_name.begin(), ::toupper);
        Redistance_Method redistance_method = CLOSEST_POINT;
        if (redistance_method_name == "CLOSEST_POINT")
          redistance_method = CLOSEST_POINT;
        else if (redistance_method_name == "FAST_MARCHING")
          redistance_method = FAST_MARCHING;
        else
          stk::RuntimeWarningAdHoc() << "Unrecognized redistance method:  " << redistance_method_name << std::endl;

        ls.set_redistance_method(redistance_method);
      }

      bool perform_initial_redistance;
      if (ls_node.get_if_present("perform_initial_redistance", perform_initial_redistance))
      {
        ls.perform_initial_redistance(perform_initial_redistance);
      }

      double initial_offset_distance = 0.0;
      if (ls_node.get_if_present("initial_offset_distance", initial_offset_distance))
      {
        ls.set_ic_offset(initial_offset_distance);
      }

      double initial_scale_factor = 1.0;
      if (ls_node.get_if_present("initial_scale_factor", initial_scale_factor))
      {
        ls.set_ic_scale(initial_scale_factor);
      }

      double narrow_band_size = 0.0;
      if (ls_node.get_if_present("narrow_band_size", narrow_band_size))
      {
        if( narrow_band_size < 0. )
        {
          stk::RuntimeDoomedAdHoc() << "Error: Narrow band size must be >= 0.\n";
        }
        ls.narrow_band_size(narrow_band_size);
      }

      const Parser::Node comp_dist_surfs = ls_node.get_sequence_if_present("compute_surface_distance");
      if ( comp_dist_surfs )
      {
        std::vector<std::string> compute_distance_surfaces;
        for ( auto && comp_dist_surf : comp_dist_surfs )
        {
          const std::string surface_name = comp_dist_surf.as<std::string>();

          ThrowErrorMsgIf( !ls.aux_meta().has_part(surface_name),
              "Could not locate a surface named " << surface_name);

          stk::mesh::Part & io_part = ls.aux_meta().get_part(surface_name);
          ThrowErrorMsgIf( ls.meta().side_rank() != io_part.primary_entity_rank(),
            "Part " << surface_name << " is not a side-rank part.");

          ls.get_compute_surface_distance_parts().push_back(&io_part);
          compute_distance_surfaces.push_back(surface_name);
        }
        ThrowErrorMsgIf( ls.get_compute_surface_distance_parts().empty(),
          "Please specify surfaces for compute surface distance.");
      }

      IC_Parser::parse(ls_node, ls);

      register_blocks_for_level_set(region, ls);
    }
  }
}

} // namespace krino
