// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_LevelSet_Parser.hpp>

#include <Akri_AuxMetaData.hpp>
#include <Akri_IC_Parser.hpp>
#include <Akri_IO_Helpers.hpp>
#include <Akri_LevelSet.hpp>
#include <Akri_Phase_Support.hpp>
#include <stk_math/StkVector.hpp>
#include <Akri_Parser.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_util/environment/RuntimeDoomed.hpp>
#include <stk_util/environment/RuntimeWarning.hpp>
#include <Akri_BoundingSurface.hpp>
#include <Akri_CSG_Parser.hpp>
#include <Akri_CSG_Tree.hpp>
#include <Akri_IC_Alg.hpp>
#include <Akri_Surface_Parser.hpp>
#include <Akri_Surface_Manager.hpp>

namespace krino {

namespace {

void
register_blocks_for_level_set(Surface_Manager & surfaceManager, Phase_Support & phaseSupport, LevelSet & ls)
{
  const std::string composite_name = ls.get_composite_name();
  if (!composite_name.empty())
  {
    LS_SideTag::declare_composite(ls.get_identifier(), surfaceManager.get_identifier(composite_name));
  }
  phaseSupport.register_blocks_for_surface(ls.get_identifier());
}

void
register_blocks_for_bounding_surface(Surface_Manager & surfaceManager, Phase_Support & phaseSupport, const std::string & boundingSurfaceName)
{
  phaseSupport.register_blocks_for_surface(surfaceManager.get_identifier(boundingSurfaceName));
}

size_t num_decimal_digits(int value)
{
  return std::to_string(value).size();
}

}

std::string parse_level_set_or_bounding_surface_name(const Parser::Node & surfaceNode)
{
  const Parser::Node nameNode = surfaceNode.get_if_present("name");
  if (!nameNode)
    return std::string();

  if (!nameNode.is_map())
  {
    std::string name;
    surfaceNode.get_if_present("name", name);
    return name;
  }

  const auto oldDoomedCount = stk::get_doomed_count();

  std::string prefix;
  if (!nameNode.get_if_present("prefix", prefix) || prefix.empty())
    stk::RuntimeDoomedAdHoc() << "Missing prefix for a name given by prefix and index.\n";

  int index = -1;
  if (!nameNode.get_if_present("index", index))
    stk::RuntimeDoomedAdHoc() << "Missing index for a name given by prefix and index.\n";
  else if (index < 0)
    stk::RuntimeDoomedAdHoc() << "The index " << index << " of the name with prefix " << prefix << " must not be negative.\n";

  int digits = 0;
  int maxIndex = -1;
  const bool hasDigits = nameNode.get_if_present("digits", digits);
  const bool hasMaxIndex = nameNode.get_if_present("max_index", maxIndex);
  size_t width = 0;
  if (hasDigits == hasMaxIndex)
  {
    stk::RuntimeDoomedAdHoc() << "A name given by prefix and index needs exactly one of digits or max_index.\n";
  }
  else if (hasDigits)
  {
    if (digits < 1)
      stk::RuntimeDoomedAdHoc() << "The digits " << digits << " of the name with prefix " << prefix << " must be positive.\n";
    else if (index >= 0 && num_decimal_digits(index) > size_t(digits))
      stk::RuntimeDoomedAdHoc() << "The index " << index << " of the name with prefix " << prefix << " has more than " << digits << " digits.\n";
    width = digits;
  }
  else
  {
    if (index >= 0 && index > maxIndex)
      stk::RuntimeDoomedAdHoc() << "The index " << index << " of the name with prefix " << prefix << " is greater than its max_index " << maxIndex << ".\n";
    width = (maxIndex >= 0) ? num_decimal_digits(maxIndex) : 0;
  }

  if (stk::get_doomed_count() > oldDoomedCount)
    return std::string();

  const std::string indexString = std::to_string(index);
  return prefix + std::string(width - indexString.size(), '0') + indexString;
}

void
LevelSet_Parser::parse(const Parser::Node & region_node, stk::mesh::MetaData & meta, const stk::diag::Timer & parentTimer)
{
  const Parser::Node ls_nodes = region_node.get_sequence_if_present("level_set_interfaces");
  if ( ls_nodes )
  {
    Surface_Manager & surfaceManager = Surface_Manager::get(meta);
    Phase_Support & phaseSupport = Phase_Support::get(meta);
    for ( auto && ls_node : ls_nodes )
    {
      const std::string ls_name = parse_level_set_or_bounding_surface_name(ls_node);
      if (ls_name.empty())
      {
        stk::RuntimeDoomedAdHoc() << "Blank or missing levelset name.\n";
      }
      LevelSet & ls = LevelSet::build(meta, ls_name, parentTimer);

      std::string distance_name;
      if (ls_node.get_if_present("distance_variable", distance_name))
      {
        ls.set_levelset_field_name(distance_name);
      }

      std::vector<std::string> interfaceVelocity;
      if (ls_node.get_if_present("interface_velocity", interfaceVelocity))
      {
        if (interfaceVelocity.size() != ls.spatial_dimension)
        {
          stk::RuntimeDoomedAdHoc() << "Expecting " << ls.spatial_dimension << " string expressions for interface_velocity for level set " << ls_name << ".\n";
        }
        ls.set_interface_velocity(interfaceVelocity);
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
        else if (redistance_method_name == "FAST_ITERATIVE")
          redistance_method = FAST_ITERATIVE;
        else
          stk::RuntimeWarningAdHoc() << "Unrecognized redistance method:  " << redistance_method_name << std::endl;

        ls.set_redistance_method(redistance_method);
      }

      std::string semilagrangianMethodName;
      if (ls_node.get_if_present("semilagrangian_algorithm", semilagrangianMethodName))
      {
        std::transform(semilagrangianMethodName.begin(), semilagrangianMethodName.end(), semilagrangianMethodName.begin(), ::toupper);
        SemiLagrangianAlgorithm algType = NON_ADAPTIVE_SINGLE_STEP;
        if (semilagrangianMethodName == "NON_ADAPTIVE_SINGLE_STEP")
          algType = NON_ADAPTIVE_SINGLE_STEP;
        else if (semilagrangianMethodName == "ADAPTIVE_SINGLE_STEP")
                  algType = ADAPTIVE_SINGLE_STEP;
        else if (semilagrangianMethodName == "ADAPTIVE_PREDICTOR_CORRECTOR")
          algType = ADAPTIVE_PREDICTOR_CORRECTOR;
        else
          stk::RuntimeWarningAdHoc() << "Unrecognized redistance method:  " << redistance_method_name << std::endl;

        ls.set_semilagrangian_algorithm(algType);
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

          STK_ThrowErrorMsgIf( !ls.aux_meta().has_part(surface_name),
              "Could not locate a surface named " << surface_name);

          stk::mesh::Part & io_part = ls.aux_meta().get_part(surface_name);
          STK_ThrowErrorMsgIf( ls.meta().side_rank() != io_part.primary_entity_rank(),
            "Part " << surface_name << " is not a side-rank part.");

          ls.get_compute_surface_distance_parts().push_back(&io_part);
          compute_distance_surfaces.push_back(surface_name);
        }
        STK_ThrowErrorMsgIf( ls.get_compute_surface_distance_parts().empty(),
          "Please specify surfaces for compute surface distance.");
      }

      CSGTree csgTree = CSG_Parser::parse_CSG_tree(ls_node, ls.meta(), ls.get_timer());
      if (!csgTree.empty())
      {
        ls.get_IC_alg().set_csg(csgTree);
      }
      else
      {
        IC_Parser::parse(ls_node, ls);
        ls.get_IC_alg().set_simple_union();
      }

      register_blocks_for_level_set(surfaceManager, phaseSupport, ls);
    }
  }
}

void
BoundingSurface_Parser::parse(const Parser::Node & region_node, stk::mesh::MetaData & meta, const stk::diag::Timer & parentTimer)
{
  const Parser::Node boundingSurfaceNodes = region_node.get_sequence_if_present("bounding_surfaces");
  if ( boundingSurfaceNodes )
  {
    Surface_Manager & surfaceManager = Surface_Manager::get(meta);
    Phase_Support & phaseSupport = Phase_Support::get(meta);
    for ( auto && boundingSurfaceNode : boundingSurfaceNodes )
    {
      const std::string boundingSurfaceName = parse_level_set_or_bounding_surface_name(boundingSurfaceNode);
      if (boundingSurfaceName.empty())
      {
        stk::RuntimeDoomedAdHoc() << "Blank or missing bounding surface name.\n";
      }

      std::shared_ptr<Surface> surface = Surface_Parser::parse(boundingSurfaceNode, meta, parentTimer);
      BoundingSurface::build(meta, boundingSurfaceName, surface);

      register_blocks_for_bounding_surface(surfaceManager, phaseSupport, boundingSurfaceName);
    }
  }
}

} // namespace krino
