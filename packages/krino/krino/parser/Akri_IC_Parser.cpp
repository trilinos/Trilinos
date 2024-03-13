// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_IC_Parser.hpp>

#include <Akri_AnalyticSurf.hpp>
#include <Akri_AuxMetaData.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_IC_Alg.hpp>
#include <Akri_LevelSet.hpp>
#include <Akri_Parser.hpp>

#include <stk_util/environment/RuntimeDoomed.hpp>
#include "Akri_Surface_Parser.hpp"

namespace krino {

namespace {

std::unique_ptr<IC_Binder>
parse_binder(const Parser::Node & ic_node)
{
  std::string binder_type;
  if (!ic_node.get_if_present("type", binder_type))
  {
    stk::RuntimeDoomedAdHoc() << "Missing type for Binder IC.\n";
  }

  double interface_size = 0.0;
  double smooth_bridge_size = 0.0;
  double smooth_bridge_offset = 0.0;
  double other_ls_scale_factor = 0.0;
  int ibinder_type = 0;
  bool root_smooth_bridge = false;

  if (binder_type == "interface")
  {
    ibinder_type = 0;
    if (!ic_node.get_if_present("interface_size", interface_size))
    {
      stk::RuntimeDoomedAdHoc() << "Missing interface_size for IC binder.\n";
    }
  }
  else if (binder_type == "smooth_bridge")
  {
    ibinder_type = 1;
    ic_node.get_if_present("smooth_bridge_size", smooth_bridge_size);
    ic_node.get_if_present("smooth_bridge_offset", smooth_bridge_offset);
    ic_node.get_if_present("other_ls_scale_factor", other_ls_scale_factor);
    ic_node.get_if_present("root_smooth_bridge", root_smooth_bridge);

    if (smooth_bridge_size < 0.0)
    {
      stk::RuntimeDoomedAdHoc() << "IC binder: Smooth bridge size should not be negative.\n";
    }

    if (smooth_bridge_size != 0.0)
    {
      if (smooth_bridge_offset <= 0.0)
      {
        stk::RuntimeDoomedAdHoc() << "IC binder: Smooth bridge offset should be greater than zero when a bridge size is specified.\n";
      }
      if (other_ls_scale_factor <= 0.0)
      {
        stk::RuntimeDoomedAdHoc() << "IC binder: Scaling of other level sets should be greater than zero when a bridge size is specified.  Typical values are O(1e2).\n";
      }
    }
  }
  else
  {
    stk::RuntimeDoomedAdHoc() << "Binder IC type should be either interface or smooth_bridge. \n";
  }

  return std::make_unique<IC_Binder>(interface_size, smooth_bridge_size, smooth_bridge_offset, other_ls_scale_factor, ibinder_type, root_smooth_bridge);
}

void
parse_composition_method(const Parser::Node & ic_node, IC_Alg& ic_alg)
{
  // This is a little strange because composition looks like an IC, but really sets a flag
  std::string composition_method;
  STK_ThrowRequire(ic_node.get_if_present("composition_method", composition_method));

  std::transform(composition_method.begin(), composition_method.end(), composition_method.begin(), ::toupper);

  if (composition_method == "MINIMUM_SIGNED_DISTANCE")
    ic_alg.set_composition_method(Composite_Surface::MINIMUM_SIGNED_DISTANCE);
  else if (composition_method == "MAXIMUM_SIGNED_DISTANCE")
    ic_alg.set_composition_method(Composite_Surface::MAXIMUM_SIGNED_DISTANCE);
  else
    stk::RuntimeDoomedAdHoc() << "Unrecognized composition_method: " << composition_method;
}

void
parse_IC(const Parser::Node & ic_node, LevelSet &ls)
{
  Surface * parsedSurface = Surface_Parser::parse(ic_node, ls.meta(), ls.get_timer());
  if (nullptr != parsedSurface)
  {
    ls.get_IC_alg().addSurface(parsedSurface);
  }
  else if ( ic_node.get_null_if_present("binder") )
  {
    ls.get_IC_alg().addCalculator(parse_binder(ic_node));
  }
  else if ( ic_node.get_scalar_if_present("composition_method") )
  {
    parse_composition_method(ic_node, ls.get_IC_alg());
  }
  else
  {
    stk::RuntimeDoomedAdHoc() << "Unrecognized Levelset IC type: " << ic_node.info();
  }
}

}

void
IC_Parser::parse(const Parser::Node & node, LevelSet & ls)
{
  const Parser::Node ic_nodes = node.get_sequence_if_present("initial_conditions");
  if ( ic_nodes )
  {
    for ( auto && ic_node : ic_nodes )
    {
      parse_IC(ic_node, ls);
    }
  }
}

} // namespace krino
