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
#include <Akri_MeshSurface.hpp>
#include <Akri_YAML_Parser.hpp>

#include <stk_util/environment/RuntimeDoomed.hpp>

namespace krino {

namespace {

Sphere *
parse_sphere(const YAML::Node & ic_node)
{
  std::string name;
  YAML_Parser::get_if_present(ic_node, "name", name);

  std::vector<double> center;
  if (YAML_Parser::get_if_present(ic_node, "center", center))
  {
    if (center.size() != 3)
    {
      stk::RuntimeDoomedAdHoc() << "Expecting 3 real values for center for IC sphere.\n";
    }
  }
  else
  {
    stk::RuntimeDoomedAdHoc() << "Missing center for IC sphere.\n";
  }

  double sign = 1.0;
  if (YAML_Parser::get_if_present(ic_node, "sign", sign))
  {
    if (sign != -1.0 && sign == 1.0)
    {
      stk::RuntimeDoomedAdHoc() << "Sign for sphere must be -1 or 1.\n";
    }
  }

  double radius = 0.0;
  if (!YAML_Parser::get_if_present(ic_node, "radius", radius))
  {
    stk::RuntimeDoomedAdHoc() << "Missing radius for IC sphere.\n";
  }

  return new Sphere(name, Vector3d(center.data()), radius, sign);
}

Ellipsoid *
parse_ellipsoid(const YAML::Node & ic_node)
{
  std::string name;
  YAML_Parser::get_if_present(ic_node, "name", name);

  std::vector<double> center;
  if (YAML_Parser::get_if_present(ic_node, "center", center))
  {
    if (center.size() != 3)
    {
      stk::RuntimeDoomedAdHoc() << "Expecting 3 real values for center for IC sphere.\n";
    }
  }
  else
  {
    stk::RuntimeDoomedAdHoc() << "Missing center for IC sphere.\n";
  }

  std::vector<double> semiaxes;
  if (YAML_Parser::get_if_present(ic_node, "semiaxes", semiaxes))
  {
    if (semiaxes.size() != 3)
    {
      stk::RuntimeDoomedAdHoc() << "Expecting 3 real values for semiaxes for IC ellipsoid.\n";
    }
  }
  else
  {
    stk::RuntimeDoomedAdHoc() << "Missing semiaxes for IC ellipsoid.\n";
  }

  double sign = 1.0;
  if (YAML_Parser::get_if_present(ic_node, "sign", sign))
  {
    if (sign != -1.0 && sign == 1.0)
    {
      stk::RuntimeDoomedAdHoc() << "Sign for ellipsoid must be -1 or 1.\n";
    }
  }

  std::vector<double> rotationVec;
  if (YAML_Parser::get_if_present(ic_node, "rotation", rotationVec))
  {
    if (semiaxes.size() != 3)
    {
      stk::RuntimeDoomedAdHoc() << "Expecting 3 real values for rotation for IC ellipsoid.\n";
    }
  }

  return new Ellipsoid(name, center, semiaxes, rotationVec, sign);
}

Plane *
parse_plane(const YAML::Node & ic_node)
{
  std::string name;
  YAML_Parser::get_if_present(ic_node, "name", name);

  std::vector<double> normal;
  if (YAML_Parser::get_if_present(ic_node, "normal", normal))
  {
    if (normal.size() != 3)
    {
      stk::RuntimeDoomedAdHoc() << "Expecting 3 real values for center for IC plane.\n";
    }
  }
  else
  {
    stk::RuntimeDoomedAdHoc() << "Missing normal for IC plane.\n";
  }

  double multiplier = 1.0;
  YAML_Parser::get_if_present(ic_node, "multiplier", multiplier);

  double sign = 1.0;
  if (YAML_Parser::get_if_present(ic_node, "sign", sign))
  {
    if (sign != -1.0 && sign == 1.0)
    {
      stk::RuntimeDoomedAdHoc() << "Sign for IC plane must be -1 or 1.\n";
    }
  }
  multiplier *= sign;

  double offset = 0.0;
  if (!YAML_Parser::get_if_present(ic_node, "offset", offset))
  {
    stk::RuntimeDoomedAdHoc() << "Missing offset for IC plane.\n";
  }

  return new Plane(name, normal.data(), offset, multiplier);
}

Cylinder *
parse_cylinder(const YAML::Node & ic_node)
{
  std::string name;
  YAML_Parser::get_if_present(ic_node, "name", name);

  std::vector<double> p1;
  if (YAML_Parser::get_if_present(ic_node, "p1", p1))
  {
    if (p1.size() != 3)
    {
      stk::RuntimeDoomedAdHoc() << "Expecting 3 real values for p1 for IC cylinder.\n";
    }
  }
  else
  {
    stk::RuntimeDoomedAdHoc() << "Missing normal for IC cylinder.\n";
  }

  std::vector<double> p2;
  if (YAML_Parser::get_if_present(ic_node, "p2", p2))
  {
    if (p2.size() != 3)
    {
      stk::RuntimeDoomedAdHoc() << "Expecting 3 real values for center for IC cylinder.\n";
    }
  }
  else
  {
    stk::RuntimeDoomedAdHoc() << "Missing p2 for IC cylinder.\n";
  }

  double radius = 0.0;
  if (!YAML_Parser::get_if_present(ic_node, "radius", radius))
  {
    stk::RuntimeDoomedAdHoc() << "Missing radius for IC cylinder.\n";
  }

  double sign = 1.0;
  if (YAML_Parser::get_if_present(ic_node, "sign", sign))
  {
    if (sign != -1.0 && sign == 1.0)
    {
      stk::RuntimeDoomedAdHoc() << "Sign for IC cylinder must be -1 or 1.\n";
    }
  }

  return new Cylinder(name, p1.data(), p2.data(), radius, sign);
}

std::unique_ptr<IC_Binder>
parse_binder(const YAML::Node & ic_node)
{
  std::string binder_type;
  if (!YAML_Parser::get_if_present(ic_node, "type", binder_type))
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
    if (!YAML_Parser::get_if_present(ic_node, "interface_size", interface_size))
    {
      stk::RuntimeDoomedAdHoc() << "Missing interface_size for IC binder.\n";
    }
  }
  else if (binder_type == "smooth_bridge")
  {
    ibinder_type = 1;
    YAML_Parser::get_if_present(ic_node, "smooth_bridge_size", smooth_bridge_size);
    YAML_Parser::get_if_present(ic_node, "smooth_bridge_offset", smooth_bridge_offset);
    YAML_Parser::get_if_present(ic_node, "other_ls_scale_factor", other_ls_scale_factor);
    YAML_Parser::get_if_present(ic_node, "root_smooth_bridge", root_smooth_bridge);

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

Faceted_Surface *
parse_facets(const YAML::Node & ic_node, const stk::diag::Timer &parent_timer)
{
  std::string surface_name;
  YAML_Parser::get_if_present(ic_node, "name", surface_name);

  std::string facet_filename;
  if (!YAML_Parser::get_if_present(ic_node, "filename", facet_filename))
  {
    stk::RuntimeDoomedAdHoc() << "Missing filename for IC facets.\n";
  }

  std::string facet_format;
  if (!YAML_Parser::get_if_present(ic_node, "format", facet_format))
  {
    stk::RuntimeDoomedAdHoc() << "Missing format for IC facets.\n";
  }

  bool scaleSpecified = false;
  double scale = 1.0;
  if (YAML_Parser::get_if_present(ic_node, "scale", scale))
  {
    scaleSpecified = true;
    if (scale <= 0.0)
    {
      stk::RuntimeDoomedAdHoc() << "Scale for IC facets must be >= 0.\n";
    }
  }

  Vector3d scaleVec{scale, scale, scale};
  std::vector<std::string> scaleComponents = {"scaleX","scaleY","scaleZ"};
  for (int i=0; i<3; i++)
  {
    if (YAML_Parser::get_if_present(ic_node, scaleComponents[i], scaleVec[i]))
    {
      if (scaleSpecified)
        stk::RuntimeDoomedAdHoc() << "Cannot specify both scale and " << scaleComponents[i] << "\n";
      if (scaleVec[i] <= 0.0)
        stk::RuntimeDoomedAdHoc() << scaleComponents[i] << " for IC facets must be >= 0.\n";
    }
  }

  double sign = 1.0;
  if (YAML_Parser::get_if_present(ic_node, "sign", sign))
  {
    if (sign != -1.0 && sign == 1.0)
    {
      stk::RuntimeDoomedAdHoc() << "Sign for IC facets must be -1 or 1.\n";
    }
  }

  std::transform(facet_format.begin(), facet_format.end(), facet_format.begin(), ::toupper);
  Faceted_Surface_From_File * surface = nullptr;

  if (facet_format == "STL")
    surface = new STLSurface(surface_name, parent_timer, facet_filename, sign, scaleVec);
  else if (facet_format == "FAC")
    surface = new FACSurface(surface_name, parent_timer, facet_filename, sign, scaleVec);
  else if (facet_format == "PLY")
    surface = new PLYSurface(surface_name, parent_timer, facet_filename, sign, scaleVec);
  else if (facet_format == "EXO")
    surface = new EXOSurface(surface_name, parent_timer, facet_filename, sign, scaleVec);
  else
    stk::RuntimeDoomedAdHoc() << "Unrecognized facet format: " << facet_format;

  return surface;
}

MeshSurface *
parse_mesh_surface(const YAML::Node & ic_node, LevelSet & ls)
{
  std::string surface_name;
  ThrowRequire(YAML_Parser::get_if_present(ic_node, "mesh", surface_name));

  double sign = 1.0;
  if (YAML_Parser::get_if_present(ic_node, "sign", sign))
  {
    if (sign != -1.0 && sign == 1.0)
    {
      stk::RuntimeDoomedAdHoc() << "Sign for IC facets must be -1 or 1.\n";
    }
  }

  ThrowErrorMsgIf(!ls.aux_meta().has_part(surface_name), "Could not locate a surface named " << surface_name);
  const stk::mesh::Part & io_part = ls.aux_meta().get_part(surface_name);

  const stk::mesh::Field<double>* coords = reinterpret_cast<const stk::mesh::Field<double>*>(&LevelSet::get_current_coordinates(ls.meta()).field());
  ThrowRequire(nullptr != coords);

  const stk::mesh::Selector surface_selector = stk::mesh::Selector(io_part);
  return new MeshSurface(ls.meta(), *coords, surface_selector, sign);
}

void
parse_composition_method(const YAML::Node & ic_node, IC_Alg& ic_alg)
{
  // This is a little strange because composition looks like an IC, but really sets a flag
  std::string composition_method;
  ThrowRequire(YAML_Parser::get_if_present(ic_node, "composition_method", composition_method));

  std::transform(composition_method.begin(), composition_method.end(), composition_method.begin(), ::toupper);

  if (composition_method == "MINIMUM_SIGNED_DISTANCE")
    ic_alg.set_composition_method(Composite_Surface::MINIMUM_SIGNED_DISTANCE);
  else if (composition_method == "MAXIMUM_SIGNED_DISTANCE")
    ic_alg.set_composition_method(Composite_Surface::MAXIMUM_SIGNED_DISTANCE);
  else
    stk::RuntimeDoomedAdHoc() << "Unrecognized composition_method: " << composition_method;
}

void
parse_IC(const YAML::Node & ic_node, LevelSet &ls)
{
  if (ic_node.Type() == YAML::NodeType::Scalar)
  {
    // support random specified as Scalar (no :)
    std::string ic_type = ic_node.as<std::string>();
    if (ic_type == "random")
    {
      Random * random = new Random(0);
      ls.get_IC_alg().addSurface(random);
    }
    else if (ic_type == "analytic_isosurface")
    {
      Analytic_Isosurface * surf = new Analytic_Isosurface();
      ls.get_IC_alg().addSurface(surf);
    }
    else
    {
      stk::RuntimeDoomedAdHoc() << "Unrecognized Levelset IC type: " << YAML_Parser::info(ic_node);
    }
    return;
  }

  if ( YAML_Parser::get_null_if_present(ic_node, "sphere") )
  {
    ls.get_IC_alg().addSurface(parse_sphere(ic_node));
  }
  else if ( YAML_Parser::get_null_if_present(ic_node, "ellipsoid") )
  {
    ls.get_IC_alg().addSurface(parse_ellipsoid(ic_node));
  }
  else if ( YAML_Parser::get_null_if_present(ic_node, "plane") )
  {
    ls.get_IC_alg().addSurface(parse_plane(ic_node));
  }
  else if ( YAML_Parser::get_null_if_present(ic_node, "cylinder") )
  {
    ls.get_IC_alg().addSurface(parse_cylinder(ic_node));
  }
  else if ( YAML_Parser::get_null_if_present(ic_node, "facets") )
  {
    ls.get_IC_alg().addSurface(parse_facets(ic_node, ls.get_timer()));
  }
  else if ( YAML_Parser::get_scalar_if_present(ic_node, "mesh") )
  {
    ls.get_IC_alg().addSurface(parse_mesh_surface(ic_node, ls));
  }
  else if ( YAML_Parser::get_null_if_present(ic_node, "random") )
  {
    int seed = 0;
    YAML_Parser::get_if_present(ic_node, "seed", seed);
    Random * random = new Random(seed);
    ls.get_IC_alg().addSurface(random);
  }
  else if ( YAML_Parser::get_null_if_present(ic_node, "binder") )
  {
    ls.get_IC_alg().addCalculator(parse_binder(ic_node));
  }
  else if ( YAML_Parser::get_scalar_if_present(ic_node, "composition_method") )
  {
    parse_composition_method(ic_node, ls.get_IC_alg());
  }
  else
  {
    stk::RuntimeDoomedAdHoc() << "Unrecognized Levelset IC type: " << YAML_Parser::info(ic_node);
  }
}

}

void
IC_Parser::parse(const YAML::Node & node, LevelSet & ls)
{
  const YAML::Node ic_nodes = YAML_Parser::get_sequence_if_present(node, "initial_conditions");
  if ( ic_nodes )
  {
    for ( auto && ic_node : ic_nodes )
    {
      parse_IC(ic_node, ls);
    }
  }
}

} // namespace krino
